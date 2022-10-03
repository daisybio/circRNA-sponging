import argparse
from email.errors import HeaderParseError

from pyliftover import LiftOver
import os
from bs4 import BeautifulSoup as bs
from selenium import webdriver
import pandas
from pathlib import Path
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select, WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.firefox.options import Options
import threading
import logging
import time
import concurrent.futures

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--organism", help="Organism in three letter code (hsa for human)", required=True)
parser.add_argument("-gv", "--genome_version", help="Used genome version", required=True)
parser.add_argument("-d", "--data_loc", help="Location of data to be converted", required=True)
# optional
parser.add_argument("-ao", "--annotated_only", help="Keep annotated circRNAs only (true/false)", default=None)
parser.add_argument("-out", "--output", help="Output directory", default="./")
parser.add_argument("-s", "--separator", help="Separator of file", default="\t")
parser.add_argument("-chrC", "--chromosome", help="Column of chromosome", default=0)
parser.add_argument("-startC", "--start", help="Column of start of position", default=1)
parser.add_argument("-endC", "--end", help="Column of end of position", default=2)
parser.add_argument("-strandC", "--strand", help="Column of strand", default=3)
parser.add_argument("-nh", "--no_header", help="Specified file has no header", action='store_true')
parser.add_argument("-off", "--offline_access", default="None", help="Location of database offline data")
parser.add_argument("-split", "--splitter", default=1000, help="Maximum number of entries before splitting into search threads")

args = parser.parse_args()

# globally save column specifications
chr_c = args.chromosome
start_c = args.start
end_c = args.end
strand_c = args.strand
organism = args.organism
splitter = args.splitter
annotated_only = True if args.annotated_only == "true" else False

# set pipeline home
pipeline_home = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_home = os.path.join(pipeline_home, "data", "circBase")


class DbOrganism:
    def __init__(self, name, version, data):
        self.name = name
        self.version = version
        self.data = data

    # return db search format
    def get_db_name(self):
        return self.name + " (" + self.version + ")"


class Database:
    def __init__(self, name, test_url, search_url):
        self.name = name
        self.test_url = test_url
        self.search_url = search_url


db_genome_versions = {
    "hsa": DbOrganism("Human", "hg19", os.path.join(data_home, "hsa_hg19_circRNA.txt")),
    "mmu": DbOrganism("Mouse", "mm9", os.path.join(data_home, "mmu_mm9_circRNA.txt")),
    "cel": DbOrganism("C.elegans", "ce6", os.path.join(data_home, "cel_ce6_Memczak2013.txt")),
    "clc": DbOrganism("L.chalumnae", "latCha1", os.path.join(data_home, "lme_latCha1_Nitsche2013.txt"))
}

# species, original and converted genome
organism = db_genome_versions[organism]
converted_genome = organism.version
original_genome = str(args.genome_version)
offline_data = organism.data

dbs = {
    "circBase": Database("circBase", "http://www.circbase.org/", "http://www.circbase.org/cgi-bin/listsearch.cgi")
}


# generate key for specific circRNA
def key_gen(c, x, y, s):
    return str(c)+":"+str(x)+"-"+str(y)+"_"+str(s)


# get key information from key_gen as bed format
def key_to_bed(k):
    s1 = k.split(":")   # "chr1":,"start-end_strand"
    c = s1[0]   # chr
    s2 = s1[1].split("_")   # "start-end", "strand"
    strand = s2[1]  # strand
    lr = s2[0].split("-")   # "start", "end"
    s = lr[0]   # start
    e = lr[1]   # end
    return " ".join([c, s, e, strand])


# format for database
def db_format(c, x, y, s):
    return "\t".join([c, x, y, s])+"\n"


# generate string to search in circBase
def tsvData(data):
    x = []
    for key in data.keys():
        if key != "header":
            x.append(key_to_bed(key))
    return x


# READ FOUND CIRC RNA AND CONVERT GENOMIC POSITIONS
def convert(c, x, y, s, converter):
    first = converter.convert_coordinate(c, int(x), s)
    second = converter.convert_coordinate(c, int(y), s)
    if len(first) == 0 or len(second) == 0:
        return None, None
    return str(first[0][1]), str(second[0][1])


# read input data, convert each position, write tmp file for db, return whole converted data
def convert_and_write(o, c, d, s, nh):
    print("lifting coordinates")
    # create LiftOver for desired genome versions
    converter = LiftOver(o, c)
    data = {}
    with open(d, "r") as file:
        i = 0
        for line in file:
            # split line according to given separator
            split = line.rstrip("\n").split(s)
            chromosome = split[chr_c]
            start_pos = split[start_c]
            end_pos = split[end_c]
            strand = split[strand_c]
            # no header option
            if i == 0:
                if not nh:
                    i += 1
                    data["header"] = split
                    continue

            # convert position
            x_converted, y_converted = convert(chromosome, start_pos, end_pos, strand, converter)
            # skip if no conversion possible
            if x_converted is None:
                print("Conversion failed for line " + str(i) + ": " + str(s).join([chromosome, start_pos, end_pos, strand]))
                continue
            # save line to data with according key
            k = key_gen(chromosome, x_converted, y_converted, strand)
            data[k] = split
            # increase counter
            i += 1
    return data

# OFFLINE ACCESS ----------------------------------------------------------


# reads the offline database
def read_db(db_loc):
    d = {}
    with open(db_loc, "r") as db:
        c = 0
        for line in db:
            split = line.rstrip("\n").split("\t")
            # save header
            if c == 0:
                d["header"] = split
                c += 1
                continue
            d[key_gen(split[0], split[1], split[2], split[3])] = split      # key = position, value = all data
            c += 1
    return d


# write annotated output circ rna
def write_mapping_file(matched_dict, unannotated_dict, db_dict, output_loc, separator):
    # build expression header
    expr_header = matched_dict["header"][:7]
    expr_header.append("circBaseID")
    expr_header.extend(matched_dict["header"][7:])
    # build annotation header
    header = matched_dict["header"][:4]
    # add circRNA type
    header.append("type")
    header.append(matched_dict["header"][5])
    # remove headers from dicts
    matched_dict.pop("header")
    # get db header
    db_header = db_dict["header"]
    id_index = db_header.index("circRNA ID")
    db_dict.pop("header")

    # extend header with converted position
    header.append(converted_genome + "_converted_pos")
    # extend rest of db data without already present information
    header.extend(db_header[3:])
    print("[final] found " + str(len(matched_dict)) + " exact matching circRNAs to submitted query of size " + str(len(unannotated_dict)))
    annotation_file = os.path.join(output_loc, "circRNAs_annotated.tsv")
    annotated_expr = os.path.join(output_loc, "circRNA_counts_annotated.tsv")
    with open(annotation_file, "w") as output, open(annotated_expr, "w") as an_expr:
        if header is not None:
            output.write(str(separator).join(header) + "\n")
            an_expr.write(str(separator).join(expr_header) + "\n")
        # filter output and merge annotations
        for k, v in matched_dict.items():
            db_info = db_dict[k]
            circBaseID = str(db_info[id_index])
            # write annotated expression
            an_expr.write(separator.join(v[:7]) + separator + circBaseID + separator + separator.join(v[7:]) + "\n")
            data = v[:4]
            # show circRNA type
            data.append(v[6])
            # include ensgid
            data.append(v[5])
            # add converted position
            data.append(k)
            data.extend(db_info[id_index:])
            output.write(separator.join(data) + "\n")
        # add unannotated circRNAs, um des Flexes Willen
        if not annotated_only:
            for k, v in unannotated_dict.items():
                an_expr.write(separator.join(v[:7]) + separator + "None" + separator + separator.join(v[7:]) + "\n")


# LAUNCH OFFLINE MODE
def offline_access(converted_circ_data, output_loc, database_loc, separator):
    db_circ_rna = read_db(database_loc)
    matched_keys = set(converted_circ_data.keys()).intersection(set(db_circ_rna.keys()))
    m = {k: converted_circ_data[k] for k in matched_keys}
    unannotated = {k: converted_circ_data[k] for k in set(converted_circ_data.keys()).difference(set(db_circ_rna.keys()))}
    write_mapping_file(m, unannotated, db_circ_rna, output_loc, separator)


# ONLINE ACCESS -------------------------------------------------------------
def read_db_data_to_dict(header, data, out="raw.matches"):
    d = {}
    with open(out, "w") as o:
        for split in data:
            o.write("\t".join(split) + "\n")
            if len(split) < 3:
                continue
            key = str(split[1]) + "_" + str(split[2])
            d[key] = split
        d["header"] = header
    return d


def read_html(response):
    soup = bs(response, "html.parser")
    table = soup.find("table")
    header = [h.text for h in table.find_all("th")]
    data = [[d.text for d in h.find_all("td")] for h in table.find_all("tr")]
    return read_db_data_to_dict(header, data)


# check connection to circBase and if it responds to search queries
def is_connected(url, options):
    driver = webdriver.Firefox(options=options)
    driver.get(url=url)
    WebDriverWait(driver, 2).until(EC.presence_of_element_located((By.TAG_NAME, 'body')))
    textbox = driver.find_element(By.ID, "searchbox")   # enter dummy search
    textbox.click()
    textbox.send_keys("hsa_circ_0142668")
    driver.find_element(By.ID, "submit").click()        # query search
    response = bs(driver.page_source, "html.parser").find("title").text  # check response
    driver.quit()
    return "Error" not in response


# single search thread
def submit(tsv_data, url, options):
    driver = webdriver.Firefox(options=options)
    driver.get(url=url)
    WebDriverWait(driver, 2).until(EC.presence_of_element_located((By.TAG_NAME, 'body')))
    # select according organism
    organism_select = Select(driver.find_element(By.ID, "organism"))
    # select for current organism by DbOrganism class
    organism_select.select_by_value(organism.get_db_name())
    # upload tmp database file
    textbox = driver.find_element(By.ID, "querybox")
    # select it
    textbox.click()
    # fill data
    for circ in tsv_data:
        textbox.send_keys(circ)
        textbox.send_keys(Keys.ENTER)
    time.sleep(0.5)
    # submit form and retrieve data
    logging.info("Submitting data")
    driver.find_element(By.ID, "submit").click()
    delay = 300  # seconds
    try:
        WebDriverWait(driver, delay).until(EC.presence_of_element_located((By.ID, 'tablesorter')))
        logging.info("Results have appeared")
    except TimeoutException:
        logging.error("Timeout: circBase did not respond within " + str(delay) + " seconds")
        exit(1)
    source = driver.page_source
    driver.quit()
    # process response
    return read_html(source)


# make request using selenium
def online_access(converted_circ_data, output_loc, url, options):
    # if more than 2500 entries are supplied, thread execute database search with splitter max splits
    tsv_data = tsvData(converted_circ_data)
    if len(tsv_data) > splitter:
        print("Splitting search into " + str(len(range(0, len(converted_circ_data), splitter))) + " parts")
        with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
            logging.basicConfig(level=logging.INFO, format='%(relativeCreated)6d %(threadName)s %(message)s')
            results = [executor.submit(submit, d) for d in [tsv_data[i:i+splitter] for i in range(0, len(converted_circ_data), splitter)]]
            db_dict = {}
            i = 0
            for f in concurrent.futures.as_completed(results):
                # each threads database search result as dict
                r = f.result()
                print(str(i) + ": found ", str(len(r.keys())) + " overlapping circRNAs")
                if len(db_dict.keys()) == 0:
                    db_dict = r         # initiate first entry
                else:
                    db_dict.update(r)   # add next entries
                i += 1
    else:
        db_dict = submit(tsv_data, url, options)
    # extract only direct matches of genomic positions and no overlaps
    direct_matches = {x: converted_circ_data[x] for x in set(converted_circ_data.keys()).intersection(set(db_dict.keys()))}
    # extract all unannotated circRNAs
    unannotated_matches = {x: converted_circ_data[x] for x in set(converted_circ_data.keys()).difference(set(db_dict.keys()))}
    print("Writing output file")
    write_mapping_file(direct_matches, unannotated_matches, db_dict, output_loc, "\t")


def main():
    # location of found circRNA
    circ_rna_loc = args.data_loc
    # output loc
    out_loc = args.output
    # check for online or offline access
    db_data = args.offline_access
    # check separator option
    separator = args.separator
    # NO HEADER OPTION
    no_header = args.no_header

    # CIRC RNA CONVERSION, TMP FILE
    converted_data = convert_and_write(original_genome, converted_genome,
                                       d=circ_rna_loc, s=separator, nh=no_header)
    db = dbs["circBase"]  # set db url
    FIREFOX_OPTS = Options()  # set firefox options
    FIREFOX_OPTS.log.level = "trace"  # Debug
    FIREFOX_OPTS.headless = True
    # CHECK DATABASE ACCESS
    if not Path(db_data).is_file() and is_connected(db.test_url, FIREFOX_OPTS):
        print("Attempting database online access")
        online_access(converted_circ_data=converted_data,
                      output_loc=out_loc,
                      url=db.search_url, options=FIREFOX_OPTS)
    # OFFLINE ACCESS
    else:
        # GIVEN PATH
        if Path(db_data).is_file():
            f = Path(db_data)
            data = db_data
        # DEFAULT PATH FROM GIT
        else:
            f = Path(offline_data)
            data = offline_data

        print("Attempting circBase offline access using " + str(f))
        offline_access(converted_circ_data=converted_data,
                       output_loc=out_loc,
                       database_loc=data,
                       separator=separator)


if __name__ == "__main__":
    main()
