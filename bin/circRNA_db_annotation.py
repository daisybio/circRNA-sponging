import argparse

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
from selenium.webdriver.chrome.options import Options
import threading
import logging
import time
import concurrent.futures

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--organism", help="Organism in three letter code (hsa for human)", required=True)
parser.add_argument("-gv", "--genome_version", help="Used genome version", required=True)
parser.add_argument("-d", "--data_loc", help="Location of data to be converted", required=True)
# optional
parser.add_argument("-out", "--output", help="Output directory", default="./")
parser.add_argument("-s", "--separator", help="Separator of file", default="\t")
parser.add_argument("-chrC", "--chromosome", help="Column of chromosome", default=0)
parser.add_argument("-startC", "--start", help="Column of start of position", default=1)
parser.add_argument("-endC", "--end", help="Column of end of position", default=2)
parser.add_argument("-strandC", "--strand", help="Column of strand", default=3)
parser.add_argument("-nh", "--no_header", help="Specified file has no header", action='store_true')
parser.add_argument("-off", "--offline_access", default="None", help="Location of database offline data")

args = parser.parse_args()

# globally save column specifications
chr_c = args.chromosome
start_c = args.start
end_c = args.end
strand_c = args.strand
organism = args.organism

# set pipeline home
pipeline_home = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class DbOrganism:
    def __init__(self, name, version):
        self.name = name
        self.version = version

    # return db search format
    def get_db_name(self):
        return self.name + " (" + self.version + ")"


db_genome_versions = {
    "hsa": DbOrganism("Human", "hg19"),
    "mmu": DbOrganism("Mouse", "mm9"),
    "cel": DbOrganism("C.elegans", "ce6"),
    "clc": DbOrganism("L.chalumnae", "latCha1")
}

# species, original and converted genome
organism = db_genome_versions[organism]
converted_genome = organism.version
original_genome = str(args.genome_version)

dbs = {
    "circBase": "http://www.circbase.org/cgi-bin/listsearch.cgi",
}


# generate key for specific circRNA
def key_gen(c, x, y, s):
    return str(c)+":"+str(x)+"-"+str(y)+"_"+str(s)


# format for database
def db_format(c, x, y, s):
    return "\t".join([c, x, y, s])+"\n"


# generate string to search in circBase
def tsvData(data):
    x = []
    for entry in data.values():
        x.append(" ".join([entry[chr_c], entry[start_c], entry[end_c], entry[strand_c]]))
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
            pos = str(split[0]) + "_" + str(split[1])
            d[pos] = split[2:]
            c += 1
    return d


# write annotated output circ rna
def write_mapping_file(matched_dict, db_dict, output_loc, separator):
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
    db_dict.pop("header")

    # extend header with converted position
    header.append(converted_genome + "_converted_pos")
    # extend rest of db data without already present information
    header.extend(db_header[3:])
    print(len(matched_dict))
    annotation_file = os.path.join(output_loc, "circRNAs_annotated.tsv")
    annotated_expr = os.path.join(output_loc, "circRNA_counts_annotated.tsv")
    with open(annotation_file, "w") as output, open(annotated_expr, "w") as an_expr:
        if header is not None:
            output.write(str(separator).join(header) + "\n")
            an_expr.write(str(separator).join(expr_header) + "\n")
        # filter output and merge annotations
        for k, v in matched_dict.items():
            db_info = db_dict[k]
            circBaseID = str(db_info[3])
            # write annotated expression
            an_expr.write(separator.join(v[:7]) + separator + circBaseID + separator + separator.join(v[7:]) + "\n")
            data = v[:4]
            # show circRNA type
            data.append(v[6])
            # include ensgid
            data.append(v[5])
            # add converted position
            data.append(k)
            data.extend(db_info[3:])
            output.write(separator.join(data) + "\n")


# LAUNCH OFFLINE MODE
def offline_access(converted_circ_data, output_loc, database_loc, separator):
    db_circ_rna = read_db(database_loc)
    matched_keys = set(converted_circ_data.keys()).intersection(set(db_circ_rna.keys()))
    m = {k: converted_circ_data[k] for k in matched_keys}
    write_mapping_file(m, db_circ_rna, output_loc, separator)


# ONLINE ACCESS -------------------------------------------------------------
def read_db_data_to_dict(header, data):
    d = {}
    for split in data:
        if len(split) < 3:
            continue
        key = str(split[1]) + "_" + str(split[2])
        d[key] = split
    d["header"] = header
    return d


def read_html(response):
    print("processing results")
    soup = bs(response, "html.parser")
    table = soup.find("table")
    header = [h.text for h in table.find_all("th")]
    data = [[d.text for d in h.find_all("td")] for h in table.find_all("tr")]
    return read_db_data_to_dict(header, data)


def submit(tsv_data):
    # get circBase url
    url = dbs["circBase"]
    # create driver and navigate to circBase list search
    options = Options()
    # set headless
    # options.headless = True
    # firefox_path = os.path.join(pipeline_home, "assets/geckodriver")
    chrome_path = "/Users/leonschwartz/Desktop/Bioinformatik/local_data/software/chromedriver"
    # driver = webdriver.Firefox(executable_path=firefox_path, options=options)
    driver = webdriver.Chrome(executable_path=chrome_path, options=options)
    driver.get(url=url)
    WebDriverWait(driver, 2).until(EC.presence_of_element_located((By.TAG_NAME, 'body')))
    # select according organism
    organism_select = Select(driver.find_element(By.ID, "organism"))
    # select for current organism by DbOrganism class
    organism_select.select_by_value(organism.get_db_name())
    # upload tmp database file
    textbox = driver.find_element(By.ID, "querybox")
    WebDriverWait(driver, 2).until(EC.element_to_be_clickable(textbox))
    # select it
    textbox.click()
    # fill data
    for circ in tsv_data:
        textbox.send_keys(circ)
        textbox.send_keys(Keys.ENTER)
    time.sleep(1)
    # submit form and retrieve data
    print("submit")
    logging.info("Submitting data")
    driver.find_element(By.ID, "submit").click()
    delay = 300  # seconds
    try:
        WebDriverWait(driver, delay).until(EC.presence_of_element_located((By.ID, 'tablesorter')))
        print("got it")
        logging.info("Results have appeared")
    except TimeoutException:
        logging.error("Timeout: cirBase did not respond within " + str(delay) + " seconds")
        exit(1)
    # process response
    return read_html(driver.page_source)


# make request using selenium
def online_access(converted_circ_data, output_loc, splitter):
    # if more than 2500 entries are supplied, thread execute database search with 2500 max splits
    tsv_data = tsvData(converted_circ_data)[1:10]
    if len(tsv_data) > splitter:
        with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
            logging.basicConfig(level=logging.INFO, format='%(relativeCreated)6d %(threadName)s %(message)s')
            results = [executor.submit(submit, driver, "\n".join(d)) for d in [tsv_data[i:i+splitter] for i in range(0, len(converted_circ_data), splitter)]]
            db_dict = {}
            for f in concurrent.futures.as_completed(results):
                # each threads database search result as dict
                r = f.result()
                # first entry
                if len(db_dict.keys()) == 0:
                    db_dict = r
                else:
                    db_dict.update(r)   # add next entries
    else:
        db_dict = submit(tsv_data)
    direct_matches = {x: converted_circ_data[x] for x in set(converted_circ_data.keys()).intersection(set(db_dict.keys()))}
    print("Writing output file")
    write_mapping_file(direct_matches, db_dict, output_loc, "\t")
    # close driver
    driver.quit()


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
    converted_data = convert_and_write(original_genome, converted_genome, d=circ_rna_loc,
                                        s=separator, nh=no_header)

    # DATABASE ACCESS
    if not Path(db_data).is_file():
        print("Attempting circBase online access")
        online_access(converted_circ_data=converted_data,
                      output_loc=out_loc,
                      splitter=2500)
    else:
        print("Using circBase offline access")
        offline_access(converted_circ_data=converted_data,
                       output_loc=out_loc,
                       database_loc=db_data,
                       separator=separator)


if __name__ == "__main__":
    main()
