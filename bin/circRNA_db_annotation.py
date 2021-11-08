import argparse

from pyliftover import LiftOver
import os
from bs4 import BeautifulSoup as bs
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.options import Options
import mygene

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--organism", help="Organism on three letter code (hsa for human)", required=True)
parser.add_argument("-gv", "--genome_version", help="Used genome version", required=True)
parser.add_argument("-d", "--data_loc", help="Location of data to be converted", required=True)
parser.add_argument("-out", "--output", help="Output file location", required=True)
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
    "clc": DbOrganism("L.chalumnae","latCha1")
}

# species, original and converted genome
organism = db_genome_versions[organism]
converted_genome = organism.version
original_genome = str(args.genome_version)

dbs = {
    "circBase": "http://www.circbase.org/cgi-bin/listsearch.cgi",
}

# columns to be annotated from circBase
annotated_columns = ['organism', 'circRNA_ID', 'genomic_length', 'spliced_length',
                     'samples', 'scores', 'repeats', 'annotation', 'best_transcript', 'gene_symbol', 'circRNA_study']


# generate key for specific circRNA
def key_gen(c, x, y, s):
    return str(c)+":"+str(x)+"-"+str(y)+"_"+str(s)


# format for database
def db_format(c, x, y, s):
    return "\t".join([c, x, y, s])+"\n"


# READ FOUND CIRC RNA AND CONVERT GENOMIC POSITIONS
def convert(c, x, y, s, converter):
    first = converter.convert_coordinate(c, int(x), s)
    second = converter.convert_coordinate(c, int(y), s)
    if len(first) == 0:
        return None, None
    return str(first[0][1]), str(second[0][1])


# read input data, convert each position, write tmp file for db, return whole converted data
def convert_and_write(o, c, d, s, tmp_out, nh):
    # create LiftOver for desired genome versions
    converter = LiftOver(o, c)
    data = {}
    with open(d, "r") as file, open(tmp_out, "w") as out_file:
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
                print("Conversion failed for line " + str(i) + ": " + str(s).join(split))
                continue
            # save line to data with according key
            data[key_gen(chromosome, x_converted, y_converted, strand)] = split
            out_file.write(db_format(chromosome, x_converted, y_converted, strand))
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


# get ensemble id for gene symbol
def get_ens_id(query, gene_symbol):
    try:
        ens_id = query.T[str(gene_symbol)].iat[2, 0]
    except (KeyError, IndexError, TypeError) as e:
        ens_id = "None"
    return ens_id


# write annotated output circ rna
def write_mapping_file(matched_dict, db_dict, output_loc, separator):
    # build gene symbol converter
    mg = mygene.MyGeneInfo()
    header = matched_dict["header"][:4]
    matched_dict.pop("header")
    db_header = db_dict["header"]
    db_dict.pop("header")

    # extend header with converted position
    header.append(converted_genome + "_converted_pos")
    # extend rest of db data without already present information
    header.extend(db_header[3:])
    # ensembl gene id
    header.append("Ensembl_gene_ID")
    with open(output_loc, "w") as output:
        if header is not None:
            output.write(str(separator).join(header) + "\n")
        # annotate matches and save in mapping file: gen.pos of or. genome -> circ_base_id
        file_data = []
        symbols = []
        # filter output and merge annotations
        for k, v in matched_dict.items():
            db_info = db_dict[k]
            data = v[:4]
            data.append(k)
            data.extend(db_info[3:])
            file_data.append(data)
            symbols.append(db_info[11])
        # add ensembl ids
        query = mg.querymany(symbols, scopes='symbol', fields='ensembl.gene', as_dataframe=True)
        for line, symbol in zip(file_data, symbols):
            ens_id = get_ens_id(query, symbol)
            line.append(str(ens_id))
            output.write(str(separator).join(line) + "\n")


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
        if len(split) == 0 or split[12] == "NA":
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


# handle html response from database
def process_db_response(response, converted_circ_data, output_loc):
    db_dict = read_html(response)
    direct_matches = {x: converted_circ_data[x] for x in set(converted_circ_data.keys()).intersection(set(db_dict.keys()))}
    write_mapping_file(direct_matches, db_dict, output_loc, "\t")


# make request using selenium
def online_access(upload_file, converted_circ_data, output_loc):
    # get circBase url
    url = dbs["circBase"]
    # create Chrome driver and navigate to circBase list search
    options = Options()
    options.add_argument("--headless")
    firefox_path = os.path.join(pipeline_home, "assets/geckodriver")
    driver = webdriver.Firefox(firefox_path, options=options)
    driver.get(url=url)
    # select according organism
    organism_select = Select(driver.find_element(By.ID, "organism"))
    # select for current organism by DbOrganism class
    organism_select.select_by_value(organism.get_db_name())
    # upload tmp database file
    driver.find_element(By.ID, "queryfile").send_keys(upload_file)
    # submit form and retrieve data
    driver.find_element(By.ID, "submit").click()
    # process response
    process_db_response(driver.page_source, converted_circ_data, output_loc)


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

    # TMP FILE LOCATION
    tmp_db_file = os.path.join(os.path.dirname(os.path.realpath(circ_rna_loc)), "tmp_converted_for_db.tsv")

    # CIRC RNA CONVERSION, TMP FILE
    converted_data = convert_and_write(original_genome, converted_genome, d=circ_rna_loc,
                                       tmp_out=tmp_db_file, s=separator, nh=no_header)

    # DATABASE ACCESS
    if db_data == "None":
        online_access(upload_file=tmp_db_file, converted_circ_data=converted_data, output_loc=out_loc)
    else:
        offline_access(converted_circ_data=converted_data, output_loc=out_loc, database_loc=db_data, separator=separator)


if __name__ == "__main__":
    main()
