import argparse
from pyliftover import LiftOver
import os
from bs4 import BeautifulSoup as bs
from selenium import webdriver
import pandas as pd
from pathlib import Path
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select, WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.firefox.options import Options
import threading
import logging
import time
import concurrent.futures

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--organism", help="Organism in three letter code (hsa for human)", required=True)
parser.add_argument("-gv", "--genome_version", help="Used genome version", required=True)
parser.add_argument("-d", "--data_loc", help="Location of data to be converted", required=True)
parser.add_argument("-m", "--meta", help="Path to meta file of circRNA expression", required=True)
# optional
parser.add_argument("-out", "--output", help="Output directory", default="./")
parser.add_argument("-s", "--separator", help="Separator of file", default="\t")
parser.add_argument("-chrC", "--chromosome", help="Column of chromosome", default=0)
parser.add_argument("-startC", "--start", help="Column of start of position", default=1)
parser.add_argument("-endC", "--end", help="Column of end of position", default=2)
parser.add_argument("-strandC", "--strand", help="Column of strand", default=3)
parser.add_argument("-ann", "--annotated_only", help="Only keep annotated circRNAs", action="store_true")
parser.add_argument("-gecko", "--gecko_binary", help="Location of geckodriver executable", default=None)
parser.add_argument("-off", "--offline_access", default="None", help="Location of database offline data")
parser.add_argument("-split", "--splitter", default=1000,
                    help="Maximum number of entries before splitting into search threads")

args = parser.parse_args()

# globally save column specifications
chr_c = args.chromosome
start_c = args.start
end_c = args.end
strand_c = args.strand
organism = args.organism
splitter = args.splitter
# filtering option
annotated_only = args.annotated_only

# read meta file
meta = pd.read_csv(args.meta, delimiter="\t", )
# set sample names
samples = meta["sample"]

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

# create converter
CONVERTER = LiftOver(original_genome, converted_genome)

dbs = {
    "circBase": Database("circBase", "http://www.circbase.org/", "http://www.circbase.org/cgi-bin/listsearch.cgi")
}


# generate key for specific circRNA
def key_gen(c, x, y, s):
    return str(c) + ":" + str(x) + "-" + str(y) + "_" + str(s)


# format for database
def db_format(c, x, y, s):
    return "\t".join([c, x, y, s]) + "\n"


# READ FOUND CIRC RNA AND CONVERT GENOMIC POSITIONS
def convert(row):
    c = row["chr"]
    x = row["start"]
    y = row["stop"]
    s = row["strand"]
    first = CONVERTER.convert_coordinate(c, int(x), s)
    second = CONVERTER.convert_coordinate(c, int(y), s)
    if len(first) == 0 or len(second) == 0:
        print("Uplift failed for circRNA: ", ":".join(row.values[0:4].astype(str)))
        return None, None
    return key_gen(c, str(first[0][1]), str(second[0][1]), s)


def set_key(row):
    return key_gen(row["chr"], row["start"], row["stop"], row["strand"])


# read input data, convert each position, write tmp file for db, return whole converted data
def lifted_data(d, s):
    print("lifting coordinates")
    # read input circRNA expression
    data = pd.read_csv(d, sep=s)
    # convert positions and change row names to converted samples
    data["converted_key"] = data.apply(convert, axis=1)
    return data.set_index("converted_key")


def annotate_expression(converted_circ_data, db_result):
    # extract direct matches of genomic positions and no overlaps
    direct_matches = set(converted_circ_data.index).intersection(set(db_result.index))
    direct_matches_df = converted_circ_data.loc[direct_matches]
    # annotate direct matches with circBase ID
    direct_matches_df.insert(loc=6, column="circBaseID", value=db_result["circRNA ID"])
    print(f"Found a total of {len(direct_matches)} exact matching circRNAs in database")

    # extract all unannotated circRNAs
    unannotated_matches = set(converted_circ_data.index) - direct_matches
    no_match_df = converted_circ_data.loc[unannotated_matches]
    # set annotation to "None"
    no_match_df.insert(loc=6, column="circBaseID", value=["None"] * len(no_match_df))

    # only include annotated circRNAs
    if annotated_only:
        print("Excluding not annotated circRNAs")
        circRNA_expression = direct_matches_df
    # concat both parts
    else:
        print("Including not annotated circRNAs")
        circRNA_expression = pd.concat([direct_matches_df, no_match_df])
    # reset row names to original genome versions, replace converted index
    circRNA_expression["ID"] = circRNA_expression.apply(set_key, axis=1)
    return circRNA_expression.set_index("ID")


# OFFLINE ACCESS ----------------------------------------------------------
# reads the offline database
def read_db(db_loc):
    # read circBase database file
    data = pd.read_csv(db_loc, sep="\t")
    # set index to key_gen
    data["key"] = data.apply(set_key, axis=1)
    return data.set_index("key")


# LAUNCH OFFLINE MODE
def offline_access(database_loc):
    return read_db(database_loc)


# ONLINE ACCESS -------------------------------------------------------------
def read_html(response):
    data = pd.read_html(response)[0]

    def key(row):
        return row["position (genome browser link)"] + "_" + row["strand"]

    data["key"] = data.apply(key, axis=1)
    return data.set_index("key")


# check connection to circBase and if it responds to search queries
def is_connected(url, options):
    print("Testing url: " + str(url) + "...")
    driver = webdriver.Firefox(options=options)
    driver.get(url=url)
    WebDriverWait(driver, 2).until(EC.presence_of_element_located((By.TAG_NAME, 'body')))
    textbox = driver.find_element(By.ID, "searchbox")  # enter dummy search
    textbox.click()
    textbox.send_keys("hsa_circ_0142668")
    driver.find_element(By.ID, "submit").click()  # query search
    print("Entering dummy search: hsa_circ_0142668")
    WebDriverWait(driver, 2).until(EC.presence_of_element_located((By.TAG_NAME, 'body')))
    response = bs(driver.page_source, "html.parser").find("title").text  # check response
    driver.quit()
    return "Error" not in response


# single search thread
def submit(bed, url, options):
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
    for coord in bed:
        textbox.send_keys(coord)
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
def online_access(converted_circ_data, url, options):
    print("Building query...")
    # if more than 2500 entries are supplied, thread execute database search with splitter max splits
    bed = list(converted_circ_data[["chr", "start", "stop", "strand"]].astype(str).apply(" ".join, axis=1))
    if len(bed) > splitter:
        print("Splitting search into " + str(len(range(0, len(converted_circ_data), splitter))) + " parts")
        with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
            logging.basicConfig(level=logging.INFO, format='%(relativeCreated)6d %(threadName)s %(message)s')
            results = [executor.submit(submit, d, url, options)
                       for d in [bed[i:i + splitter]
                                 for i in range(0, len(converted_circ_data), splitter)
                                 ]
                       ]
            db_result = None
            i = 0
            for thread in concurrent.futures.as_completed(results):
                # each threads database search result as dict
                df = thread.result()
                print(str(i) + ": found ", str(len(df)) + " overlapping circRNAs")
                if db_result is None:
                    db_result = df  # initiate first data frame
                else:
                    db_result = pd.concat([db_result, df])  # add other frames
                i += 1
    else:
        db_result = submit(bed, url, options)
    return db_result


def main():
    # location of found circRNA
    circ_rna_loc = args.data_loc
    # output loc
    out_loc = args.output
    # check for online or offline access
    db_data = args.offline_access
    # check separator option
    separator = args.separator

    # CIRC RNA CONVERSION
    converted_data = lifted_data(d=circ_rna_loc, s=separator)
    db = dbs["circBase"]  # set db url
    FIREFOX_OPTS = Options()  # set firefox options
    FIREFOX_OPTS.log.level = "trace"  # Debug
    FIREFOX_OPTS.headless = True
    # add geckodriver binary location
    if args.gecko_binary is not None:
        from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
        FIREFOX_OPTS.binary = FirefoxBinary(args.gecko_binary)

    # CHECK DATABASE ACCESS
    if not Path(db_data).is_file() and is_connected(db.test_url, FIREFOX_OPTS):
        print("Connected to database")
        database_result = online_access(converted_circ_data=converted_data,
                                        url=db.search_url, options=FIREFOX_OPTS)
    # OFFLINE ACCESS
    else:
        # GIVEN PATH
        if Path(db_data).is_file():
            database_loc = db_data
        # DEFAULT PATH FROM GIT
        else:
            database_loc = offline_data

        print("Attempting circBase offline access using " + str(Path(database_loc)))
        database_result = offline_access(database_loc=database_loc)

    # remove duplicate matches
    database_result = database_result[~database_result.index.duplicated(keep='first')]
    # annotate with database ID
    circRNA_expression = annotate_expression(converted_circ_data=converted_data, db_result=database_result)
    # remove index label
    circRNA_expression.index.name = None
    # save annotated expression
    circRNA_expression.to_csv(os.path.join(out_loc, "circRNA_counts_annotated.tsv"), sep="\t")
    # save database output
    database_result.to_csv(os.path.join(out_loc, "circRNAs_annotated.tsv"), sep="\t")


if __name__ == "__main__":
    main()
