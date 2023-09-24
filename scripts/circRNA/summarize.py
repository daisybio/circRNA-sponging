#!/usr/bin/env python

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Summarize circRNA counts')

parser.add_argument('--input', '-i', type=argparse.FileType('r'), 
                    required=True, nargs='+', help='Input files')
parser.add_argument('--output', '-o', type=argparse.FileType('w'), 
                    required=True, help='Output file')

args = parser.parse_args()

grouping_cols = ["chr", "start", "stop", "strand", "gene_symbol", "type"]

def process_sample_name(sample_file):
    return sample_file.name.split('/')[-1].split('.')[0], sample_file

def read_df(sample_name: str, sample_file: str):
    df = pd.read_csv(sample_file, sep='\t', 
                     usecols=[0,1,2,5,12,13,14], 
                     header=None, 
                     names=["chr", "start", "stop", "strand", sample_name, "type", "gene_symbol"])

    # Drop lines wich chromosomes containing '_'
    df = df[~df['chr'].str.contains('_')]

    # Group by chr, start, stop, strand, gene_symbol and type, take max of counts
    # TODO: Check if summation is better
    df = df.groupby(grouping_cols, as_index=False).max()

    return df

samples = [process_sample_name(sample) for sample in args.input]

dataframes = { sample_name: read_df(sample_name, sample_file) 
              for sample_name, sample_file in samples }

merged = pd.concat(dataframes).groupby(grouping_cols, as_index=False).sum()

merged.to_csv(args.output, sep='\t', index=False)