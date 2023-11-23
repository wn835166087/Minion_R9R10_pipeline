#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO

def extract_reads(input_fasta, output_fasta, read_ids):
    read_id_set = set(read_ids)
    with open(output_fasta, 'w') as output_handle:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            if record.id in read_id_set:
                SeqIO.write(record, output_handle, 'fasta')


def main(dasotpt_path):
    source_files = sorted(Path(dasotpt_path).glob('*_DASTool_contig2bin.tsv'))

    for file in source_files:
        sample_id = file.name.split("_")[0:2]
        df = pd.read_csv(file, sep='\t', names=["contigs", "bin_name"])
        dfs = dict(tuple(df.groupby('bin_name')))
        for bin_name, bin_df in dfs.items():
            contig_ls = bin_df['contigs'].tolist()
            polished_contigs = dasotpt_path + sample_id[0] + "_" + sample_id[1] + "_" + 'consensus.fasta'
            bin_fasta = dasotpt_path + sample_id[0] + "_" + sample_id[1] + "_" + bin_name + '.fasta'
            extract_reads(polished_contigs, bin_fasta, contig_ls)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract reads from DASTool output.")
    parser.add_argument("dasotpt_path", type=str, help="Path to DASTool output directory.")
    args = parser.parse_args()

    main(args.dasotpt_path)
