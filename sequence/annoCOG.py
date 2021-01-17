# -*- coding: utf-8 -*-
# @Time : 2021/1/17 19:34
# @Author : Zhongyi Hua
# @FileName: rhizobium_codon.py
# @Usage:
# @Note:
# @E-mail: njbxhzy@hotmail.com

import argparse
import os
from tempfile import mktemp
import pandas as pd
from Bio.Blast.Applications import NcbiblastpCommandline


def annoCOG(_seq_path, _result_path, _db_path, _gi2cogid, _coguid2cog, _tmp_path, _threads=4):
    cog_blast = NcbiblastpCommandline(
            query=_seq_path,
            db=_db_path,
            dust='no',
            outfmt='\"6 qacc sacc evalue sgi\"',
            num_threads=_threads,
            evalue='1e-3',
            out=_tmp_path,
            max_hsps=1,
            max_target_seqs=5)
    cog_blast()
    # parse result
    gi2COGID = pd.read_csv(_gi2cogid, sep=",", usecols=[0, 6], names=["gi", "COGID"])
    COGID2COG = pd.read_csv(_coguid2cog, sep="\t", usecols=[0, 1],
                            names=["COGID", "COG"])
    cog_blast_result = pd.read_csv("COG.blast.tsv", sep="\t", names=["query", "subject", "evalue", "gi"])
    COG_blast_mapped = pd.merge(pd.merge(cog_blast_result, gi2COGID, on="gi", how="left"),
                                COGID2COG,
                                on="COGID",
                                how="left")
    COG_blast_count_filter = pd.merge(
        COG_blast_mapped.groupby("query")["COGID"].value_counts().to_frame("COG_count").reset_index(),
        COGID2COG,
        on="COGID",
        how="left")
    COG_blast_count_filter = COG_blast_count_filter[COG_blast_count_filter['COG_count'] > 2]
    COG_blast_count_filter["COG2"] = COG_blast_count_filter["COG"].apply(lambda x: " ".join(list(x)))
    COG_blast_2 = COG_blast_count_filter.drop('COG2', axis=1).\
        join(COG_blast_count_filter["COG2"].str.split(' ', expand=True).stack().reset_index(level=1, drop=True).rename(
            'COG2'))
    COG_blast_2["species"] = COG_blast_2["query"].apply(lambda x: x[0:4].replace("_", ""))
    COG_blast_2 = COG_blast_2.merge((1 / COG_blast_2["query"].value_counts()).to_frame("weight"),
                                    left_on='query',
                                    right_index=True)
    COG_blast_2 = COG_blast_2[["query", "COG_count", "COG2", "species", "weight"]]
    COG_blast_2.to_csv(_result_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="This is the script for annotating COG")
    parser.add_argument('-i', '--input_fasta', required=True,
                        help='<filepath>  The protein fasta file')
    parser.add_argument('-d', '--database', required=True,
                        help='<filepath>  The COG BLAST database')
    parser.add_argument('-f', '--file1', required=True,
                        help='<filepath>  The cog2003-2014.csv')
    parser.add_argument('-F', '--file2', required=True,
                        help='<filepath>  cognames2003-2014.tab')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='<int> Number of threads')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<filepath>  The output path')
    args = parser.parse_args()
    tmp_path = mktemp()
    annoCOG(args.input_fasta, args.output_file, args.database, args.file1, args.file2, tmp_path, args.threads)
    os.remove(tmp_path)


if __name__ == '__main__':
    main()
