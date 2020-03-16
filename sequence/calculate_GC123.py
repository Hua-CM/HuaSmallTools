# -*- coding: utf-8 -*-
# @Time : 2020/3/16 15:12
# @Author : Zhongyi Hua
# @FileName: calculate_GC123.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

from Bio.SeqUtils import GC123
from Bio import SeqIO
import pandas as pd


def calculate123(seq_path):
    result_df = pd.DataFrame(columns=["seqid", "GC", "GC12", "GC1", "GC2", "GC3"])
    seq_list = SeqIO.parse(seq_path, "fasta")
    for record in seq_list:
        tmp_dict = dict()
        tmp_dict["seqid"] = record.id
        tmp_dict["GC"], tmp_dict["GC1"], tmp_dict["GC2"], tmp_dict["GC3"] = GC123(record.seq)
        tmp_dict["GC12"] = (tmp_dict["GC1"] + tmp_dict["GC2"])/2
        result_df = result_df.append(tmp_dict, sort=False, ignore_index=True)
    return result_df


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This is the script for calculate GC123")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<file path>  The cds sequence in fasta format')
    parser.add_argument('-t', '--out_file', required=True,
                        help='<file path>  A readable table')
    args = parser.parse_args()
    gc_df = calculate123(args.input_file)
    gc_df.to_csv(args.out_file, sep="\t", index=False)
