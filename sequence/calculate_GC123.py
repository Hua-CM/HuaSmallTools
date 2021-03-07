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
from sequence.calculate_motif import motif123
import argparse
from os.path import split as splt


def calculate123(seq_path):
    seq_list = SeqIO.parse(seq_path, "fasta")
    dict_list = []
    for record in seq_list:
        tmp_dict = dict()
        tmp_dict["GC"], tmp_dict["GC1"], tmp_dict["GC2"], tmp_dict["GC3"] = GC123(record.seq)/100
        tmp_dict['slen'] = len(record.seq)
        dict_list.append(tmp_dict)
    result_df = pd.DataFrame(dict_list)
    result_dict = {'seqid': splt(seq_path)[1],
                   'GC': sum(result_df.GC * result_df.slen) / result_df.slen.sum()/100,
                   'GC1': sum(result_df.GC1 * result_df.slen) / result_df.slen.sum()/100,
                   'GC2': sum(result_df.GC2 * result_df.slen) / result_df.slen.sum()/100,
                   'GC3': sum(result_df.GC3 * result_df.slen) / result_df.slen.sum()/100}
    return result_dict


def calculateATCG3(seq_path):
    result_df = pd.DataFrame(columns=["seqid", "A3", "T3", "G3", "C3"])
    seq_list = SeqIO.parse(seq_path, "fasta")
    for record in seq_list:
        tmp_dict = dict()
        tmp_dict["seqid"] = record.id
        tmp_dict["A3"] = motif123(record.seq, ["A"])[3]
        tmp_dict["T3"] = motif123(record.seq, ["T"])[3]
        tmp_dict["C3"] = motif123(record.seq, ["C"])[3]
        tmp_dict["G3"] = motif123(record.seq, ["G"])[3]
        result_df = result_df.append(tmp_dict, sort=False, ignore_index=True)
    result_df.iloc[:, 1:] = result_df.iloc[:, 1:] / 100
    return result_df


def parse_args():
    parser = argparse.ArgumentParser(description="This is the script for calculate GC123")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<file path>  A file containing the path of the cds sequence in fasta format')
    parser.add_argument('-a', '--ATCG', default=False, action='store_true',
                        help='calculate A3, T3, G3, C3 as well')
    parser.add_argument('-t', '--out_file', required=True,
                        help='<file path>  A readable table')
    args = parser.parse_args()
    return args


def main(args):
    result_list = []
    with open(args.input_file) as f:
        file_list = [_.strip() for _ in f.read().split('\n')]
    for file in file_list:
        try:
            result_list.append(calculate123(file))
        except:
            continue
    GC_DF = pd.DataFrame(result_list)
    if args.ATCG:
        result_list2 = []
        for file in file_list:
            result_list2.append(calculateATCG3(file))
        GC_DF2 = pd.DataFrame(result_list2)
        GC_DF = pd.merge(GC_DF, GC_DF2, on='seqid')
    GC_DF.to_csv(args.out_file, sep="\t", index=False)


if __name__ == '__main__':
    main(parse_args())
