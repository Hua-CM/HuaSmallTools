# -*- coding: utf-8 -*-
# @Time : 2020/3/13 22:02
# @Author : Zhongyi Hua
# @FileName: trim_align.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

from Bio import AlignIO


def trim_align(alignment, mismatch):
    for col in range(alignment.get_alignment_length()):
        if alignment[:, col].count("-")/alignment.__len__() < mismatch:
            start_position = col
            break
    for col in range(alignment.get_alignment_length()-1, -1, -1):
        if alignment[:, col].count("-")/alignment.__len__() < mismatch:
            end_position = col
            break
    return alignment[:, start_position:end_position]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This is the script for simple trim MSA file")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<file_path> The MSA file in fasta format')
    parser.add_argument('-m', '--miss_match', default=0, type=float,
                        help='<float> The miss match value used for trimming')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<file_path> The result MSA file in format')
    args = parser.parse_args()
    tmp_align = AlignIO.read(args.input_file, "fasta")
    tmp_align = trim_align(tmp_align, args.miss_match)
    AlignIO.write(tmp_align, args.output_file, "fasta")
