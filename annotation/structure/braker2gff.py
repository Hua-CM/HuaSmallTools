# -*- coding: utf-8 -*-
# @Time    : 2022/3/4 15:05
# @Author  : Zhongyi Hua
# @File    : braker2gff.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import pandas as pd
from collections import defaultdict
from copy import deepcopy
from random import randint
import argparse


def main(args):
    gtf = pd.read_table(args.input,
                        header=None,
                        names=("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"))

    res_lst = []
    gene_num = 0

    for _, _row in gtf.iterrows():
        if _row['source'] == 'AUGUSTUS' and _row["type"] in {'CDS', 'exon', 'transcript', 'gene'}:
            tmp_dict = defaultdict()
            for _head in ["seqid", "source", "type", "start", "end", "score", "strand", "phase"]:
                tmp_dict[_head] = _row[_head]
            if _row["type"] == 'gene':
                tmp_dict['attr'] = 'ID=' + _row['attr']
            elif _row["type"] == 'transcript':
                tmp_dict['type'] = 'mRNA'
                tmp_dict['attr'] = 'ID={};Parent={}'.format(_row['attr'], _row['attr'].split('.')[0])
            else:
                _tmp_dict = {_.strip().split(' ')[0]: _.strip().split(' ')[1].strip('"') for _ in _row["attr"].split(';')[:-1]}
                tmp_dict['attr'] = 'ID={}_{};Parent={}'.\
                    format(_tmp_dict['transcript_id'].split('_')[-1] + _row["type"], randint(1, 100), _tmp_dict['transcript_id'].split('_')[-1])
            res_lst.append(tmp_dict)

    gff = pd.DataFrame(res_lst)
    gff.to_csv(args.output, sep='\t', index=False, header=False)


def parseArgs():
    parser = argparse.ArgumentParser(description='Design Primer using Primer3 for AS-PCR')
    parser.add_argument('-i', '--input', required=True,
                        help='<file path The input braker.gtf>')
    parser.add_argument('-o', '--output', required=True,
                        help='<file path>  The output path')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main(parseArgs())
