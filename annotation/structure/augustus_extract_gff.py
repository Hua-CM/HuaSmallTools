# -*- coding: utf-8 -*-
# @Time    : 2021/9/15 19:07
# @Author  : Zhongyi Hua
# @File    : augustus_extract_gff.py
# @Usage   : remove similar proteins in Augustus training dataset
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import argparse
import copy
from BCBio import GFF


def getArgs():
    parser = argparse.ArgumentParser(description="Sort chromosome by their length. !!Must use fasta without wrap, \
                                                  otherwise it will be very very slow.")
    parser.add_argument('lst', help="gene list")
    parser.add_argument('gff', help="gff file path")
    parser.add_argument('output', help="output path")
    args = parser.parse_args()
    return args


def main(args):
    with open(args.lst) as f_in:
        keeplist = f_in.read().strip().split('\n')
    gff3_all = [_ for _ in GFF.parse(args.gff)]
    gff3_out = []
    for _scaffold in gff3_all:
        _scaffold_out = copy.copy(_scaffold)
        _scaffold_out.features = []
        for _gene in _scaffold.features:
            for _subfeature in _gene.sub_features:
                if _subfeature.type == 'mRNA' and _subfeature.id in keeplist:
                    _scaffold_out.features.append(_gene)
        gff3_out.append(_scaffold_out)
    GFF.write(gff3_out, open(args.output, 'w'))


if __name__ == '__main__':
    main(getArgs())
