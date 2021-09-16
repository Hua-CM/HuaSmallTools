# -*- coding: utf-8 -*-
# @Time    : 2021/9/15 22:17
# @Author  : Zhongyi Hua
# @File    : rm_similarity_protein.py
# @Usage   : remove similar proteins in Augustus training dataset
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import argparse
from BCBio import GFF
import copy


def getArgs():
    parser = argparse.ArgumentParser(description="Sort chromosome by their length. !!Must use fasta without wrap, \
                                                  otherwise it will be very very slow.")
    parser.add_argument('blast', help="blast.out")
    parser.add_argument('gff', help="gff file path")
    parser.add_argument('output', help="output path")
    args = parser.parse_args()
    return args


def get_rm_list(blast_out):
    with open(blast_out) as f_in:
        pairs = f_in.read().strip().split('\n')
    keep_set = set()
    rm_set = set()
    for pair in pairs:
        pair = pair.split('\t')
        bool01 = pair[0] in keep_set
        bool02 = pair[0] in rm_set
        bool11 = pair[1] in keep_set
        bool12 = pair[1] in rm_set
        if bool01 + bool02 + bool11 + bool12 == 0:
            keep_set.add(pair[0])
            rm_set.add(pair[1])
        elif bool01 + bool11 == 2:
            keep_set.remove(pair[0])
            rm_set.add(pair[0])
        elif bool02 + bool12 == 2:
            continue
        elif bool01 + bool12 == 2:
            continue
        elif bool02 + bool11 == 2:
            continue
        elif bool01 + bool02 == 0 and bool12:
            keep_set.add(pair[0])
        elif bool01 + bool02 == 0 and bool12:
            rm_set.add(pair[0])
        elif bool11 + bool12 == 0 and bool01:
            rm_set.add(pair[1])
        elif bool11 + bool12 == 0 and bool02:
            keep_set.add(pair[1])
    return list(rm_set)


def main(args):
    rm_list = get_rm_list(args.blast)
    gff3_all = [_ for _ in GFF.parse(args.gff)]
    gff3_out = []
    for _scaffold in gff3_all:
        _scaffold_out = copy.copy(_scaffold)
        _scaffold_out.features = []
        for _gene in _scaffold.features:
            rm_mark = False
            for _subfeature in _gene.sub_features:
                if _subfeature.type == 'mRNA' and _subfeature.id in rm_list:
                    rm_mark = True
            if rm_mark:
                continue
            else:
                _scaffold_out.features.append(_gene)
        gff3_out.append(_scaffold_out)
    GFF.write(gff3_out, open(args.output, 'w'))


if __name__ == '__main__':
    main(getArgs())
