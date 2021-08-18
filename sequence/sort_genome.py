# -*- coding: utf-8 -*-
# @Time    : 2021/8/17 14:10
# @Author  : Zhongyi Hua
# @File    : sort_genome.py
# @Usage   : sort chromosome by their lengths
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import pandas as pd
import argparse
import string
import random
from collections import OrderedDict, defaultdict


def get_fasta_len(_fasta, chr_num: int):
    with open(_fasta, 'r') as _f_in:
        _fasta = list(line.rstrip() for line in _f_in)
    fasta_dict = defaultdict(str)
    for _line in _fasta:
        if _line.startswith('>'):
            scaffold_id = _line.split()[0].replace('>', '')
            continue
        fasta_dict[scaffold_id] += _line
    len_dict = {_key: len(_value) for _key, _value in fasta_dict.items()}
    ordered_dict = OrderedDict(zip(OrderedDict(sorted(len_dict.items(), key=lambda x: x[1], reverse=True)).keys(),
                                   ['chr' + str(_) for _ in range(1, chr_num+1)] + ['scaffold' + str(_) for _ in range(1, len(len_dict)-chr_num+1)])
                               )
    return fasta_dict, ordered_dict


def id_generate(size=6, chars=string.ascii_uppercase):
    return ''.join(random.choice(chars) for _ in range(size))


def sort_genome(args):
    fa_dict, order_dict = get_fasta_len(args.fasta, args.chr)
    change_dict1 = defaultdict(str)
    change_dict2 = defaultdict(str)
    for raw_, ordered_ in order_dict.items():
        tmp_id = id_generate()
        change_dict1[raw_] = tmp_id
        change_dict2[tmp_id] = ordered_
    for old_key, new_key in change_dict1.items():
        fa_dict[new_key] = fa_dict.pop(old_key)
    for old_key, new_key in change_dict2.items():
        fa_dict[new_key] = fa_dict.pop(old_key)
    # gff
    pd_gff = pd.read_table(args.gff,
                           comment='#',
                           names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    pd_gff = pd_gff[pd_gff['seqid'].isin(list(change_dict1.keys()))]
    pd_gff['seqid'].replace(change_dict1, inplace=True)
    pd_gff['seqid'].replace(change_dict2, inplace=True)
    return fa_dict, pd_gff


def getArgs():
    parser = argparse.ArgumentParser(description="Sort chromosome by their length. !!Must use fasta without wrap, \
                                                  otherwise it will be very very slow.")
    parser.add_argument('fasta',
                        help="fasta file name")
    parser.add_argument('gff',
                        help="gff file name")
    parser.add_argument('-p', '--prefix', type=str, default="output",
                        help="prefix for output")
    parser.add_argument('-c', '--chr', type=int, default=0,
                        help="chromosome number. Default: 0 (scaffold level genome)")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    Args = getArgs()
    fasta, gff = sort_genome(Args)
    gff.to_csv(Args.prefix + '.gff', sep='\t', index=False, header=False)
    with open(Args.prefix + '.fasta', 'w') as f:
        f.write('\n'.join(['>' + _key + '\n' + _value for _key, _value in fasta.items()]))
        f.write('\n')
