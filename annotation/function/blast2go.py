# -*- coding: utf-8 -*-
# @Time : 2022/1/31 14:46
# @Author : Zhongyi Hua
# @FileName: blast2go.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

from Bio.Blast.Applications import NcbiblastpCommandline
import argparse
import gzip
import json
from collections import defaultdict


def blast_wrapper(_query, _db, _out, _threads, _top):
    cline = NcbiblastpCommandline(
        query=_query,
        db=_db,
        evalue=0.001,
        max_hsps=1,
        num_threads=_threads,
        max_target_seqs=_top,
        out=_out,
        outfmt="6 qseqid sseqid")
    cline()


def parse_result(_blast, _map):
    res_dict = defaultdict()
    res_list = []
    with gzip.open(_map, 'r') as f_in:
        map_dict = json.loads(f_in.read())
    with open(_blast, 'r') as f_in:
        for line in f_in.readlines():
            _query, _suject = line.split()
            if go := map_dict.get(_suject):
                res_dict.setdefault(_query, set()).add(set(go))
    for _gene, _go_set in res_dict.items():
        res_list += [_gene + '\t' + _go for _go in list(_go_set)]
    with open(_blast, 'w') as f_out:
        f_out.write('\n'.join(res_list))


def parse_args():
    parser = argparse.ArgumentParser(description='Blast2GO')
    parser.add_argument('-i', '--input', required=True,
                        help='<File path> The input protein sequences in fasta format')
    parser.add_argument('-d', '--db', required=True,
                        help='<File path> The BLAST proteome database')
    parser.add_argument('-m', '--map', required=True,
                        help='<File path> The idmapping file in json.gz format')
    parser.add_argument('-@', '--threads', default=4, type=int,
                        help='<int> Threads number Default: 4')
    parser.add_argument('-t', '--top', default=20, type=int,
                        help='<int> Using X top hit for annotation. Default: 20')
    parser.add_argument('-o', '--output', default='GO.annotation',
                        help='<File path> The output results')
    args = parser.parse_args()
    return args


def main(args):
    blast_wrapper(args.input, args.db, args.output, args.threads, args.top)
    parse_result(args.output, args.map)


if __name__ == '__main__':
    main(parse_args())
