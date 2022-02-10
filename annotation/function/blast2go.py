# -*- coding: utf-8 -*-
# @Time : 2022/1/31 14:46
# @Author : Zhongyi Hua
# @FileName: blast2go.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import argparse
import gzip
import json
import subprocess as sp
import os
from collections import defaultdict


def mmseqs_wrapper(_query, _db, _threads, _tmp, _top):
    # using 'tmp' as tmp
    sp.run(['mmseqs', 'createdb', _query, os.path.join(_tmp, 'queryDB')])
    sp.run(['mmseqs', 'search', os.path.join(_tmp, 'queryDB'), _db, os.path.join(_tmp, 'resultDB'), _tmp,
            '-s', '7',
            '--threads', str(_threads),
            '-e', '1e-3',
            '--max-seqs', str(_top)])
    sp.run(['mmseqs', 'convertalis',
            os.path.join(_tmp, 'queryDB'),
            _db,
            os.path.join(_tmp, 'resultDB'),
            os.path.join(_tmp, 'tmp.out'),
            '--format-mode', '0', '--format-output', 'query,target'])


def parse_result(_mmeseq, _map, _out):
    res_dict = defaultdict(set)
    res_list = []
    with gzip.open(_map, 'r') as f_in:
        map_dict = json.loads(f_in.read())
    with open(_mmeseq, 'r') as f_in:
        for line in f_in.readlines():
            _query, _subject = line.strip().split()
            if go := map_dict.get(_subject):
                res_dict[_query] = res_dict[_query].union(set(go))
    for _gene, _go_set in res_dict.items():
        res_list += [_gene + '\t' + _go for _go in list(_go_set)]
    with open(_out, 'w') as f_out:
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
    tmp_dir = os.path.join(os.getcwd(), 'tmp')
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    mmseqs_wrapper(args.input, args.db, args.threads, tmp_dir, args.top)
    parse_result(os.path.join(tmp_dir, 'tmp.out'), args.map, args.output)


if __name__ == '__main__':
    main(parse_args())
