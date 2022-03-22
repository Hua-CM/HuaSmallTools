# -*- coding: utf-8 -*-
# @Time : 2022/2/7 23:50
# @Author : Zhongyi Hua
# @FileName: gene_rename_PASA.py
# @Usage: tidy gene name after PASA update
# @Note:
# @E-mail: njbxhzy@hotmail.com

import argparse
import re
from BCBio import GFF
from collections import defaultdict


def sort_seq(seqrecord_lst, _chr, _scaffold):
    def num_sort(seq_record):
        return list(map(int, re.findall(r'\d+', seq_record.id)))[0]

    _chr_lst = []
    _scaffold_lst = []
    for _ in seqrecord_lst:
        if re.search(rf'^{_chr}\d+', _.id):
            _chr_lst.append(_)
        elif re.search(rf'^{_scaffold}\d+', _.id):
            _scaffold_lst.append(_)
    _chr_lst.sort(key=num_sort)
    _scaffold_lst.sort(key=num_sort)
    return _chr_lst, _scaffold_lst


def rename(_gff, _out, _chr, _scaffold, _prefix):
    ele_dict = {'exon': 'exon', 'CDS': 'cds', 'three_prime_UTR': 'utr3p', 'five_prime_UTR': 'utr5p'}

    _gff = [_ for _ in GFF.parse(_gff)]
    _chr_lst, _scaffold_lst = sort_seq(_gff, _chr, _scaffold)

    def per_seq(_seqrecord, _gene_count, _chr_num, type):
        gene_pre = '%s%02dG%05d' if type == _chr else '%s%02dS%05d'
        _seqrecord.features.sort(key=lambda x: x.location.start)
        for _gene in _seqrecord.features:
            _gene_count += 1
            mrna_count = 0
            _gene.id = gene_pre % (_prefix, _chr_num, _gene_count)
            _gene.qualifiers['ID'] = [_gene.id]
            _gene.qualifiers.pop('Name', None)
            for _mRNA in _gene.sub_features:
                mrna_count += 1
                _mRNA.qualifiers.pop('Name', None)
                _mRNA.id = _gene.id + '.%d' % mrna_count
                _mRNA.qualifiers['Parent'] = [_gene.id]
                _mRNA.qualifiers['ID'] = [_mRNA.id]
                ele_count = defaultdict(int)
                for _ele in _mRNA.sub_features:
                    ele_count[_ele.type] += 1
                    _ele.id = '%s.%s%d' % (_mRNA.id, ele_dict.get(_ele.type), ele_count[_ele.type])
                    _ele.qualifiers['Parent'] = [_mRNA.id]
                    _ele.qualifiers['ID'] = [_ele.id]
                _mRNA.sub_features.sort(key=lambda x: x.location.start)
        return _gene_count

    gene_count = 0
    for _idx, _seq in enumerate(_chr_lst):
        gene_count = per_seq(_seq, gene_count, _idx + 1, _chr)
    for _idx, _seq in enumerate(_scaffold_lst):
        _seq_num = int(re.findall(r'\d+', _seq.id)[0])
        gene_count = per_seq(_seq, gene_count, _seq_num, _scaffold)
    GFF.write(_chr_lst + _scaffold_lst, open(_out, 'w'), include_fasta=False)


def main():
    parser = argparse.ArgumentParser(description="Convert gene id in gff file from EVM into readble format")
    parser.add_argument('gff', help="gff file")
    parser.add_argument('out', help="output file")
    parser.add_argument('-c', '--chr', default='chr', help="chromosome prefix")
    parser.add_argument('-s', '--scaffold', default='scaffold', help="scaffold prefix")
    parser.add_argument('-p', '--prefix', default='Gene', help="gene prefix")
    args = parser.parse_args()

    rename(args.gff, args.out, args.chr, args.scaffold, args.prefix)


if __name__ == '__main__':
    main()
