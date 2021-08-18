# -*- coding: utf-8 -*-
# @Time    : 2021/8/17 21:13
# @Author  : Zhongyi Hua
# @File    : gene_rename.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import argparse
import pandas as pd
import re
from collections import defaultdict
from os import remove


class GFF:
    def __init__(self, gff, _chr, _scaf, _prefix, _out):
        self.chr = _chr
        self.scaf = _scaf
        self.prefix = _prefix
        self.gff = self.__parse_gff__(gff)
        self.__dialect__ = self.__record_name__()
        self.__rename__()
        self.__rewrite__(_out)

    def __parse_gff__(self, _gff):
        gff_new = pd.read_table(_gff,
                                header=None,
                                names=("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"))
        gff_new[['seqid', 'type']] = gff_new[['seqid', 'type']].astype('category')
        # feature order
        list_custom = ['gene', 'mRNA', 'five_prime_UTR', 'exon', 'CDS', 'three_prime_UTR']
        # get correct chromosome order
        list_seq = list(gff_new['seqid'].unique())
        ## chromosome first
        chromo = [int(re.sub(self.chr, '', _)) for _ in list_seq if re.search('^' + self.chr, _)]
        chromo.sort()
        chromo = [self.chr + str(_) for _ in chromo]
        ## scaffold second
        scaffold = [int(re.sub(self.scaf, '', _)) for _ in list_seq if re.search('^' + self.scaf, _)]
        scaffold.sort()
        scaffold = [self.scaf + str(_) for _ in scaffold]
        list_seq = chromo + scaffold
        # sort
        gff_new['seqid'].cat.reorder_categories(list_seq, inplace=True)
        gff_new['type'].cat.reorder_categories(list_custom, inplace=True)
        gff_new.sort_values(['seqid', 'start', 'type'], inplace=True)
        gff_new.to_csv(_gff+'.tmp', sep='\t', index=False, header=False)
        # parse
        with(open(_gff+'.tmp')) as f:
            gff_dialect = f.readlines()
        result_list = []
        for _record in gff_dialect:
            try:
                _record = defaultdict(str,
                                      zip(("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"),
                                           _record.strip().split('\t')))
                _record['attr'] = {_.split('=')[0]: _.split('=')[1] for _ in _record['attr'].strip(';').split(';')}
                result_list.append(_record)
            except:
                continue
        remove(_gff+'.tmp')
        return result_list

    def __record_name__(self):
        gene_count = 1
        gene_dict = defaultdict(str)
        mrna_dict = defaultdict(str)
        exon_dict = defaultdict(str)
        cds_dict = defaultdict(str)
        utr3_dict = defaultdict(str)
        utr5_dict = defaultdict(str)
        for _record in self.gff:
            if _record['type'] == 'gene':
                if _record['seqid'].startswith(self.chr):
                    gene_dict[_record['attr']['ID']] = \
                        self.prefix + _record['seqid'].replace(self.chr, '').zfill(2) + 'G' + str(gene_count).zfill(5)
                elif _record['seqid'].startswith(self.scaf):
                    gene_dict[_record['attr']['ID']] = \
                        self.prefix + _record['seqid'].replace(self.scaf, '').zfill(2) + 'S' + str(gene_count).zfill(5)
                gene_count += 1
            elif _record['type'] == 'mRNA':
                mrna_dict[_record['attr']['ID']] = _record['attr']['Parent']
            elif _record['type'] == 'exon':
                exon_dict[_record['attr']['ID']] = _record['attr']['Parent']
            elif _record['type'] == 'CDS':
                cds_dict[_record['attr']['ID']] = _record['attr']['Parent']
            elif _record['type'] == 'five_prime_UTR':
                utr5_dict[_record['attr']['ID']] = _record['attr']['Parent']
            elif _record['type'] == 'three_prime_UTR':
                utr3_dict[_record['attr']['ID']] = _record['attr']['Parent']
        return gene_dict, mrna_dict, exon_dict, cds_dict, utr3_dict, utr5_dict

    def __rename__(self):
        count_dict = defaultdict(int)
        count_dict2 = defaultdict(int)
        count_dict3 = defaultdict(int)
        count_dict4 = defaultdict(int)
        parent_dict = defaultdict(str)
        for _key in self.__dialect__[1].keys():
            gene_id = self.__dialect__[0][self.__dialect__[1][_key]]
            count_dict[gene_id] += 1
            self.__dialect__[1][_key] = gene_id + '.' + str(count_dict[gene_id])
            parent_dict[_key] = gene_id
        for _key in self.__dialect__[2].keys():
            mrna_id = self.__dialect__[1][self.__dialect__[2][_key]]
            count_dict[mrna_id] += 1
            self.__dialect__[2][_key] = mrna_id + '.exon.' + str(count_dict[mrna_id])
            parent_dict[_key] = mrna_id
        for _key in self.__dialect__[3].keys():
            mrna_id = self.__dialect__[1][self.__dialect__[3][_key]]
            count_dict2[mrna_id] += 1
            self.__dialect__[3][_key] = mrna_id + '.cds.' + str(count_dict2[mrna_id])
            parent_dict[_key] = mrna_id
        for _key in self.__dialect__[4].keys():
            mrna_id = self.__dialect__[1][self.__dialect__[4][_key]]
            count_dict2[mrna_id] += 1
            self.__dialect__[4][_key] = mrna_id + 'utr3p1' + str(count_dict3[mrna_id])
            parent_dict[_key] = mrna_id
        for _key in self.__dialect__[5].keys():
            mrna_id = self.__dialect__[1][self.__dialect__[5][_key]]
            count_dict2[mrna_id] += 1
            self.__dialect__[5][_key] = mrna_id + 'utr5p' + str(count_dict4[mrna_id])
            parent_dict[_key] = mrna_id
        self.__dialect__ = ({_key: _value for dic in self.__dialect__ for _key, _value in dic.items()}, parent_dict)

    def __rewrite__(self, _out):
        rewrite_list = []
        for _record in self.gff:
            if _record['attr'].get('Parent'):
                _record['attr']['Parent'] = self.__dialect__[1][_record['attr']['ID']]
            _record['attr']['ID'] = self.__dialect__[0][_record['attr']['ID']]
            _record['attr'] = ';'.join([_key + '=' + _value for _key, _value in _record['attr'].items()])
            rewrite_list.append('\t'.join(_record.values()))
        with open(_out, 'w') as f_out:
            f_out.write('\n'.join(rewrite_list))
            f_out.write('\n')


def getArgs():
    parser = argparse.ArgumentParser(description="Convert gene id in gff file from EVM into readble format")
    parser.add_argument('gff', help="gff file")
    parser.add_argument('-p', '--prefix', type=str, required=True, help="gff file name")
    parser.add_argument('-c', '--chromosome', dest="chr", type=str, default='chr',
                        help="chromosome sequence prefix. Default:chr")
    parser.add_argument('-s', '--scaffold', type=str, default='scaffold',
                        help="chromosome sequence prefix. Default:scaffold")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    Args = getArgs()
    GFF(Args.gff, Args.chr, Args.scaffold, Args.prefix, 'rename.gff')
