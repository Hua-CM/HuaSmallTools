# -*- coding: utf-8 -*-
# @Time    : 2022/7/11 10:19
# @Author  : Zhongyi Hua
# @File    : statistic.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

from BCBio import GFF
from collections import defaultdict
from numpy import median, mean


def per_gene(gene):
    """
    :param gene: Bio.SeqFeature.SeqFeature
    :return:
    """
    mRNA_length_lst = []
    exon_count_lst = []
    exon_length_lst = []
    intron_count_lst = []
    intron_length_lst = []
    protein_length_lst = []  # based on the sum of CDS
    for mRNA in gene.sub_features:
        exon_count = 0
        _exon_length_lst = []
        _intron_length_lst = []
        CDS_length = 0
        exon_end = 0  # for intron
        mRNA.sub_features.sort(key=lambda x: x.location.start)
        for ele in mRNA.sub_features:
            if ele.type == 'exon':
                _exon_length_lst.append(len(ele))
                _intron_length_lst.append(ele.location.start - exon_end)
                exon_end = ele.location.end
                exon_count += 1
            elif ele.type == 'CDS':
                CDS_length += len(ele)
        # remove the first
        _intron_length_lst = _intron_length_lst[1:]
        intron_count = exon_count - 1
        # output
        mRNA_length_lst.append(len(mRNA))
        protein_length_lst.append(CDS_length/3)
        exon_length_lst += _exon_length_lst
        exon_count_lst.append(exon_count)
        intron_count_lst.append(intron_count)
        intron_length_lst += _intron_length_lst
    return [len(gene)], mRNA_length_lst, exon_length_lst, exon_count_lst, intron_length_lst, intron_count_lst, protein_length_lst


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        description='This is the script for counting gene, transcript, exon, intron, and their length')
    parser.add_argument('-g', '--gff', required=True,
                        help='<file_path>  gff file path')
    parser.add_argument('-o', '--out', required=True,
                        help='<file_path>  output file path')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    out_dict = defaultdict(list)

    for scaffold in GFF.parse(args.gff):
        for gene in scaffold.features:
            l1, l2, l3, l4, l5, l6, l7 = per_gene(gene)
            out_dict['gene_length'] += l1
            out_dict['mRNA_length'] += l2
            out_dict['exon_length'] += l3
            out_dict['exon_count'] += l4
            out_dict['intron_length'] += l5
            out_dict['intron_count'] += l6
            out_dict['protein_length'] += l7
    gene_mean, gene_median = mean(out_dict['gene_length']), median(out_dict['gene_length'])
    mRNA_mean, mRNA_median = mean(out_dict['mRNA_length']), median(out_dict['mRNA_length'])
    exon_mean, exon_median = mean(out_dict['exon_length']), median(out_dict['exon_length'])
    exon_ct_mean, exon_ct_median = mean(out_dict['exon_count']), median(out_dict['exon_count'])
    intron_mean, intron_median = mean(out_dict['intron_length']), median(out_dict['intron_length'])
    protein_mean, protein_median = mean(out_dict['protein_length']), median(out_dict['protein_length'])
    print('length (bp) of:\taverage\tmedian')
    print(f'gene\t{gene_mean:.2f}\t{gene_median:.2f}')
    print(f'transcript\t{mRNA_mean:.2f}\t{mRNA_median:.2f}')
    print(f'exon\t{exon_mean:.2f}\t{exon_median:.2f}')
    print(f'intron\t{intron_mean:.2f}\t{intron_median:.2f}')
    print(f'protein\t{protein_mean:.2f}\t{protein_median:.2f}')
    print(f'exon per gene\t{exon_ct_mean:.2f}\t{exon_ct_median:.2f}')
    print(f'total gene:\t{len(out_dict["gene_length"])}')


if __name__ == '__main__':
    main()
