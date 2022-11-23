# -*- coding: utf-8 -*-
# @Time : 2022/2/23 0:58
# @Author : Zhongyi Hua
# @FileName: parse_privateSNP.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import pandas as pd
import numpy as np
import argparse
import vcf


def find_marker1(group_df: np.array, other_df: np.array):
    snp_tmp = []
    # round1 calculate and record
    for _idx in range(len(group_df)):
        _var_ = [_ for _ in group_df[_idx] if not _ == -1]
        _mean = np.mean(_var_)
        failure = (other_df[_idx, :] == _mean) * 1
        if all(_var_ == _mean):
            snp_tmp.append(
                (_idx,
                 np.sum(failure),
                 np.array(np.where(failure)[0].tolist())
                 )
            )
    snp_tmp.sort(key=lambda x: x[1])
    return snp_tmp


def find_marker2(group_df: np.array, other_df: np.array, putative: list, order, markers):
    order = order if len(markers) == 1 else 0
    next_rest_snp = np.array([_[0] for _ in putative[order+1:]])
    next_rest_sample = putative[order][2]
    print(next_rest_sample)
    _tmp = []
    for _idx in next_rest_snp:
        geno = [_ for _ in group_df[_idx] if not _ == -1][0]
        failure = (other_df[_idx, next_rest_sample] == geno) * 1
        _tmp.append(
            (_idx,
             np.sum(failure),
             next_rest_sample[np.where(failure)[0].tolist()]
             )
        )
    _tmp.sort(key=lambda x: x[1])
    markers.append(_tmp[0][0])
    if (_tmp[0][2].size == 0) or (len(markers) >= 10):
        return markers
    else:
        return find_marker2(group_df, other_df, _tmp, 0, markers)


def find_marker2(group_df: np.array, other_df: np.array, snp_list):
    """
    :param snp_list: return from find_marker1
    :return:
    """
    snp_tmp2 = []
    batch_list = [snp_list[i:i + 100] for i in range(0, len(snp_list), 100)]
    snp_index = []
    # set 100 as as a batch step
    """
    # strategy one
    # the order may be very important ?
    other_index_set = set(_ for _ in range(other_df.shape[1]))
    for _snp_tmp_lst in batch_list:
        for _var_idx, _num, _fail in _snp_tmp_lst:
            snp_index.append(_var_idx)
            other_index_set.difference_update(_fail)
            if not other_index_set:
                return snp_index
    """
    """
    # strategy two: tweak strategy one to check pair-by-pair
    # very very slow, not a good idea
    for _var_idx, _num, _fail in snp_tmp:
        for _var_idx2, _num2, _fail2 in snp_tmp:
            genotype_group = set(list(zip(group_df[_var_idx], group_df[_var_idx2])))
            genotype_other = list(zip(other_df[_var_idx], other_df[_var_idx2]))
            if not genotype_group & set(genotype_other):
                return _var_idx, _var_idx2
            else:
                snp_tmp2.append(
                    (
                        (_var_idx, _var_idx2),
                        sum(genotype_other == genotype_group),
                        np.where(genotype_other == genotype_group)
                    )
                )
    for _var_idx, _num, _fail in snp_tmp2:
        pass
    """
    # just another strategy (search specific SNP that could solve rest samples)
    for _var_idx, _num, _fail in snp_tmp:
        # group_geno = [_ for _ in group_df[_var_idx] if not _ == -1][0]
        for _var_idx2, _num2, _fail2 in snp_tmp:
            if not _fail & _fail2:
                return _var_idx, _var_idx2


class VCF2Matrix:
    # too slow
    def __init__(self, vcf_path):
        self.vcf = vcf.Reader(filename=vcf_path)
        self.genotype = {'0/0':0, '0/1':1, '1/1':2}
        self.matrix = []
        self.sample = []
        self.snp = []

    def generate_ele(self):
        self.sample = self.vcf.samples
        for record in self.vcf:
            if record.is_snp and len(record.ALT) == 1:
                self.snp.append([record.CHROM, record.POS, record.REF, record.ALT])
                self.matrix.append(np.array(self.genotype.get(record.samples[0].gt_nums, -1)))
        return self.matrix, self.sample, self.snp


def main(args):
    snp_matrix = np.loadtxt(args.input, dtype=int)
    sample_info = pd.read_table(args.input + '.indv', header=None, names=['sample'])
    snp_info = pd.read_table(args.input + '.pos', header=None, sep='\s+', names=['seqid', 'postion', 'ref', 'alt'], engine='python')
    group_info = pd.read_table(args.group, header=None, names=['sample', 'group'])
    target = args.target # target group
    # remove snp has snp with in 20 bp
    del_lst = [cur for cur in range(len(snp_info) - 1) if snp_info.iloc[cur, 1] + 20 > snp_info.iloc[cur + 1, 1]]
    del_lst += list(np.array(del_lst) + 1)
    del_lst = list(set(del_lst))
    select_lst = list(set(range(len(snp_info))) - set(del_lst))
    snp_matrix = snp_matrix[select_lst]
    snp_info = snp_info.iloc[select_lst, :]
    # split dataset
    sample_lst = group_info.loc[group_info['group'] == target, 'sample'].to_list()
    sample_index = sample_info[sample_info['sample'].isin(sample_lst)].index.to_list()
    other_index = sample_info[~(sample_info['sample'].isin(sample_lst))].index.to_list()
    group_df = snp_matrix[:, sample_index]
    other_df = snp_matrix[:, other_index]
    #rest_snp = np.array([_ for _ in range(group_df.shape[0])])
    #rest_sample = np.array([_ for _ in range(other_df.shape[1])])
    snp_tmp = find_marker1(group_df, other_df)
    marker_lst = []
    for marker_num in range(10):
        marker_pair = []
        marker_pair.append(snp_tmp[marker_num][0])
        snp_tmp2 = find_marker1(group_df, other_df[:, snp_tmp[marker_num][2]])
        for marker_num2 in range(10):
            marker_pair.append(
                find_marker2(group_df, other_df, snp_tmp, marker_num2, [snp_tmp2[marker_num2][2]])
            )
        marker_lst.append(marker_pair)
    #snp_info.iloc[marker1,:]

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Filter ')
    parser.add_argument('-i', '--input', required=True,
                        help='<file path>  snp_matrix from vcftools. The *.012 file'
                             '(*.012.indv and *.012.pos must be in the same directory)')
    parser.add_argument('-g', '--group', required=True,
                        help='<file path> group file, two columns: sample/group')
    parser.add_argument('-t', '--target', required=True,
                        help='<character> the group you want to identify')
    parser.add_argument('-o', '--output', required=True,
                        help='<file path>  The output path')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main(parseArgs())