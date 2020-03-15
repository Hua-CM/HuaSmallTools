# -*- coding: utf-8 -*-
# @Time : 2020/3/12 17:40
# @Author : Zhongyi Hua
# @FileName: calculate_psy_chemi.py
# @Usage:
# @Note:
# @E-mail: njbxhzy@hotmail.com

import pandas as pd
from Bio.SeqUtils import ProtParam
from Bio import SeqIO


def calculate_property(seq_path):
    seq_fasta = SeqIO.parse(seq_path, "fasta")
    result_primary_feature = pd.DataFrame(columns=["SeqID",
                                                   "molecular_weight",
                                                   "instability_index",
                                                   "GRAVY",
                                                   "theoretical_pI"])
    func_dict = {
        "molecular_weight": ProtParam.ProteinAnalysis.molecular_weight,
        "instability_index": ProtParam.ProteinAnalysis.instability_index,
        "GRAVY": ProtParam.ProteinAnalysis.gravy,
        "theoretical_pI": ProtParam.ProteinAnalysis.isoelectric_point}
    for seq in seq_fasta:
        protein_seq = str(seq.seq).strip("*")
        protein_result = ProtParam.ProteinAnalysis(protein_seq)
        tmp_dict = {"SeqID": seq.id}
        for key, Prot_func in func_dict.items():
            try:
                tmp_dict[key] = Prot_func(protein_result)
            except BaseException:
                tmp_dict[key] = "NA"
        result_primary_feature = result_primary_feature.append(
            tmp_dict, ignore_index=True)
    return result_primary_feature


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This is the script for batch calculating protein physical-chemistry \
                                                  properties")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<filepath> The fasta file')
    parser.add_argument('-o', '--outpur_file', required=True,
                        help='<filepath> The result table')
    args = parser.parse_args()
    result_table = calculate_property(args.input_file)
    result_table.to_csv(args.output_file, sep="\t", index=None)
