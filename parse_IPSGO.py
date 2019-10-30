# -*- coding: utf-8 -*-
# @Time : 2019/10/30 11:12
# @Author : Zhongyi Hua
# @FileName: parse_IPSGO.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com
import argparse
import pandas as pd


def get_gene2go(input_file):
    """
    parse InterProscan result file
    :param input_file: InterProscan result file
    :return: gene-go id mapping DataFrame. columns ["gene","GO term"]
    """
    pas = pd.read_csv(input_file, sep="\t", names=[i for i in range(1, 15)]).dropna()
    GO_list = pd.DataFrame(columns=["gene", "GO_num"])
    for _ in range(pas.__len__()):
        trans, GO_num = pas.iloc[_, 0], pas.iloc[_, 13]
        for i in GO_num.split("|"):
            tmp_dict = {"gene": trans, "GO_num": i}
            GO_list = GO_list.append(tmp_dict, ignore_index=True)
    return GO_list


# Combine trans2GO DataFrame with GO annotaion
def get_GOlist(gene2go, go_file):
    """
    Combine gene-go mapping and go-term mapping results.
    :param gene2go: DataFrame from get get_gene2go
    :param go_file: GO file path. This file should be parsed by parse_go_obobfile.py
    :return:
    """
    GO_df = pd.read_csv(go_file, sep="\t", names=["GO_num", "GO_description", "GO_level"])
    trans2go = pd.merge(GO_df, gene2go, on="GO_num")
    return trans2go


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This is the script for parsing GO database OBO file to a more \
                                                 user-friendly format")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<filepath>  The GO database file in OBO format')
    parser.add_argument('-g', '--go_file', required=True,
                        help='<filepath>  The GO database file(come from parse_go_obofile.py)')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<filepath>  The result path')
    args = parser.parse_args()
    GO_list = get_GOlist(get_gene2go(args.input_file), args.go_file)
    GO_list = GO_list[["gene", "GO_num", "GO_level", "GO_description"]]
    GO_list.to_csv(args.output_file, sep="\t", index=0)
