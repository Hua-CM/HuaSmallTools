# -*- coding: utf-8 -*-
# @Time : 2019/10/25 10:46
# @Author : Zhongyi Hua
# @FileName: parse_eggNOG.py
# @Usage: """non-module organisms are always lack of KEGG and GO information. A convenient method is use eggNOG server
#            to annotate. This script is used for parse the eggNOG result to a user friendly format which could be
#            directly applied in clusterProfile R package"""
# @E-mail: njbxhzy@hotmail.com
import pandas as pd
import numpy as np
import argparse
import requests
import re
import os
from pathlib import Path

headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:65.0) Gecko/20100101 Firefox/65.0"}


def get_html(url):
    resp = requests.get(url, headers=headers)
    return resp.text


def own_mapdb(org_list: list):
    """
    Pick map num based on user prefer
    :param org_list:organism name list (KEGG abb). eg. ["ath","osa","dosa",.....]
    :return: own_maplist: The picked map num: eg. [map0001,map0002,.....,map1111]
    """
    if (len(org_list) == 1) and org_list[0] in ['Animals', 'Plants', 'Fungi']:
        tmp = get_html("https://rest.kegg.jp/list/organism")
        tmp = list(map(lambda x: x.split("\t"), tmp.strip().split("\n")))
        org_list = [line[1] for line in tmp if line[3].split(';')[1] == org_list[0]]        
    if not org_list == ["all"]:
        map_set = set()
        for org in org_list:
            tmp = get_html("https://rest.kegg.jp/list/pathway/"+org)
            tmp = list(map(lambda x: x.split("\t"), tmp.strip().split("\n")))
            tmp_df = pd.DataFrame(tmp, columns=["pathway", "description"]).dropna()
            map_set = map_set.union(set(tmp_df.pathway.apply(lambda x: re.search("[0-9]{5}", str(x)).group())))
        map_list = list(map(lambda x: "map"+x, list(map_set)))
        return map_list
    else:
        return None


def get_KEGGmap(own_maplist:list):
    """
    INPUT:own_maplist: The picked map num: eg. [map0001,map0002,.....,map1111]
    :return: KO2map_df:columns: [KO,pathway]
             KEGGmap_df: columns: [pathway, description]
    """
    # KEGG map description
    KEGGmap = get_html("https://rest.kegg.jp/list/pathway/")
    KEGGmap = list(map(lambda x: x.split("\t"), KEGGmap.split("\n")))
    KEGGmap_df = pd.DataFrame(KEGGmap, columns=["pathway", "description"]).dropna()
    KEGGmap_df["pathway"] = KEGGmap_df.pathway.apply(lambda x: x.strip("path:"))
    # KEGG KO2map
    KO2map = get_html("https://rest.kegg.jp/link/pathway/ko")
    KO2map = list(map(lambda x: x.split("\t"), KO2map.split("\n")))
    KO2map_df = pd.DataFrame(KO2map, columns=["ko", "pathway"]).dropna()
    KO2map_df["pathway"] = KO2map_df.pathway.apply(lambda x: x.strip("path:"))
    KO2map_df["KO"] = KO2map_df.ko.apply(lambda x: x.strip("ko:"))
    if own_maplist:
        KEGGmap_df = KEGGmap_df[KEGGmap_df["pathway"].isin(own_maplist)]
        KO2map_df = KO2map_df[KO2map_df["pathway"].isin(own_maplist)]
    return KO2map_df, KEGGmap_df


def parse_KO_f1(query_line):
    """
    :param query_line: a line in eggNOG annotation file which contains a query result
    :return:  a dict for parse_KO function. eg. {"gene":gene name,"KO": KO}
    """
    KO_list = [i for i in map(lambda x: x.lstrip("ko:"), query_line["KEGG_ko"].split(","))]
    return {"gene": query_line["Query"], "KO": KO_list}


def parse_KO(df4parse, org_list):
    """
    parse the KEGG part in eggNOG annotation file
    :param df4parse: the pd.Dataframe object directly comes from the eggNOG annotation file(skip the comment lines, of course)
    :return:the parsed KEGG annotation in pd.Dataframe.
    """
    KO2map, KEGGmap = get_KEGGmap(own_mapdb(org_list))
    gene2KO = df4parse[["Query", "KEGG_ko"]].dropna().apply(parse_KO_f1, axis=1, result_type="expand")
    gene2KO = pd.DataFrame({'gene': gene2KO.gene.repeat(gene2KO["KO"].str.len()), 'KO': np.concatenate(gene2KO["KO"].values)})
    gene2map = pd.merge(gene2KO, KO2map, on="KO", how="left")
    gene2map = pd.merge(gene2map, KEGGmap, on="pathway", how="left")
    return gene2map.dropna()


def parse_GO_f1(query_line):
    """
    :param query_line: a line in eggNOG annotation file which contains a query result
    :return:  a dict for parse_KO function. eg. {"gene":gene name,"GO": KO}
    """
    GO_list = [i for i in query_line["GOs"].split(",")]
    return {"gene": query_line["Query"], "GO": GO_list}


def parse_GO(df4parse, GO_path):
    """
    parse the GO part in eggNOG annotation file
    :param df4parse: the pd.Dataframe object directly comes from the eggNOG annotation file(skip the comment lines, of course)
    :return:the parsed GO annotation in pd.Dataframe.
    """
    gene2GO = df4parse[["Query", "GOs"]].dropna().apply(parse_GO_f1, axis=1, result_type="expand")
    gene2GO = pd.DataFrame(
        {'gene': gene2GO.gene.repeat(gene2GO["GO"].str.len()), 'GO': np.concatenate(gene2GO["GO"].values)})
    go_df = pd.read_table(GO_path)
    go_df.drop(columns=['Description'], inplace=True)
    gene2GO = pd.merge(gene2GO, go_df, on="GO", how="left")
    return gene2GO


def check_column(_file):
    with open(_file) as f_in:
        while True:
            line = f_in.readline()
            if not line.startswith("#"):
                items = line.split('\t')
                for _idx in range(len(items)):
                    if items[_idx].startswith("GO:"):
                        go_idx = _idx
                    elif items[_idx].startswith("ko:"):
                        kegg_idx = _idx
            if 'go_idx' in locals() and 'kegg_idx' in locals():
                return 0, go_idx, kegg_idx


def main(input_file:Path, out_dir: Path, go_file: Path, org_list:list):
    file4parse = pd.read_csv(input_file,
                             sep="\t",
                             comment="#",
                             usecols=check_column(input_file),
                             names=['Query', 'GOs', "KEGG_ko"])
    file4KO = parse_KO(file4parse, org_list)[["gene", "KO", "pathway", "description"]]
    file4GO = parse_GO(file4parse, go_file)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    file4KO.to_csv(os.path.join(out_dir, "KOannotation.tsv"), sep="\t", index=False)
    file4GO.to_csv(os.path.join(out_dir, "GOannotation.tsv"), sep="\t", index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This is the script for parsing eggNOGv2 annotation file to the file \
                                                  which ClusterProfile need")
    parser.add_argument('-i', '--input_file', required=True, default=Path,
                        help='<filepath>  The eggNOG annotation file')
    parser.add_argument('-g', '--go_file', required=True, default=Path,
                        help='<filepath>  The GO database file(come from parse_go_obofile.py)')
    parser.add_argument('-O', '--org_list', default="all", type=str,
                        help='<organism name>  The reference organisms in KEGG database(use KEGG abb. please) or one of Animals/Plants/Fungi. If not \
                        specified, it will use all pathway.')
    parser.add_argument('-o', '--output_directory', required=True, default=Path,
                        help='<dirpath>  the directory for results')
    args = parser.parse_args()
    main(args.input_file, args.output_directory, args.go_file, args.org_list.split(","))
