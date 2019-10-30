# -*- coding: utf-8 -*-
# @Time : 2019/10/20 16:56
# @Author : Zhongyi Hua
# @FileName: parse_KEGGkolist.py
# @Usage: This is the script for parsing KOALA web server results with pathway \
#         information. The result file could be used in clusterProfiler
# @Note:
# @E-mail: njbxhzy@hotmail.com
import requests
import pandas as pd
import re
import argparse
headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:65.0) Gecko/20100101 Firefox/65.0"}


def get_html(url):
    resp = requests.get(url, headers=headers)
    return resp.text


def own_map2decription(org_list):
    """
    Use REST API to generate map info in orgnism list
    :param org_list: organism name list (KEGG abb). eg. ["ath","osa","dosa",.....]
    :return: The pathways and their info in organism list. columns["pathway" , "description"]
    """
    tmp = get_html("http://rest.kegg.jp/list/pathway/")
    tmp = map(lambda x: str.split(x, "\t"), tmp.split("\n"))
    map2info_df = pd.DataFrame(tmp, columns=["pathway", "description"]).dropna()
    map2info_df["pathway"] = map2info_df.pathway.apply(lambda x: str.strip(x, "path:"))
    if org_list:
        map_set = set()
        for org in org_list:
            tmp = get_html("http://rest.kegg.jp/list/pathway/"+org)
            tmp = map(lambda x: str.split(x, "\t"), tmp.split("\n"))
            tmp_df = pd.DataFrame(tmp, columns=["pathway", "description"]).dropna()
            map_set = map_set.union(set(tmp_df.pathway.apply(lambda x: re.search("[0-9]{5}", str(x)).group())))
        map_list = list(map(lambda x: "map"+x, list(map_set)))
        map2info_df = map2info_df[map2info_df["pathway"].isin(map_list)]
    return map2info_df


def ko2map():
    """
    Use REST API to generate ko-map id mapping file.
    :return:The ko-map id mapping file. columns["pathway" , "ko"]
    """
    tmp = get_html("http://rest.kegg.jp/link/pathway/ko")
    tmp = map(lambda x: str.split(x, "\t"), tmp.split("\n"))
    ko2map_df = pd.DataFrame(tmp, columns=["ko", "pathway"]).dropna()
    ko2map_df["pathway"] = ko2map_df.pathway.apply(lambda x: str.strip(x, "path:"))
    ko2map_df["ko"] = ko2map_df.ko.apply(lambda x: str.strip(x, "ko:"))
    return ko2map_df


def gene2map(input_file, output_file, ko2map_df, map2info_df):
    """
    main function
    :param input_file: input file path
    :param output_file: output file path
    :param ko2map_df: ko-map id mapping
    :param map2info_df: map-description (mapinfo)
    :return: The file clusterProfiler needed. columns["gene","pathway" , "description"]
    """
    gene2ko = pd.read_csv(input_file, sep="\t", names=["gene", "ko"]).dropna()
    gene2ko["ko"] = gene2ko.ko.apply(lambda x: str.strip(x, "ko:"))
    gene2map_df = pd.merge(gene2ko, ko2map_df, on="ko", how="left").dropna()
    gene2map_df = pd.merge(gene2map_df, map2info_df, on="pathway", how="right").dropna()
    gene2map_df = gene2map_df[["gene", "pathway", "description"]]
    gene2map_df.to_csv(output_file, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This is the script for parsing KOALA web server results with pathway \
                                                  information. The result file could be used in clusterProfiler ")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<filepath> The BlastKOALA web server result file')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<output path> The output file')
    parser.add_argument('-O', '--org_list', nargs="?", const=1, default="all",
                        help='<organism name>  The reference organisms in KEGG database(use KEGG abb. please).If not \
                        specified, it will use all pathway.')
    args = parser.parse_args()
    mapinfo = own_map2decription(args.org_list.split(","))
    gene2map(args.input_file, args.output_file, ko2map(), mapinfo)
