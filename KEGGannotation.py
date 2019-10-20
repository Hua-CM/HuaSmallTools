# -*- coding: utf-8 -*-
# @Time : 2019/10/20 16:56
# @Author : Zhongyi Hua
# @FileName: KEGGannotation.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com
import requests
import pandas as pd
import re
headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:65.0) Gecko/20100101 Firefox/65.0"}


def get_html(url):
    resp = requests.get(url, headers=headers)
    return resp.text


def own_mapdb(org_list):
    map_set = set()
    for org in org_list:
        tmp = get_html("http://rest.kegg.jp/list/pathway/"+org)
        tmp = map(lambda x: str.split(x, "\t"), tmp.split("\n"))
        tmp_df = pd.DataFrame(tmp, columns=["pathway", "description"]).dropna()
        map_set = map_set.union(set(tmp_df.pathway.apply(lambda x: re.search("\d{5}", str(x)).group())))
    map_list = list(map(lambda x: "map"+x, list(map_set)))
    return map_list


def own_map2decription(own_maplist):
    tmp = get_html("http://rest.kegg.jp/list/pathway/")
    tmp = map(lambda x: str.split(x, "\t"), tmp.split("\n"))
    tmp_df = pd.DataFrame(tmp, columns=["pathway", "description"]).dropna()
    tmp_df["pathway"] = tmp_df.pathway.apply(lambda x: str.strip(x, "path:"))
    tmp_df = tmp_df[tmp_df["pathway"].isin(own_maplist)]
    tmp_df.to_csv(r"E:\Tianma_data\genome\KEGGpathwayInfo", sep="\t", index=False)
    return tmp_df


def ko2map(own_maplist):
    tmp = get_html("http://rest.kegg.jp/link/pathway/ko")
    tmp = map(lambda x: str.split(x, "\t"), tmp.split("\n"))
    tmp_df = pd.DataFrame(tmp, columns=["ko", "pathway"]).dropna()
    tmp_df["pathway"] = tmp_df.pathway.apply(lambda x: str.strip(x, "path:"))
    tmp_df["ko"] = tmp_df.ko.apply(lambda x: str.strip(x, "ko:"))
    tmp_df = tmp_df[tmp_df["pathway"].isin(own_maplist)]
    tmp_df.to_csv(r"E:\Tianma_data\genome\KO2pathway.map", sep="\t", index=False)
    return tmp_df


def gene2map(ko2map_df):
    gene2ko = pd.read_csv(r"E:\Tianma_data\genome\gastrodia_ko.txt", sep="\t", names=["gene", "ko"]).dropna()
    gene2ko["ko"] = gene2ko.ko.apply(lambda x: str.strip(x, "ko:"))
    gene2map_df = pd.merge(gene2ko, ko2map_df, on="ko", how="left").dropna()
    gene2map_df = gene2map_df[["pathway", "gene"]]
    return gene2map_df
