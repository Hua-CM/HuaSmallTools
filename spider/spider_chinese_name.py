# -*- coding: utf-8 -*-
# @Time : 2020/2/23 15:43
# @Author : Zhongyi Hua
# @FileName: spider_chinese_name.py
# @Usage: This script is used for spider Chinese name of plants from iplant.com
# @Note:
# @E-mail: njbxhzy@hotmail.com

import requests
import re
from pyquery import PyQuery as pq


def get_html(url):
    my_header = {
        'User-Agent': r"Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:73.0) Gecko/20100101 Firefox/73.0"}
    return requests.get(url, headers=my_header).text


def get_chinese_name(latin_name):
    query_url = "http://www.iplant.cn/info/" + latin_name
    a = pq(get_html(query_url))
    try:
        chinese_name_inner = re.match("[\u4e00-\u9fa5]+", a('head>title').text()).group(0)
        if chinese_name_inner == "植物智":
            return
        return chinese_name_inner
    except AttributeError:
        synom_name = a('.infomore>a').text()
        if not synom_name == '':
            query_url = "http://www.iplant.cn/info/" + synom_name
            a = pq(get_html(query_url))
            try:
                chinese_name_inner = re.match("[\u4e00-\u9fa5]+", a('head>title').text()).group(0)
                return chinese_name_inner + "\t" + synom_name
            except AttributeError:
                return None
        else:
            return


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="This is the script for spider Chinese name from iplants website.")
    parser.add_argument('-i', '--file_path', required=True,
                        help='<file_path>  A txt file contains one latin name per line')
    parser.add_argument('-o', '--result_path', required=True,
                        help='<file_path>  The result file')
    args = parser.parse_args()
    result_list = []
    with open(args.file_path, "r") as query_file:
        while True:
            query_name = query_file.readline().strip()
            if query_name:
                chinese_name = get_chinese_name(query_name)
                result_list.append(f"{query_name}\t{chinese_name}")
            else:
                break
    with open(args.result_path, "w") as result_file:
        result_file.write("\n".join(result_list))
        result_file.close()
