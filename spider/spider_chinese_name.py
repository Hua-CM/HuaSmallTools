# -*- coding: utf-8 -*-
# @Time : 2020/2/23 15:43
# @Author : Zhongyi Hua
# @FileName: spider_chinese_name.py
# @Usage: This script is used for spider Chinese name of plants from iplant.com
# @Note:
# @E-mail: njbxhzy@hotmail.com

import requests
import re
import argparse
from tqdm import tqdm
from pyquery import PyQuery as Pq
from time import sleep


def get_html(url):
    my_header = {
        'User-Agent': r'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:73.0) Gecko/20100101 Firefox/73.0'}
    return requests.get(url, headers=my_header).text


def get_chinese_name(latin_name, retry_num=0):
    retry_num += 1
    if retry_num > 3:
        return
    query_url = 'http://www.iplant.cn/info/' + latin_name
    try:
        a = Pq(get_html(query_url))
        try:
            chinese_name_inner = re.match('[\u4e00-\u9fa5]+', a('head>title').text()).group(0)
            if chinese_name_inner == "植物智":
                return
            return chinese_name_inner
        except AttributeError:
            synom_name = a('.infomore .spantxt').text()
            if not synom_name == '':
                return re.search('[\u4e00-\u9fa5]+', synom_name).group(0) + \
                       '\t' + \
                       re.search('[a-zA-Z\s]+', synom_name).group(0).strip()
            else:
                return
    except TimeoutError:
        sleep(1)
        get_chinese_name(latin_name, retry_num=retry_num)


def get_latin_name(chinese_name, retry_num=0):
    retry_num += 1
    if retry_num > 3:
        return
    query_url = 'http://www.iplant.cn/info/' + chinese_name
    try:
        a = Pq(get_html(query_url))
        latin_name = a('#sptitlel.infolatin').text()
        if latin_name == '':
            try:
                return re.search('[a-zA-Z\s]+', a('.infomore>a').text()).group(0).strip() + \
                       '\t' + \
                       re.search('[\u4e00-\u9fa5]+', a('.infomore>a').text()).group(0)
            except AttributeError:
                return
        return latin_name
    except TimeoutError:
        sleep(1)
        get_chinese_name(chinese_name, retry_num=retry_num)


def parse_args():
    parser = argparse.ArgumentParser(
        description='This is the script for spider Chinese name from iplants website.')
    parser.add_argument('-i', '--file_path', required=True,
                        help='<file_path>  A txt file contains one latin name per line')
    parser.add_argument('-m', '--model', default='latin', choices=['latin', 'chinese'],
                        help='<latin|chinese> "latin" for latin to chinese, '
                             '"chinese" for chinese to latin. DEFAULT:latin')
    parser.add_argument('-o', '--result_path', required=True,
                        help='<file_path>  The result file')
    args = parser.parse_args()
    return args


def main(args):
    func = get_chinese_name if args.model == 'latin' else get_latin_name
    result_list = []
    with open(args.file_path, 'r', encoding='utf-8') as query_file:
        query_list = query_file.read().split('\n')
    for query_name in tqdm(query_list):
        chinese_name = func(query_name)
        result_list.append(f'{query_name}\t{chinese_name}')
    with open(args.result_path, 'w', encoding='utf-8') as result_file:
        result_file.write('\n'.join(result_list))
        result_file.close()


if __name__ == '__main__':
    main(parse_args())
