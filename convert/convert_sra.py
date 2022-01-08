# -*- coding: utf-8 -*-
# @Time : 2020/9/15 10:45
# @Author : Zhongyi Hua
# @FileName: convert_sra.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import pandas as pd
import requests
from pyquery import PyQuery as PQ
from tqdm import tqdm
import re

headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:65.0) Gecko/20100101 Firefox/65.0"}


def get_html(url):
    resp = requests.get(url, headers=headers)
    return resp


def parse_ncbi(_input_file, _output_file):
    with open(_input_file, "r") as num_list:
        srx_par = re.compile(r"[DES]RX(\d{6,8})")
        srs_par = re.compile(r"[DES]RS(\d{6,8})")
        _result_list = []
        for line in tqdm(num_list):
            try:
                web_parser = PQ(get_html("https://www.ncbi.nlm.nih.gov/sra/?term=" + line.strip()).content, parser='html')
                _srp = web_parser('div #ResultView>div:eq(2)>span>div>a:eq(1)').text()
                _srx = re.search(srx_par, web_parser('p').text()).group()
                _srs = re.search(srs_par, web_parser('div #ResultView>div:eq(3)>span>div').text()).group()
                _org = web_parser('div #ResultView>div:eq(3)>div').text().replace('Organism: ', '')
                _srr, _spots, _bases, _size, *_ = web_parser('div #ResultView>table>tbody').text().split('\n')
                _result_dict = {'organism': _org, 'srp_acc': _srp, 'srx_acc': _srx, 'srs_acc': _srs, 'srr_acc': _srr,
                                'spots': _spots, 'bases(G)': _bases, 'size': _size}
                _result_list.append(_result_dict)
            except:
                continue
        pd.DataFrame(_result_list).to_csv(_output_file, sep='\t', index=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This is the script for get detailed information of SRR/SRX "
                                                 "accession,")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<filepath>  The file contain SRX/SRR num, One number per line please')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<filepath>  The output path')
    args = parser.parse_args()
    parse_ncbi(args.input_file, args.output_file)
