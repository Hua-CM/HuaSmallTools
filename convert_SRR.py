# -*- coding: utf-8 -*-
# @Time : 2019/10/26 11:08
# @Author : Zhongyi Hua
# @FileName: convert_SRR.py
# @Usage:
# @Note:
# @E-mail: njbxhzy@hotmail.com
import requests
import argparse
import re
headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:65.0) Gecko/20100101 Firefox/65.0"}


def get_html(url):
    resp = requests.get(url, headers=headers)
    return resp.text


def convert_srr(input_file, output_file):
    with open(input_file, "r") as num_list:
        a = re.compile(r"[D,E,S]RR(\d{6,8})")
        for line in num_list:
            tmp_text = get_html("https://www.ncbi.nlm.nih.gov/sra/?term=" + line.strip())
            b = re.search(a, tmp_text).group()
            c = "/".join(["anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra", b[0:3], b[0:6], b, b+".sra"])
            with open(output_file, "a") as f:
                f.writelines(line.strip()+"\t"+b+"\t"+c+"\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This is the script for convert SRX num to SRR num, further more, it \
                                                  will generate the website for aspera software accordingly")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<filepath>  The file contain SRX num, One number per line please')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<filepath>  The output path')
    args = parser.parse_args()
    convert_srr(args.input_file, args.output_file)

