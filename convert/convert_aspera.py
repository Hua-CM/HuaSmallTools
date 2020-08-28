# -*- coding: utf-8 -*-
# @Time : 2019/10/26 11:08
# @Author : Zhongyi Hua
# @FileName: convert_aspera.py
# @Usage:
# @Note:
# @E-mail: njbxhzy@hotmail.com
import requests
import argparse
import re
from itertools import chain
headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:65.0) Gecko/20100101 Firefox/65.0"}


def get_html(url):
    resp = requests.get(url, headers=headers)
    return resp.text


def convert_ncbi(input_file, output_file):
    with open(input_file, "r") as num_list:
        a = re.compile(r"[DES]RR(\d{6,8})")
        for line in num_list:
            tmp_text = get_html("https://www.ncbi.nlm.nih.gov/sra/?term=" + line.strip())
            b = re.search(a, tmp_text).group()
            c = "/".join(["anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra", b[0:3], b[0:6], b, b+".sra"])
            with open(output_file, "a") as f:
                f.writelines(line.strip()+"\t"+b+"\t"+c+"\n")


def convert_ena(input_file, output_file):
    with open(input_file, "r") as num_list:
        for sra_number in num_list:
            srr_list = get_html(f"https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={sra_number.strip()}&result=read_run&fields=fastq_ftp")
            srr_list = list(chain.from_iterable([x.split('\t')[1].split(";") for x in srr_list.split("\n")[1:-1]]))
            srr_list = [re.sub("ftp\\.sra\\.ebi\\.ac\\.uk", "era-fasp@fasp.sra.ebi.ac.uk:", x) for x in srr_list]
            with open(output_file, "a") as f:
                f.writelines("\n".join([sra_number.strip() + "\t" + x for x in srr_list[:-1]]) + "\n")


def convert_ddbj(input_file, output_file):
    with open(input_file, "r") as num_list:
        a = re.compile(r"[DES]RR(\d{6,8})")
        for line in num_list:
            tmp_text = requests.post('https://ddbj.nig.ac.jp/DRASearch/', {'acc': line.strip()}, headers=headers).text
            b = re.search(a, tmp_text).group()
            c = "/".join(['anonftp@ascp.ddbj.nig.ac.jp:/ddbj_database/dra/fastq', b[0:3], b[0:6], b, b+'.sra'])
            with open(output_file, 'a') as f:
                f.writelines(line.strip()+'\t'+b+'\t'+c+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This is the script for convert SRX num to SRR num, further more, it \
                                                  will generate the URL for aspera software accordingly (SRR accession \
                                                  is also acceptable for generating aspera links")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<filepath>  The file contain SRX/SRR num, One number per line please')
    parser.add_argument('-t', '--type', default="NCBI", choices=['NCBI', 'ENA', 'DDBJ'],
                        help='<NCBI|ENA|DDBJ> Which source you want to download data? Default is NCBI')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<filepath>  The output path')
    args = parser.parse_args()
    if args.type == 'NCBI':
        convert_ncbi(args.input_file, args.output_file)
    elif args.type == 'ENA':
        convert_ena(args.input_file, args.output_file)
    elif args.type == 'DDBJ':
        convert_ddbj(args.input_file, args.output_file)
