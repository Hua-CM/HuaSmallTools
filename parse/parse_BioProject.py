# -*- coding: utf-8 -*-
# @Time : 2019/9/30 11:03
# @Author : Zhongyi Hua
# @FileName: parse_BioProject.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import re
import pandas as pd


class Biosample:
    def __init__(self, biosample_str):
        try:
            self.biosample_id = re.search(r"SAMN[0-9]{8}", str(biosample_str)).group()
        except AttributeError:
            self.biosample_id = "NA"
        try:
            self.sra_id = re.search(r"SRS[0-9]{7}", str(biosample_str)).group()
        except AttributeError:
            self.sra_id = "NA"
        try:
            self.host = re.search("host=(.*?)\\\\n", str(biosample_str)).group()[6:-3]
        except AttributeError:
            self.host = "NA"
        try:
            self.source = re.search("source=(.*?)\\\\n", str(biosample_str)).group()[8:-3]
        except AttributeError:
            self.source = "NA"
        try:
            self.sample_name = re.search("Sample name:(.*?);", str(biosample_str)).group()[13:-1]
        except AttributeError:
            self.sample_name = "NA"


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This is the script for extracting Biosample information from \
                                                  Bioproject txt. Since Bio could not parse its gbk result correctly")
    parser.add_argument('-i', '--input_txt', required=True,
                        help='<file_path> The Bioproject txt downloaded from NCBI')
    parser.add_argument('-o', '--output_tsv', required=True,
                        help='<file_path> The result table')
    args = parser.parse_args()
    parse_results = pd.DataFrame(columns=["biosample_id", "sra_id", "host", "source", "sample_name"])
    with open(args.input_txt_file, 'r', encoding='utf8') as f:
        cont = True
        li = []
        while cont:
            cont = f.readline()
            li.append(cont)
            if cont == '\n':
                print(li)
                if li == ["\n"]:
                    continue
                else:
                    parse_results = parse_results.append(vars(Biosample(li)), ignore_index=True)
                    li = []
    parse_results.to_csv(args.output_tsv, sep="\t", index=None)
