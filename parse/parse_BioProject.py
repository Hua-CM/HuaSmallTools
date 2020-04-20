# -*- coding: utf-8 -*-
# @Time : 2019/9/30 11:03
# @Author : Zhongyi Hua
# @FileName: parse_BioProject.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import re
import pandas as pd

biosample_pattern = {"biosample": "BioSample:(.*?)[;\n]",
                     "sra_id": "SRA:(.*?)[;\n]",
                     "host": "host=(.*?)\n",
                     "source": "source=(.*?)\n",
                     "sample_name": "Sample name:(.*?)[;\n]",
                     "organism": "Organism:(.*?)[;\n]"
                     }


def biosample(biosample_record):
    biosample_str = "".join(biosample_record)
    biosample_dict = {}
    for key, pattern in biosample_pattern.items():
        try:
            biosample_dict[key] = re.search(pattern, biosample_str).group(1).strip()
        except:
            biosample_dict[key] = "NA"
    return biosample_dict


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This is the script for extracting Biosample information from \
                                                  Bioproject txt. Since Bio could not parse its gbk result correctly")
    parser.add_argument('-i', '--input_txt', required=True,
                        help='<file_path> The Bioproject txt downloaded from NCBI')
    parser.add_argument('-o', '--output_tsv', required=True,
                        help='<file_path> The result table')
    args = parser.parse_args()
    parse_results = pd.DataFrame(columns=["biosample", "sra_id", "host", "source", "sample_name", "organism"])
    with open(args.input_txt, 'r', encoding='utf8') as f:
        cont = True
        li = []
        while cont:
            cont = f.readline()
            li.append(cont)
            if cont == '\n':
                if li == ["\n"]:
                    continue
                else:
                    parse_results = parse_results.append(biosample(li), ignore_index=True)
                    li = []
    parse_results.to_csv(args.output_tsv, sep="\t", index=None)
