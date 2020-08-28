# -*- coding: utf-8 -*-
# @Time : 2019/10/21 10:49
# @Author : Zhongyi Hua
# @FileName: parse_go_obofile.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com
import pandas as pd
import argparse
import re


def main(input_path, output_path):
    """
    :param input_path: The go OBO file path
    :param output_path: result file path.
    :return: None
    """
    go_annotation = pd.DataFrame(columns=["GO", "Description", "level"])
    with open(input_path, 'r') as f:
        cont = True
        cont_num = 0
        append_dict = {}
        dict_list = []
        while cont:
            cont = f.readline()
            cont_num += 1
            if cont == "\n":
                dict_list.append(append_dict)
                append_dict = {}
                cont_num = 0
                continue
            else:
                try:
                    if cont_num == 2:
                        a = re.search("GO:[0-9]{7}", cont).group()
                    elif cont_num == 3:
                        b = re.search("name: .*", cont).group()[6:]
                    elif cont_num == 4:
                        c = re.search("namespace: .*", cont).group()[11:]
                        append_dict = {"GO": a, "Description": b, "level": c}
                    elif re.match("alt_id", cont).group():
                        append_dict2 = {"GO": re.search("GO:[0-9]{7}", cont).group(), "Description": b, "level": c}
                        go_annotation = go_annotation.append(append_dict2, ignore_index=True)
                except AttributeError:
                    if cont_num == 2:
                        append_dict.update({"GO": None})
                    elif cont_num == 3:
                        append_dict.update({"Description": None})
                    elif cont_num == 4:
                        append_dict.update({"level": None})
        go_annotation = pd.DataFrame(dict_list)
        go_annotation = go_annotation[["GO", "Description", "level"]]
        go_annotation.dropna().to_csv(output_path, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This is the script for parsing GO database OBO file to a more \
                                                 user-friendly format")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<filepath>  The GO database file in OBO format')
    parser.add_argument('-o', '--output_file', required=True,
                        help='<filepath>  The result path')
    args = parser.parse_args()
    main(args.input_file, args.output_file)
