# -*- coding: utf-8 -*-
# @Time : 2019/10/21 10:49
# @Author : Zhongyi Hua
# @FileName: parse_go_obofile.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com
import pandas as pd
import re

go_annotation = pd.DataFrame(columns=["GO", "Description", "Space"])
with open(r'F:\Tianma\transcriptome\go.obo', 'r', encoding='utf8') as f:
    cont = True
    cont_num = 0
    append_dict = {}
    while cont:
        cont = f.readline()
        cont_num += 1
        if cont == "\n":
            go_annotation = go_annotation.append(append_dict, ignore_index=True)
            append_dict = {}
            cont_num = 0
            continue
        else:
            try:
                if cont_num == 2:
                    append_dict.update({"GO": re.search("GO:[0-9]{7}", cont).group()})
                elif cont_num == 3:
                    append_dict.update({"Description": re.search("name: .*", cont).group()[6:]})
                elif cont_num == 4:
                    append_dict.update({"Space": re.search("namespace: .*", cont).group()[11:]})
            except AttributeError:
                if cont_num == 2:
                    append_dict.update({"GO": None})
                elif cont_num == 3:
                    append_dict.update({"Description": None})
                elif cont_num == 4:
                    append_dict.update({"Space": None})
    go_annotation.dropna().to_csv(r'F:\Tianma\transcriptome\go-basic.txt', sep="\t", index=False)

