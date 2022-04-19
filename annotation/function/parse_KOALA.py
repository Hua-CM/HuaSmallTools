# -*- coding: utf-8 -*-
# @Time    : 2022/3/9 10:47
# @Author  : Zhongyi Hua
# @File    : parse_KOALA.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import pandas as pd

with open(r'E:\MHJ\02.annotation\KEGG\MHJ.KOALA', 'r') as f_in:
    lines = [_.split()[0:5] for _ in f_in.read().split('\n')]

ko_df = pd.DataFrame(lines, columns=['gene', 'KO', 'threshhold', 'score', 'evalue'])
ko_df.drop_duplicates(subset=['gene'], keep='first', inplace=True)

