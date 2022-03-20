# -*- coding: utf-8 -*-
# @Time : 2021/12/8 11:07
# @Author : Zhongyi Hua
# @FileName: genewise2gff.py
# @Usage:
# @Note:
# @E-mail: njbxhzy@hotmail.com

import sys
import re

if len(sys.argv) < 2:
    sys.exit()

filepath = sys.argv[1]
fp = open(filepath, 'r')

choose = 0
for line in fp:
    if line.startswith("#"):
        continue
    elif re.search(r"\tmatch\t", line):
        choose = choose + 1
        score = line.strip().split('\t')[5]
    elif re.search(r"\tcds\t", line):
        CDS_List = line.strip().split('\t')
        CDS_List[2] = "nucleotide_to_protein_match"
        CDS_List[7] = '.'
        CDS_List[5] = score
        CDS_start = CDS_List[0].split(':')[1].split('-')[0]
        CDS_List[0] = CDS_List[0].split(':')[0]
        CDS_List[3] = str(int(CDS_start) + int(CDS_List[3]) - 1)
        CDS_List[4] = str(int(CDS_start) + int(CDS_List[4]) - 1)
        if CDS_List[6] == '-':
            CDS_List[3], CDS_List[4] = CDS_List[4], CDS_List[3]
        CDS_List[8] = "ID=match{:06d};Target=Dosnt_Matter{:06d}".format(choose, choose)
        print("\t".join(CDS_List))
    else:
        continue
fp.close()
