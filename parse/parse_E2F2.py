# -*- coding: utf-8 -*-
# @Time    : 2021/12/17 10:03
# @Author  : Zhongyi Hua
# @File    : parse_E2F2.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import pandas as pd
import argparse
from os.path import splitext


def parse_e2p2(e2p2, map1, map2):
    e2p2 = pd.read_table(e2p2, header=None, names=['geneid', 'EF'], comment="#").dropna()
    map1 = pd.read_table(map1, header=None, names=['EF', 'EC'])
    map2 = pd.read_table(map2)
    # split map1 to EF-Reacion ID and EF-EC
    map_reaction = map1[map1.EC.str.contains('RXN')].rename(columns={'EC': 'Reaction-id'})
    map_enzyme = map1[~map1.EC.str.contains('RXN')]
    map_enzyme.loc['EC'] = 'EC-' + map_enzyme['EC']
    # unstack e2p2
    e2p2 = pd.DataFrame(e2p2.EF.str.split('|').tolist(), index=e2p2.geneid).stack()
    e2p2 = e2p2.reset_index([0, 'geneid'])
    e2p2.columns = ['geneid', 'EF']
    # merge
    map_reaction = map_reaction.merge(map2, how='left').dropna()
    map_enzyme = map_enzyme.merge(map2, how='left').dropna()
    map1 = pd.concat([map_reaction, map_enzyme]).reset_index()
    map1.drop_duplicates(subset=['EF', 'Reaction-id', 'Pathway-id', 'Pathway-name', 'EC'], inplace=True)
    e2p2 = e2p2.merge(map1, how='inner')
    return e2p2


def parseArg():
    parser = argparse.ArgumentParser(description="Convert E2P2 result into table-format")
    parser.add_argument('result', help="E2P2 result")
    parser.add_argument('map1', help="E2P2 to Reaction ID/Enzyme")
    parser.add_argument('map2', help="Reaction ID/Enzyme to Pathway")
    args = parser.parse_args()
    return args


def main():
    args = parseArg()
    result = parse_e2p2(args.result, args.map1, args.map2)
    result.to_csv('_parsed'.join(splitext(args.result)), sep='\t', index=False)


if __name__ == '__main__':
    main()
