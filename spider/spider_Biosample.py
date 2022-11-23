# -*- coding: utf-8 -*-
# @Time : 2019/9/30 11:03
# @Author : Zhongyi Hua
# @FileName: parse_BioProject.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import pandas as pd
import requests
import xml.etree.ElementTree as ET
import argparse
from time import sleep

from spider_sra import parse_XML

headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:65.0) Gecko/20100101 Firefox/65.0"}
URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def get_html(url):
    resp = requests.get(url, headers=headers)
    return resp


def parse_biosample(query_lst):
    """Parse a list of Biosample accession (e.g. SAMD00009395)

    Args:
        _query_list (List[str]): A list of BioSample accession.

    Returns:
        pd.DataFrame: The result DataFrame with the following information
        'Biosample', 'taxid', 'taxonomy', 'longitude', 'latitude', 'description'
    """
    url = URL + 'efetch.fcgi?db=BioSample&id=' + ','.join(query_lst) + '&rettype=docsum&version=2.0'
    _efetch_result = ET.fromstring(get_html(url).text)
    _result_lst = []
    for _bsample_item in _efetch_result[0]:
        try:
            if _bsample_item.tag == 'DbBuild':
                continue
            bsample_res_dict = {}
            for _ele in _bsample_item:
                if _ele.tag=='Accession':
                    bsample_res_dict['Biosample'] = _ele.text
                elif _ele.tag=='Title':
                    bsample_res_dict['Description'] = _ele.text
                elif _ele.tag=='Taxonomy':
                    bsample_res_dict['Taxid'] = _ele.text
                elif _ele.tag=='Organism':
                    bsample_res_dict['Taxonomy'] = _ele.text
                elif _ele.tag == "SampleData":
                    tmp_res1 = parse_XML(_ele.text)[0].value
            # Not all samples has latitude and longitude
            latitude, longitude = '', ''
            for _ele in tmp_res1:
                if _ele.tag == 'Attributes':
                    for _subele in _ele.value:
                        if _subele.attr.get('attribute_name') == "lat_lon":
                            _tmp_str = _subele.value.split()
                            latitude, longitude = ' '.join(_tmp_str[0:2]), ' '.join(_tmp_str[2:])
                        bsample_res_dict['Longitude'] = longitude
                        bsample_res_dict['Latitude'] = latitude
            _result_lst.append(bsample_res_dict)
        except:
            print(_bsample_item.text)
    res_df = pd.DataFrame(_result_lst)
    return res_df


def convert_EMBL(accession_lst):
    """_summary_
    Since there were both SAMD00064000 and SAMN00064000 in Biosample database.
    The Biosample from EMBL and DDBJ databases must use uid (not accession). This function convert
    accession ids to uids. 

    Args:
        accession_lst (_type_): A list of accession ids.
    Return:
        list: A list of uids
    """
    url = URL + 'esearch.fcgi?db=Biosample&term=SAMD00064000|SAMD00070131' + '|'.join(accession_lst)
    _esearch_result = ET.fromstring(get_html(url).text)
    _result_lst = []
    for _item in _esearch_result:
        if _item.tag == 'IdList':
            for _id in _item:
                _result_lst.append(_id.text)
    return _result_lst


def parse_args():
    """Parse arguments
    """
    parser = argparse.ArgumentParser(description="This is the script for get detailed information of BioSample")
    parser.add_argument('-i', '--input', required=True,
                        help='<file_path> The file contain Biosample num, One number per line please')
    parser.add_argument('-o', '--output', required=True,
                        help='<file_path> The result table')
    args = parser.parse_args()
    return args

def main(args):
    with open(args.input, "r") as f_in:
        tmp_lst = f_in.read().strip().splitlines()
    # convert to uid first
    convert_lst = []
    acc_lst = []
    for acc in tmp_lst:
        if acc.startswith('SAMD') or acc.startswith('SAMEA'):
            convert_lst.append(acc)
        else:
            acc_lst.append(acc)
    batch_list = [convert_lst[i:i + 20] for i in range(0, len(convert_lst), 20)]
    converted_lst = []
    for query_lst in batch_list:
        converted_lst += convert_EMBL(query_lst)
        sleep(0.1)
    acc_lst += converted_lst

    batch_lst = [acc_lst[i:i + 20] for i in range(0, len(acc_lst), 20)]
    df_lst = []
    for query_lst in batch_lst:
        df_lst.append(parse_biosample(query_lst))
        sleep(0.1)
    res_df = pd.concat(df_lst)
    res_df = res_df[['Biosample', 'Taxid',	'Taxonomy',	'Longitude', 'Latitude', 'Description']]
    res_df.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main(parse_args())
