# -*- coding: utf-8 -*-
# @Time : 2020/9/15 10:45
# @Author : Zhongyi Hua
# @FileName: convert_sra.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import pandas as pd
import requests
from pyquery import PyQuery as PQ
import xml.etree.ElementTree as ET
from tqdm import tqdm
import re
import argparse
from time import sleep


headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:65.0) Gecko/20100101 Firefox/65.0"}
URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'



def get_html(url):
    resp = requests.get(url, headers=headers)
    return resp

# This is an outdate function
def parse_ncbi(_input_file, _output_file):
    """_summary_

    Args:
        _input_file (_type_): _description_
        _output_file (_type_): _description_
    """
    SRX_PAR = re.compile(r"[DES]RX(\d{6,8})")
    SRS_PAR = re.compile(r"[DES]RS(\d{6,8})")
    with open(_input_file, "r") as num_list:
        _result_list = []
        for line in tqdm(num_list):
            try:
                web_parser = PQ(get_html("https://www.ncbi.nlm.nih.gov/sra/?term=" + line.strip()).content, parser='html')
                _srp = web_parser('div #ResultView>div:eq(2)>span>div>a:eq(1)').text()
                _srx = re.search(SRX_PAR, web_parser('p').text()).group()
                _srs = re.search(SRS_PAR, web_parser('div #ResultView>div:eq(3)>span>div').text()).group()
                _org = web_parser('div #ResultView>div:eq(3)>div').text().replace('Organism: ', '')
                _srr, _spots, _bases, _size, *_ = web_parser('div #ResultView>table>tbody').text().split('\n')
                _result_dict = {'organism': _org, 'srp_acc': _srp, 'srx_acc': _srx, 'srs_acc': _srs, 'srr_acc': _srr,
                                'spots': _spots, 'bases(G)': _bases, 'size': _size}
                _result_list.append(_result_dict)
            except:
                continue
        pd.DataFrame(_result_list).to_csv(_output_file, sep='\t', index=False)


class ExpXMLElement:
    def __init__(self, tag, attr=None, value=None) -> None:
        """Create a class for EXPERIMENT_PACKAGE_SET XML information storage

        Args:
            key (str): The element key
            attr (dict, optional): The informaiton in '<>'. Defaults to None.
            value (str, List[ExpXMLElement], optional): The informaiton between '<></>'. Defaults to None.
        """
        self.tag = tag
        self.attr = attr
        self.value = value


def parse_XML(ExpXML :str):
    """
    Parse the EXPERIMENT_PACKAGE_SET XML.
    Args:
        ExpXML (str): The EXPERIMENT_PACKAGE_SET XML
    Return:
        _res_lst (lst): A list of ExpXMLElement}

    """
    ExpXML = ExpXML.strip()
    _res_lst = []
    _idx = 0
    while _idx < len(ExpXML) - 1:
        _key, _value = ('', '')
        while ExpXML[_idx] != '>':
            _key += ExpXML[_idx]
            _idx += 1
        # for the ">"
        _key += ExpXML[_idx]
        _idx += 1
        _key = _key.strip(' ') # The key with <>
        _real_key = _key.split(' ')[0].strip('<>')
        tmp_ele =  ExpXMLElement(_real_key)
        # Parse attributes first
        if ' ' in _key:
            _attrs_lst = _key.split('" ')
            _attrs_lst[0] = _attrs_lst[0].replace(_real_key+ ' ','')
            _attrs_dict = {}
            for _attr in _attrs_lst:
                __ = _attr.split('=')
                _attrs_dict[__[0].strip('<>"')]  = __[1].strip('"/')
            tmp_ele.attr = _attrs_dict
        if '/>' not in _key:
        # This is a item with text, find another <>
            end_key = '</' + _real_key + '>'
            while ExpXML[_idx: _idx+len(end_key)] != end_key:
                _value += ExpXML[_idx]
                _idx += 1
            if '<' in _value:
                # This is a nested item
                tmp_ele.value = parse_XML(_value.strip())
            else:
                # This is a single item
                tmp_ele.value = _value
            _idx += len(end_key)
        _res_lst.append(tmp_ele)
    return _res_lst


def convert_base(base: str):
    """Convert raw nucleotide base to readable format

    Args:
        base (str): a numeric string

    Returns:
        str : a string with Gb, Mb or kb
    """
    base_dict = {
     0: 'b',
     1: 'kb',
     2: 'Mb',
     3: 'Gb'
    }
    base = int(base)
    div_time = 0
    while (base > 1000) and (div_time < 3):
        div_time += 1
        base /= 1000
    out_base = str(round(base, 3)) + base_dict.get(div_time)
    return out_base


def convert_size(size: str):
    """Convert raw computer storage size to readable format

    Args:
        size (str): a numeric string

    Returns:
        str : a string with GB, MB or KB
    """
    size_dict = {
     0: 'B',
     1: 'KB',
     2: 'MB',
     3: 'GB'
    }
    size = int(size) / 8
    div_time = 0
    while (size > 1000) and (div_time < 3):
        div_time += 1
        size /= 1000
    out_size = str(round(size,3)) + size_dict.get(div_time)
    return out_size


def parse_sra(_query_list):
    """Parse a list of SRA accession

    Args:
        _query_list (List[str]): A list of SRR/SRX accession.

    Returns:
        pd.DataFrame: The result DataFrame with the following information
        'ScientificName', 'taxid', 'SRP', 'SRX', 'SRS', 'SRR', 'Bioproject', 'Biosample', 'Bases', 'Size'
    """
    url = URL + 'efetch.fcgi?db=sra&id=' + ','.join(_query_list) + '&rettype=docsum&version=2.0'
    _efetch_result = ET.fromstring(get_html(url).text)
    _result_lst = []
    for _sra_item in _efetch_result:
        try:
            sra_res_dict = {}
            tmp_res1 = None
            tmp_res2 = None
            for _ele in _sra_item:
                if _ele.attrib.get('Name') == "ExpXml":
                    tmp_res1 = parse_XML(_ele.text)
                if _ele.attrib.get('Name') == "Runs":
                    tmp_res2 = parse_XML(_ele.text)[0]
            if not tmp_res1:
                continue
            for _expele in tmp_res1:
                if _expele.tag == 'Summary':
                    for _sub_ele in _expele.value:
                        if _sub_ele.tag == 'Platform':
                            sra_res_dict['Platform'] =  _sub_ele.attr.get('instrument_model')
                        elif _sub_ele.tag == 'Statistics':
                            sra_res_dict['Bases'] =  convert_base(_sub_ele.attr.get('total_bases'))
                            sra_res_dict['Size']  =  convert_size(_sub_ele.attr.get('total_size'))
                elif _expele.tag == 'Organism':
                    sra_res_dict['taxid'] = _expele.attr.get('taxid')
                    sra_res_dict['ScientificName'] = _expele.attr.get('ScientificName')
                elif _expele.tag == 'Bioproject':
                    sra_res_dict['Bioproject'] = _expele.value
                elif _expele.tag == 'Biosample':
                    sra_res_dict['Biosample'] = _expele.value
                elif _expele.tag == 'Sample':
                    sra_res_dict['SRS'] = _expele.attr.get('acc')
                elif _expele.tag == 'Experiment':
                    sra_res_dict['SRX'] = _expele.attr.get('acc')
                elif _expele.tag == 'Study':
                    sra_res_dict['SRP'] = _expele.attr.get('acc')
            sra_res_dict['SRR'] = tmp_res2.attr.get('acc')
            _result_lst.append(sra_res_dict)
        except:
            print(_sra_item.text)
    res_df = pd.DataFrame(_result_lst)
    res_df = res_df[['ScientificName', 'taxid', 'SRP', 'SRX', 'SRS', 'SRR', 'Bioproject', 'Biosample', 'Bases', 'Size']]
    return res_df

def parse_args():
    """_summary_
    Parse command line arguments
    Returns:
        args : the arguments
    """
    parser = argparse.ArgumentParser(description="This is the script for get detailed information of SRR/SRX "
                                                 "accession,")
    parser.add_argument('-i', '--input', required=True,
                        help='<filepath>  The file contain SRX/SRR num, One number per line please')
    parser.add_argument('-o', '--output', required=True,
                        help='<filepath>  The output path')
    args = parser.parse_args()
    return args


def main(args):
    df_lst = []
    with open(args.input, "r") as f_in:
        acc_lst = f_in.read().strip().splitlines()
    batch_list = [acc_lst[i:i + 100] for i in range(0, len(acc_lst), 100)]
    for query_lst in batch_list:
        df_lst.append(parse_sra(query_lst))
        sleep(1)
    res_df = pd.concat(df_lst)
    res_df.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main(parse_args())

    
