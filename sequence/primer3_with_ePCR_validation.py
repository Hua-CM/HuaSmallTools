# -*- coding: utf-8 -*-
# @Time : 2022/1/7 0:48
# @Author : Zhongyi Hua
# @FileName: primer3_with_ePCR_validation.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com
import pandas as pd
from Bio import SeqIO
from tempfile import mktemp
import subprocess as sp
import argparse
import os
from collections import deque


def del_dir(_dir):
    for r, d, f in os.walk(_dir):
        for files in f:
            os.remove(os.path.join(r, files))
        os.removedirs(r)


def get_seq(_seq, _start, _end, _product):
    """
    :param _seq: Bio.SeqRecord.SeqRecord
    :param _start: int
    :param _end: int
    :param _product: the max product size,int
    :return: seq
    """
    trim_start = _start - _product - 1 if _start - _product - 1 > 0 else 0
    trim_end = _end + _product if _end + _product < _seq.__len__() else _seq.__len__()
    template_seq = _seq[trim_start: trim_end]
    template_seq.id = _seq.id + '_{}_{}'.format(_start, _end)
    return template_seq


class Primer3:
    def __init__(self, _infodf, _tmp, _primer_bin):
        """
        :param _infodf: pd.Dataframe contain AT LEAST five columns: 'seqid', 'start', 'end', 'minsize', 'maxsize'
        :param _tmp:
        :param _primer_bin:
        """
        self.infodf = _infodf
        self.bin = _primer_bin
        self.p3in = mktemp(dir=_tmp)
        self.p3out = mktemp(dir=_tmp)

    @staticmethod
    def readfile(_file):
        with open(_file) as f_in:
            return f_in.read().splitlines()

    def prepare_p3in(self, _fasta, _add):
        """
        _fasta: a dict from SeqIO.to_dict
        """
        _write_list = []
        for _item in self.infodf.itertuples():
            _seq = get_seq(_fasta[_item.seqid], _item.start, _item.end, _item.maxsize)
            _sequence_template = 'SEQUENCE_TEMPLATE=' + str(_seq.seq)
            _write_list.append('SEQUENCE_ID=' + _seq.id)
            _write_list.append(_sequence_template)
            _write_list.append('SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,{},{},{}'.
                               format(_item.maxsize, _item.maxsize+_item.end-_item.start, _item.maxsize)
                               )
            _write_list.append('PRIMER_PRODUCT_SIZE_RANGE={}-{}'.format(_item.minsize, _item.maxsize))
            if _add:
                _write_list += self.readfile(_add)
            _write_list.append('=')
        with open(self.p3in, 'a') as f2:
            f2.write('=\n')
            f2.write('\n'.join(_write_list))
            f2.write("\n")

    def run_primer3(self, _config):
        # check whether product size set in config file
        pass
        if _config:
            sp.run([self.bin, '--p3_settings_file=%s' % _config, '--output=%s' % self.p3out, self.p3in])
        else:
            sp.run([self.bin,  '--output=%s' % self.p3out, self.p3in])

    def p3out2sts(self):
        _file = deque(self.readfile(self.p3out))
        _module = []
        _record_list = []
        for _line in _file:
            _module.append(_line)
            if _line == '=':
                _model_dict = dict(_.split('=') for _ in _module)
                if 'SEQUENCE_ID' in _model_dict:
                    if not _model_dict['PRIMER_PAIR_NUM_RETURNED'] == '0':
                        for _num in range(int(_model_dict['PRIMER_PAIR_NUM_RETURNED'])):
                            _record = {'ID': _model_dict['SEQUENCE_ID'],
                                       'GenomePosition': _model_dict['SEQUENCE_ID'].split('_')[1].strip(':.'),
                                       'ProductSize2': _model_dict['PRIMER_PAIR_{}_PRODUCT_SIZE'.format(_num)],  # avoid a merge function
                                       'Forward': _model_dict['PRIMER_LEFT_{}_SEQUENCE'.format(_num)],
                                       'Forward_position': _model_dict['PRIMER_LEFT_{}'.format(_num)].split(',')[0],
                                       'Forward_length': _model_dict['PRIMER_LEFT_{}'.format(_num)].split(',')[1],
                                       'Forward_TM': _model_dict['PRIMER_LEFT_{}_TM'.format(_num)],
                                       'Reverse': _model_dict['PRIMER_RIGHT_{}_SEQUENCE'.format(_num)],
                                       'Reverse_position': _model_dict['PRIMER_RIGHT_{}'.format(_num)].split(',')[0],
                                       'Reverse_length': _model_dict['PRIMER_RIGHT_{}'.format(_num)].split(',')[1],
                                       'Reverse_TM': _model_dict['PRIMER_RIGHT_{}_TM'.format(_num)],
                                       }
                            _record_list.append(_record)
            else:
                continue
        tb_p3out = pd.DataFrame(_record_list)
        _size_dict = {_['ID']: _['maxsize'] for _ in self.infodf.to_dict(orient='records')}
        _start_dict = {_['ID']: _['start'] for _ in self.infodf.to_dict(orient='records')}
        # fix primer position
        tb_p3out['Forward_position2'] = tb_p3out.apply(
            lambda x: int(x['Forward_position']) + _start_dict.get(x['ID']) - _size_dict.get(x['ID']), axis=1)
        tb_p3out['Reverse_position2'] = tb_p3out.apply(
            lambda x: int(x['Reverse_position']) + _start_dict.get(x['ID']) - _size_dict.get(x['ID']), axis=1)
        # generate sts
        product_dict = {_['ID']: _['ProductSize'] for _ in self.infodf.to_dict(orient='records')}
        seqid = ''
        primer_count = 0
        sts_list = []
        for _idx, _item in tb_p3out.iterrows():
            if seqid == _item['ID']:
                primer_count += 1
            else:
                seqid = _item['ID']
                primer_count = 0
            sts_list.append({'primerid': '{}_primer{}'.format(seqid, primer_count),
                             'Forward': _item['Forward'],
                             'Reverse': _item['Reverse'],
                             'size': product_dict.get(seqid)})
        tb_sts = pd.DataFrame(sts_list)
        return tb_p3out, tb_sts


class rePCR:
    def __init__(self, _sts, _fahash,  _tmp, _repcr_bin):
        """
        :param _sts: IMPORTANCE The *.sts file for re-PCR must contain the product size!!!
        """
        self.sts = _sts
        self.fahash = _fahash
        self.bin = _repcr_bin
        self.out = mktemp(dir=_tmp)

    def validate(self):
        sp.run('{} -S {} -g 1 -n 1 -o {} {}'.format(self.bin, self.fahash, self.out, self.sts), shell=True)
        # parse results
        tb_repcr_result = pd.read_table(self.out,
                                        comment='#',
                                        names=['primerid', 'seqid', 'strand', 'start', 'end', 'mism', 'gaps', 'len'])
        # Scheme I
        """
        repcr_dict_list = tb_repcr_result.groupby('primerid').size().reset_index(name='count').to_dict('records')
        seqid_dict = defaultdict()
        primer_list = []
        for repcr_dict in repcr_dict_list:
            seqid_dict.setdefault(re.sub('_primer[0-9]+', '', repcr_dict.get('primerid')), []).append(repcr_dict)
        for _seqid, _repcr_list in seqid_dict.items():
            min_primer = min(_repcr_list, key=lambda d: d['count'])
            primer_list.append({'ID': _seqid, 'primerid': min_primer.get('primerid'), 'count': min_primer.get('count')})
        return pd.DataFrame(primer_list)
        """
        # Scheme II
        return tb_repcr_result.groupby('primerid').size().reset_index(name='count')


def main(args):
    FASTA = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))
    INFO = pd.read_table(args.info,
                         header=None,
                         names=['seqid', 'start', 'end', 'minsize', 'maxsize'],
                         dtype={'seqid': 'str', 'start': 'int64', 'end': 'int64', 'minsize': 'int64', 'maxsize': 'int64'})
    INFO['ID'] = INFO['seqid'] + '_' + INFO['start'].map(str) + '_' + INFO['end'].map(str)
    INFO['ProductSize'] = INFO['minsize'].map(str) + '-' + INFO['maxsize'].map(str)
    # Set temporary directory
    if (not args.tmp == '') and (not os.path.exists(args.tmp)):
        os.mkdir(args.tmp)
    elif args.tmp == '':
        args.tmp = mktemp()
    primer3 = Primer3(INFO, args.tmp, args.primer3)
    primer3.prepare_p3in(FASTA, args.add)
    primer3.run_primer3(args.config)
    tb_p3out, tb_sts = primer3.p3out2sts()
    sts_file_path = mktemp(dir=args.tmp)
    tb_sts.to_csv(sts_file_path, sep='\t', index=False, header=False)
    repcr = rePCR(sts_file_path, args.fahash, args.tmp, args.epcr)
    tb_validate = repcr.validate()
    tb_final = tb_validate.merge(tb_sts, how='left').merge(tb_p3out, how='left').merge(INFO, how='left')
    tb_final[['seqid', 'start', 'end', 'ProductSize2', 'Forward', 'Forward_position2',
              'Forward_length', 'Forward_TM', 'Reverse', 'Reverse_position2', 'Reverse_length', 'Reverse_TM',
              'count']].to_csv(args.output, sep='\t', index=False)
    #del_dir(args.tmp)


def parseArgs():
    parser = argparse.ArgumentParser(description='Design Primer using Primer3 and validate the specificity using ePCR')
    parser.add_argument('-i', '--info', required=True,
                        help='<file path>  A tab-separated file, containing five column to show the target region: '
                             'seqid/start/end/min product size/max product size')
    parser.add_argument('-s', '--seq', dest='fasta', required=True,
                        help='<file path> seq path (in fasta format)')
    parser.add_argument('-c', '--config',
                        help='<file path> specify a Primer3 config (product sizes specified here make no sense)')
    parser.add_argument('-a', '--add',
                        help='<file path> Other restrictions add to *.p3in file. Write in Primer3 config format')
    parser.add_argument('-e', '--epcr', default='re-PCR',
                        help='<file path> The re-PCR bin path (if it not in PATH)')
    parser.add_argument('-d', '--database', dest='fahash',
                        help='<file path> fahash database')
    parser.add_argument('--primer3', default='primer3_core',
                        help='<file path> The primer3_core bin path (if it not in PATH)')
    parser.add_argument('--tmp', default='',
                        help='<directory path> The temporary directory Default: system temporary directory')
    parser.add_argument('-o', '--output', required=True,
                        help='<file path>  The output path')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main(parseArgs())
