# -*- coding: utf-8 -*-
# @Time    : 2022/2/21 21:24
# @Author  : Zhongyi Hua
# @File    : primer3_ASPCR.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import argparse
import os
import pandas as pd
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tempfile import mktemp
from primer3_with_ePCR_validation import del_dir, get_seq


class Primer3SNP:
    c_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # s_dict = {'A': 'G', 'T': 'C', 'C': 'T', 'G': 'A'}
    m_dict = {'A': 'C', 'T': 'G', 'C': 'A', 'G': 'T'}
    w_dict = {'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G'}
    mismatch_dict = {frozenset(('A', 'G')): 's', frozenset(('C', 'T')): 's',
                     frozenset(('C', 'A')): 'm', frozenset(('G', 'T')): 'm',
                     frozenset(('A', 'A')): 'w', frozenset(('T', 'T')): 'w',
                     frozenset(('C', 'C')): 'w', frozenset(('G', 'G')): 'w'}

    def __init__(self, _infodf, _tmp, _primer_bin):
        """
        :param _infodf: pd.Dataframe contain AT LEAST five columns: 'seqid', 'pos', 'ref', 'alt'
        :param _tmp:
        :param _primer_bin:
        """
        self.infodf = _infodf
        self.bin = _primer_bin
        self.p3in = mktemp(dir=_tmp)
        self.p3in2 = mktemp(dir=_tmp)  # use SNP as the right primer
        self.p3out = mktemp(dir=_tmp)
        self.p3out2 = mktemp(dir=_tmp)  # use SNP as the right primer

    @staticmethod
    def readfile(_file):
        with open(_file) as f_in:
            return f_in.read().split('\n')

    def _import_mismatch(self, _primer: str, mut: str, _type) -> str:
        if _type == 'left':
            if self.mismatch_dict.get(frozenset((_primer[-1], mut))) == 's':
                mismatch_dict = self.w_dict
            elif self.mismatch_dict.get(frozenset((_primer[-1], mut))) == 'm':
                mismatch_dict = self.m_dict
            _primer2 = list(_primer)
            _primer2[-3] = mismatch_dict.get(self.c_dict.get(_primer2[-3]))
        else:
            if self.mismatch_dict.get(frozenset((_primer[0], mut))) == 's':
                mismatch_dict = self.w_dict
            elif self.mismatch_dict.get(frozenset((_primer[0], mut))) == 'm':
                mismatch_dict = self.m_dict
            _primer2 = list(_primer)
            _primer2[2] = mismatch_dict.get(self.c_dict.get(_primer2[2]))
        return ''.join(_primer2)

    def prepare_p3in(self, _fasta):
        # _fasta: a dict from SeqIO.to_dict
        _write_list = []
        for _item in self.infodf.itertuples():
            _seq_ref = get_seq(_fasta[_item.seqid], _item.pos, _item.pos, 500+100)
            # prepare mutation seq
            _seq_lst = list(_seq_ref)
            _seq_lst[500+100] = _item.alt
            _seq_alt = SeqRecord(''.join(_seq_lst), id=_seq_ref.id)
            # Ref record
            _sequence_ref = 'SEQUENCE_TEMPLATE=' + str(_seq_ref.seq)
            # Alt record
            _sequence_alt = 'SEQUENCE_TEMPLATE=' + str(_seq_alt.seq)
            for _suffix, _suffix2, _seq in [(_item.ref, _item.alt, _sequence_ref), (_item.alt, _item.ref, _sequence_alt)]:
                _write_list.append('SEQUENCE_ID=' + _seq_ref.id + '_' + _suffix + '_' + _suffix2)
                _write_list.append(_seq)
                # IMPORTANCE!: It seems that Primer3 using the computational index! so 601 is 600
                _write_list.append('SEQUENCE_FORCE_LEFT_END=600')
                _write_list.append('PRIMER_PRODUCT_SIZE_RANGE=300-500')
                _write_list.append('=')
        with open(self.p3in, 'w') as f2:
            f2.write('\n'.join(_write_list))
            f2.write('\n')
        _write_list2 = []
        for _record in _write_list:
            if _record == 'SEQUENCE_FORCE_LEFT_END=600':
                _write_list2.append('SEQUENCE_FORCE_RIGHT_START=600')
            else:
                _write_list2.append(_record)
        with open(self.p3in2, 'w') as f2:
            f2.write('\n'.join(_write_list2))
            f2.write('\n')

    def run_primer3(self, _config):
        if _config:
            sp.run([self.bin, '--p3_settings_file=%s' % _config, '--output=%s' % self.p3out, self.p3in])
            sp.run([self.bin, '--p3_settings_file=%s' % _config, '--output=%s' % self.p3out2, self.p3in2])
        else:
            sp.run([self.bin,  '--output=%s' % self.p3out, self.p3in])
            sp.run([self.bin, '--output=%s' % self.p3out2, self.p3in2])

    def parse_p3out(self):
        tb_p3out = pd.DataFrame()
        for _type, p3out in [('left', self.p3out), ('right', self.p3out2)]:
            _file = self.readfile(p3out)
            _module = []
            _record_list = []
            for _line in _file:
                _module.append(_line)
                if _line == '=':
                    _model_dict = dict(_.split('=') for _ in _module)
                    if 'SEQUENCE_ID' in _model_dict:
                        if not _model_dict['PRIMER_PAIR_NUM_RETURNED'] == '0':
                            for _num in range(int(_model_dict['PRIMER_PAIR_NUM_RETURNED'])):
                                _chr, _genome_pos, _, _dectect, _mutant = _model_dict['SEQUENCE_ID'].split('_')
                                _record = {'ID': _chr,
                                           'GenomePosition': _genome_pos,
                                           'Detect': _dectect,
                                           'Mutant': _mutant,
                                           'ProductSize': _model_dict['PRIMER_PAIR_{}_PRODUCT_SIZE'.format(_num)],
                                           'Forward': self._import_mismatch(_model_dict['PRIMER_LEFT_{}_SEQUENCE'.format(_num)], _mutant, _type) if _type == 'left' else _model_dict['PRIMER_LEFT_{}_SEQUENCE'.format(_num)],
                                           'ForwardPosition': _model_dict['PRIMER_LEFT_{}'.format(_num)].split(',')[0],
                                           'ForwardLength': _model_dict['PRIMER_LEFT_{}'.format(_num)].split(',')[1],
                                           'ForwardTM': _model_dict['PRIMER_LEFT_{}_TM'.format(_num)],
                                           'Reverse': _model_dict['PRIMER_RIGHT_{}_SEQUENCE'.format(_num)] if _type == 'left' else self._import_mismatch(_model_dict['PRIMER_RIGHT_{}_SEQUENCE'.format(_num)], _mutant, _type),
                                           'ReversePosition': _model_dict['PRIMER_RIGHT_{}'.format(_num)].split(',')[0],
                                           'ReverseLength': _model_dict['PRIMER_RIGHT_{}'.format(_num)].split(',')[1],
                                           'ReverseTM': _model_dict['PRIMER_RIGHT_{}_TM'.format(_num)],
                                           }
                                _record_list.append(_record)
                else:
                    continue
            tb_p3out = tb_p3out.append(pd.DataFrame(_record_list))
        return tb_p3out


def main(args):
    FASTA = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))
    INFO = pd.read_table(args.info,
                         header=None,
                         names=['seqid', 'pos', 'ref', 'alt', 'exp'],
                         dtype={'seqid': 'str', 'pos': 'int64', 'ref': 'str', 'alt': 'str', 'exp': 'str'})
    if (not args.tmp == '') and (not os.path.exists(args.tmp)):
        os.mkdir(args.tmp)
    elif args.tmp == '':
        args.tmp = mktemp()
    # main
    primer3 = Primer3SNP(INFO, args.tmp, args.primer3)
    primer3.prepare_p3in(FASTA)
    primer3.run_primer3(args.config)
    tb_p3out = primer3.parse_p3out()
    INFO.rename(columns={'seqid': 'ID', 'pos': 'GenomePOsition', 'ref': 'Ref', 'alt': 'Alt'}, inplace=True)
    tb_p3out.merge(INFO, how='left').to_csv(args.output, sep='\t', index=False)
    del_dir(args.tmp)


def parseArgs():
    parser = argparse.ArgumentParser(description='Design Primer using Primer3 for AS-PCR')
    parser.add_argument('-i', '--info', required=True,
                        help='<file path>  A tab-separated file, containing at least four column to show the target region: '
                             'seqid/position/ref/alt/....(other notes, for example, "this is wilt type")')
    parser.add_argument('-s', '--seq', dest='fasta', required=True,
                        help='<file path> seq path (in fasta format)')
    parser.add_argument('--primer3', default='primer3_core',
                        help='<file path> The primer3_core bin path (if it not in PATH)')
    parser.add_argument('--config',
                        help='<file path> Primer3 config path')
    parser.add_argument('--tmp', default='',
                        help='<directory path> The temporary directory Default: system temporary directory')
    parser.add_argument('-o', '--output', required=True,
                        help='<file path>  The output path')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main(parseArgs())
