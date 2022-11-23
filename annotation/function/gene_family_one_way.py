# -*- coding: utf-8 -*-
# @Time    : 2021/3/23 10:37
# @Author  : Zhongyi Hua
# @File    : gene_family_one_way.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import pandas as pd
import os
import shutil
import argparse
import tempfile
from Bio.Blast.Applications import NcbiblastpCommandline, NcbimakeblastdbCommandline
from Bio import SeqIO

pd.set_option('precision', 2)


def del_directory(_dir):
    for r_, _d_, f_ in os.walk(_dir):
        for files in f_:
            os.remove(os.path.join(r_, files))
        os.removedirs(r_)


class GeneFamily:
    def __init__(self, seed, database, tmp_dir):
        self.seed = seed
        self.db = database
        self.tmp_dir = tmp_dir

    def blast_iden(self, threads=2):
        shutil.copyfile(self.db, os.path.join(self.tmp_dir, 'database.fasta'))
        makeblastdb_cline = NcbimakeblastdbCommandline(
            dbtype='prot',
            input_file=os.path.join(self.tmp_dir, 'database.fasta')
        )
        blastp_cline = NcbiblastpCommandline(
            query=self.seed,
            db=os.path.join(self.tmp_dir, 'database.fasta'),
            evalue='1e-5',
            outfmt="6 qacc sacc qlen slen length pident evalue",
            max_hsps=1,
            num_threads=threads,
            out=os.path.join(self.tmp_dir, 'blast.tbl')
        )
        makeblastdb_cline()
        blastp_cline()
        blast_result = pd.read_table(os.path.join(self.tmp_dir, 'blast.tbl'),
                                     header=None,
                                     names=['qacc', 'sacc', 'qlen', 'slen', 'length', 'pident', 'evalue'])
        blast_result = blast_result[
            (blast_result['pident'] > 40) & (blast_result['length'] / blast_result['slen'] > 0.2) & (blast_result['length'] / blast_result['qlen'] > 0.2)]
        blast_result.to_csv(os.path.join(self.tmp_dir, 'blast2.tbl'), sep='\t', index=False)
        seq_list = [_ for _ in SeqIO.parse(self.db, 'fasta') if _.id in blast_result['sacc'].to_list()]
        SeqIO.write(seq_list, os.path.join(self.tmp_dir, 'subgenes.fasta'), 'fasta')

    @staticmethod
    def test_list(lst1, lst2):
        start = 0
        try:
            for item in lst1:
                start = lst2.index(item, start) + 1
        except ValueError:
            return False
        return True

    @staticmethod
    def _cal_level_(q_dom, s_dom):
        if q_dom == s_dom == []:
            return 'C'
        if q_dom == s_dom:
            return 'A'
        if GeneFamily.test_list(s_dom, q_dom):
            return 'B'
        elif len(set(s_dom).intersection(set(q_dom))) / len(q_dom) >= 0.5:
            return 'D'
        else:
            return 'E'

    def pfam_iden(self, _pfam_db, threads=2):
        os.system('pfam_scan.pl -fasta ' +
                  self.seed +
                  ' -cpu ' + str(threads) +
                  ' -dir ' + _pfam_db +
                  ' -outfile ' + os.path.join(self.tmp_dir, 'query.tbl'))
        os.system('pfam_scan.pl -fasta ' +
                  os.path.join(self.tmp_dir, 'subgenes.fasta') +
                  ' -cpu ' + str(threads) +
                  ' -dir ' + _pfam_db +
                  ' -outfile ' + os.path.join(self.tmp_dir, 'sub.tbl'))
        blast_result = pd.read_table(os.path.join(self.tmp_dir, 'blast2.tbl'))
        query_hmm = pd.read_table(os.path.join(self.tmp_dir, 'query.tbl'),
                                  comment='#',
                                  sep='\\s+',
                                  usecols=[0, 3, 4, 5, 6],
                                  names=['seqid', 'start', 'end', 'domainacc', 'domainname'])
        subject_hmm = pd.read_table(os.path.join(self.tmp_dir, 'sub.tbl'),
                                    comment='#',
                                    sep='\\s+',
                                    usecols=[0, 3, 4, 5, 6],
                                    names=['seqid', 'start', 'end', 'domainacc', 'domainname'])
        query_hmm['temp'] = query_hmm['domainacc'] + ':' + query_hmm['domainname']
        subject_hmm['temp'] = subject_hmm['domainacc'] + ':' + subject_hmm['domainname']
        query_dict = query_hmm.sort_values(by=['seqid', 'start']). \
            groupby('seqid'). \
            agg({'temp': list}). \
            to_dict()['temp']
        subject_dict = subject_hmm.sort_values(by=['seqid', 'start']). \
            groupby('seqid'). \
            agg({'temp': list}). \
            to_dict()['temp']
        new_list = []
        for _idx, _row in blast_result.iterrows():
            new_list.append({
                'query': _row['qacc'],
                'subject': _row['sacc'],
                'Qoverlap': round(float(_row['length']) / float(_row['qlen']),2),
                'Soverlap': round(float(_row['length']) / float(_row['slen']),2),
                'identity': _row['pident'],
                'E-value': _row['evalue'],
                'QPFAM': ' '.join(query_dict.get(_row['qacc'], [])),
                'SPFAM': ' '.join(subject_dict.get(_row['sacc'], [])),
                'level': GeneFamily._cal_level_(query_dict.get(_row['qacc'], []), subject_dict.get(_row['sacc'], []))
            })
        return pd.DataFrame(new_list)


def parse_args():
    parser = argparse.ArgumentParser(
        description='This is the script for one-direction gene family identification')
    parser.add_argument('-s', '--seed', required=True,
                        help='<file_path>  A fasta file containing seed')
    parser.add_argument('-d', '--database', required=True,
                        help='<file_path>  A fasta file containing subject sequences')
    parser.add_argument('-p', '--pfam', required=True,
                        help='<dir_path>  The directory containing Pfam-A.hmm')
    parser.add_argument('-o', '--output', required=True,
                        help='<file_path>  output file path')
    parser.add_argument('-@', '--threads', default=2,
                        help='<NUM> number of parallel CPU workers (default:2)')
    parser.add_argument('-t', '--tmp', default='',
                        help='<Path> The temporary directory path')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    if args.tmp:
        tmp_dir = args.tmp
    else:
        tmp_dir = tempfile.mktemp()
    os.mkdir(tmp_dir)
    ins = GeneFamily(args.seed, args.database, tmp_dir)
    ins.blast_iden(args.threads)
    result_df = ins.pfam_iden(args.pfam, args.threads)
    result_df.to_csv(args.output, sep='\t', index=False)
    #del_directory(tmp_dir)


if __name__ == '__main__':
    main()
