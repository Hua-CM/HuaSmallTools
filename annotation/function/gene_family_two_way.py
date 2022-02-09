# -*- coding: utf-8 -*-
# @Time    : 2022/2/9 14:09
# @Author  : Zhongyi Hua
# @File    : gene_family_two_way.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import pandas as pd
import os
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


class GeneFamilyTwoDirect:
    def __init__(self, query, seed, ref, domain, pfam, tmp_dir):
        self.query = query
        self.ref = ref
        self.seed = seed
        self.domain = domain
        self.pfam = pfam
        self.tmp_dir = tmp_dir

    def _forward_blast_iden(self, threads):
        all_proteins = SeqIO.to_dict(SeqIO.parse(self.ref, 'fasta'))
        seed_proteins = [all_proteins.get(_) for _ in self.seed]
        SeqIO.write(seed_proteins, os.path.join(self.tmp_dir, 'seed.fasta'), 'fasta')
        del all_proteins, seed_proteins
        os.symlink(os.path.abspath(self.query), os.path.join(self.tmp_dir, 'query.fasta'))
        makeblastdb_cline = NcbimakeblastdbCommandline(
            dbtype='prot',
            input_file=os.path.join(self.tmp_dir, 'query.fasta')
        )

        blastp_cline = NcbiblastpCommandline(
            query=os.path.join(self.tmp_dir, 'seed.fasta'),
            db=os.path.join(self.tmp_dir, 'query.fasta'),
            evalue='1e-3',
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
        seq_list = [_ for _ in SeqIO.parse(self.query, 'fasta') if _.id in blast_result['sacc'].to_list()]
        SeqIO.write(seq_list, os.path.join(self.tmp_dir, 'putative.fasta'), 'fasta')

    def _reverse_blast_iden(self, threads):
        os.symlink(os.path.abspath(self.ref), os.path.join(self.tmp_dir, 'ref.fasta'))
        makeblastdb_cline = NcbimakeblastdbCommandline(
            dbtype='prot',
            input_file=os.path.join(self.tmp_dir, 'ref.fasta')
        )
        blastp_cline = NcbiblastpCommandline(
            query=os.path.join(self.tmp_dir, 'putative.fasta'),
            db=os.path.join(self.tmp_dir, 'ref.fasta'),
            evalue='1e-3',
            outfmt="6 qacc sacc qlen slen length pident evalue",
            max_hsps=1,
            num_threads=threads,
            out=os.path.join(self.tmp_dir, 'blast2.tbl')
        )
        makeblastdb_cline()
        blastp_cline()
        blast_result = pd.read_table(os.path.join(self.tmp_dir, 'blast2.tbl'),
                                     header=None,
                                     names=['qacc', 'sacc', 'qlen', 'slen', 'length', 'pident', 'evalue'])
        seq_idx_lst = [_[0] for _ in enumerate(blast_result['sacc'].to_list()) if _[1] in self.seed]
        _tmp_lst = blast_result['qacc'].to_list()
        seq_lst = [_tmp_lst[_] for _ in seq_idx_lst]
        return set(seq_lst)

    def _pfam_iden(self, threads):
        os.system('pfam_scan.pl -fasta ' +
                  os.path.join(self.tmp_dir, 'putative.fasta') +
                  ' -cpu ' + str(threads) +
                  ' -dir ' + self.pfam +
                  ' -outfile ' + os.path.join(self.tmp_dir, 'pfam.tbl'))
        query_hmm = pd.read_table(os.path.join(self.tmp_dir, 'pfam.tbl'),
                                  comment='#',
                                  sep='\\s+',
                                  usecols=[0, 3, 4, 5, 6],
                                  names=['seqid', 'start', 'end', 'domainacc', 'domainname'])
        seq_idx_lst = [_[0] for _ in enumerate(query_hmm['domainacc'].to_list()) if _[1].split('.')[0] in self.domain]
        _tmp_lst = query_hmm['seqid'].to_list()
        seq_lst = [_tmp_lst[_] for _ in seq_idx_lst]
        return set(seq_lst)

    def identify(self, threads):
        self._forward_blast_iden(threads)
        set1 = self._reverse_blast_iden(threads)
        set2 = self._pfam_iden(threads)
        gene_set = set1 & set2
        seq_lst = [_ for _ in SeqIO.parse(self.query, 'fasta') if _.id in gene_set]
        return seq_lst


def parse_args():
    parser = argparse.ArgumentParser(
        description='This is the script for one-direction gene family identification. Ensure pfam_scan.pl in your PATH')
    parser.add_argument('-q', '--query', required=True,
                        help='<file_path>  A fasta file containing query species protein sequences')
    parser.add_argument('-s', '--seed', required=True,
                        help='<file_path>  A file containing seed accession')
    parser.add_argument('-r', '--reference', dest='ref', required=True,
                        help='<file_path>  A fasta file reference species protein sequences')
    parser.add_argument('-d', '--domain', required=True,
                        help='<character>  The Pfam domain accession, '
                             'if there are multiple accessions, separate by coma.')
    parser.add_argument('-p', '--pfam', required=True,
                        help='<dir_path>  The directory containing Pfam-A.hmm')
    parser.add_argument('-o', '--output', required=True,
                        help='<file_path>  output file path')
    parser.add_argument('-t', '--threads', default=4,
                        help='<NUM> number of parallel CPU workers. Default: 4')
    args = parser.parse_args()
    return args


def main(args):
    tmp_dir = tempfile.mktemp()
    os.mkdir(tmp_dir)
    try:
        with open(args.seed) as f_in:
            seed_lst = [_.strip() for _ in f_in.readlines()]
        domain_lst = args.domain.split(',')
        gfins = GeneFamilyTwoDirect(args.query, seed_lst, args.ref, domain_lst, args.pfam, tmp_dir)
        seq_lst = gfins.identify(args.threads)
        SeqIO.write(seq_lst, args.output, 'fasta')
    except Exception as e:
        print(e)
    finally:
        del_directory(tmp_dir)


if __name__ == '__main__':
    main(parse_args())
