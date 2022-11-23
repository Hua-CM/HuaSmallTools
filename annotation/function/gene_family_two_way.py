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
        if self.seed:
            seed_proteins = [all_proteins.get(_) for _ in self.seed]
        else:  # do not need accessions, the ref fasta file is the ref sequences
            seed_proteins = all_proteins
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
            max_target_seqs=1,
            num_threads=threads,
            out=os.path.join(self.tmp_dir, 'reverse_blast.tbl')
        )
        makeblastdb_cline()
        blastp_cline()
        blast_result = pd.read_table(os.path.join(self.tmp_dir, 'reverse_blast.tbl'),
                                     header=None,
                                     names=['qacc', 'sacc', 'qlen', 'slen', 'length', 'pident', 'evalue'])
        seq_idx_lst = [_[0] for _ in enumerate(blast_result['sacc'].to_list()) if _[1] in self.seed]
        _tmp_lst = blast_result['qacc'].to_list()
        seq_lst = [_tmp_lst[_] for _ in seq_idx_lst]
        return set(seq_lst), blast_result[blast_result['qacc'].isin(seq_lst)]

    def _pfam_iden(self, threads):
        # make a temporary file
        with open(os.path.join(self.tmp_dir, 'pfam_acc'), 'w') as f_in:
            f_in.write('\n'.join(self.domain))
        # fetch PF family
        os.system(f"hmmfetch -o {os.path.join(self.tmp_dir, 'query.hmm')} -f {self.pfam} {os.path.join(self.tmp_dir, 'pfam_acc')}")
        os.system(f"hmmpress {os.path.join(self.tmp_dir, 'query.hmm')}")
        # hmmscan
        os.system(f"hmmscan --cut_ga --notextw --tblout {os.path.join(self.tmp_dir, 'pfam.tbl')} "
                  f"{os.path.join(self.tmp_dir, 'query.hmm')} "
                  f"{os.path.join(self.tmp_dir, 'putative.fasta')} > /dev/null")
        hmm_result = pd.read_table(os.path.join(self.tmp_dir, 'pfam.tbl'),
                                   comment='#',
                                   sep='\\s+',
                                   usecols=[0, 1, 2, 4],
                                   names=['domainname', 'domainacc', 'seqid', 'evalue'])
        seq_lst = hmm_result['seqid'].to_list()
        return set(seq_lst), hmm_result[hmm_result['seqid'].isin(seq_lst)]

    def identify(self, threads):
        if self.ref:  # need blast
            self._forward_blast_iden(threads)
            set1, blast_result = self._reverse_blast_iden(threads)
        else:  # do not need blast, make a file for Pfam search
            os.symlink(os.path.join(os.getcwd(), self.query), os.path.join(self.tmp_dir, 'putative.fasta'))
        if self.domain:  # need Pfam
            set2, hmm_result = self._pfam_iden(threads)
        local_vars = locals()
        if 'set1' in local_vars and 'set2' in local_vars:
            gene_set = set1 & set2
        elif 'set1' in local_vars:
            gene_set = set1
        else:
            gene_set = set2
        seq_lst = [_ for _ in SeqIO.parse(self.query, 'fasta') if _.id in gene_set]
        return_dict = {'seq_lst': seq_lst}
        for _key in ['blast_result', 'hmm_result']:
            if local_vars.get(_key) is not None:
                return_dict[_key] = local_vars.get(_key)
        return return_dict


def parse_args():
    parser = argparse.ArgumentParser(
        prog="GeneFamily",
        formatter_class=argparse.RawTextHelpFormatter,
        description=''' 
        -------------------------------------------------------------------------------------------------------
        %(prog)s 
        Author:  Zhongyi Hua
        Version: v2.0
        This is the script for bidirection gene family identification. Ensure HMMER in your PATH and tweak your 
        Pfam-A database before index using the following commands: sed -i "/^ACC/s/\.[0-9]\+//g" Pfam-A.hmm
        Example: python gene_family_two_way.py -q query.fa -d PF00931 -p /home/database/pfam  -o NB-ARC -@ 10
        --------------------------------------------------------------------------------------------------------
        ''')

    parser.add_argument('-q', '--query', required=True,
                        help='<file_path>  A fasta file containing query species protein sequences')
    parser.add_argument('-s', '--seed', default="",
                        help='<file_path>  A file containing seed accession')
    parser.add_argument('-r', '--reference', dest='ref', default="",
                        help='<file_path>  A fasta file reference species protein sequences')
    parser.add_argument('-d', '--domain', default="",
                        help='<character>  The Pfam domain accession, '
                             'if there are multiple accessions, separate by coma.')
    parser.add_argument('-p', '--pfam', default="",
                        help='<file_path>  The tweaked Pfam-A.hmm (see help, please index it first)')
    parser.add_argument('-o', '--output', required=True,
                        help='<file_path>  output file prefix')
    parser.add_argument('-@', '--threads', default=4,
                        help='<NUM> number of parallel CPU workers. Default: 4')
    parser.add_argument('-t', '--tmp',
                        help='<dir> tmp directory. If will be kept after identification to provide more information')
    args = parser.parse_args()
    return args


def main(args):
    if args.tmp:
        tmp_dir = args.tmp
    else:
        tmp_dir = tempfile.mktemp()
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    try:
        if args.seed:
            with open(args.seed) as f_in:
                seed_lst = [_.strip() for _ in f_in.read().strip().splitlines()]
        else:
            seed_lst = []
        if args.domain:
            domain_lst = args.domain.split(',')
        else:
            domain_lst = []
        gfins = GeneFamilyTwoDirect(args.query, seed_lst, args.ref, domain_lst, args.pfam, tmp_dir)
        res_dict = gfins.identify(args.threads)
        SeqIO.write(res_dict['seq_lst'], args.output + '.fasta', 'fasta')
        if res_dict.get('hmm_result') is not None:
            res_dict.get('hmm_result').to_csv(args.output + '_hmm.txt', sep='\t', index=False)
        if res_dict.get('blast_result') is not None:
            res_dict.get('blast_result').to_csv(args.output + '_blast.txt', sep='\t', index=False)
    except Exception as e:
        print(e)
    finally:
        if not args.tmp:
            del_directory(tmp_dir)

if __name__ == '__main__':
    main(parse_args())
