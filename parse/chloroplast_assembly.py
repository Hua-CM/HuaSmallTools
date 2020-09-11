# -*- coding: utf-8 -*-
# @Time : 2020/5/1 17:40
# @Author : Zhongyi Hua
# @FileName: chloroplast_assembly.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


def extract_seq(seq_dict, contig, start, end):
    if start < end:
        seq = seq_dict[contig].seq[start-1:end]
    else:
        seq = seq_dict[contig].seq[end-1:start][::-1]
    return seq


def merge_contigs(mummer_df, seq_dictG):
    result_seq = ""
    last_end = 0
    for i in range(len(mummer_df)):
        contig_seq = extract_seq(seq_dictG, mummer_df.iloc[i, 10], mummer_df.iloc[i, 2], mummer_df.iloc[i, 3])
        if mummer_df.iloc[i, 0] > last_end:
            contact_seq = "N" * (mummer_df.iloc[i, 0]-last_end-1)
            result_seq = result_seq + contact_seq + contig_seq
        elif mummer_df.iloc[i, 0] <= last_end:
            contact_seq = contig_seq[last_end-mummer_df.iloc[i, 0]+1:]
            result_seq = result_seq + contact_seq
        last_end = mummer_df.iloc[i, 1]
    return SeqRecord(id="putative_genome",
                     seq=result_seq,
                     description="")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This is the script assists chloroplast assembly")
    parser.add_argument('-c', '--contigs_fasta', required=True,
                        help='<file_path> The contigs file in fasta format')
    parser.add_argument('-m', '--mummer', required=True,
                        help='<file_path> The nucmer result *.1coords')
    parser.add_argument('-o', '--putative_genome', required=True,
                        help='<file_path> The result putative genomes with gaps in fasta format')
    args = parser.parse_args()
    seqdict = {}
    for x in SeqIO.parse(args.contigs_fasta, "fasta"):
        seqdict[x.id] = x
    mummer = pd.read_csv(args.mummer, sep="\t", header=None)
    putative_genome = merge_contigs(mummer, seqdict)
    SeqIO.write(putative_genome, args.putative_genome, "fasta")
