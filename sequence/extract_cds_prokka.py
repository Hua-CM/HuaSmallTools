# -*- coding: utf-8 -*-
# @Time : 2020/3/14 22:51
# @Author : Zhongyi Hua
# @FileName: extract_cds_prokka.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gffutils


def extract_seq(seq, feature, use_strand=True):
    seq = seq.seq[feature.start - 1:feature.end]
    if use_strand:
        if feature.strand in ["+", "."]:
            return seq
        elif feature.strand == "-":
            return seq.reverse_complement()


def extract_cds(seq_path, gff_path, out_path):
    cds_seq = []
    genome_seq = {x.id: x for x in SeqIO.parse(seq_path, "fasta")}
    genome_gff = gffutils.create_db(gff_path, ':memory:', merge_strategy="create_unique", keep_order=True)
    for ele_cds in genome_gff.features_of_type(featuretype="CDS", order_by="start"):
            cds_seq.append(SeqRecord(id=ele_cds.id,
                                     seq=extract_seq(genome_seq[ele_cds.seqid], ele_cds),
                                     description=""))
    SeqIO.write(cds_seq, out_path, "fasta")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This is the script for extracting cds from PROKKA result. Since Bio \
                                                  could not parse its gbk result correctly")
    parser.add_argument('-i', '--seq_file', required=True,
                        help='<file_path> The genome sequence file in fasta format')
    parser.add_argument('-g', '--gff_file', required=True,
                        help='<file_path> The genome annotation file in gff format')
    parser.add_argument('-o', '--cds_file', required=True,
                        help='<file_path> The cds sequence file in fasta format')
    args = parser.parse_args()
    extract_cds(args.seq_file, args.gff_file, args.cds_file)
