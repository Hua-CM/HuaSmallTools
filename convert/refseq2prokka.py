# -*- coding: utf-8 -*-
# @Time    : 2022/4/20 15:23
# @Author  : Zhongyi Hua
# @File    : refseq2prokka.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import argparse
from BCBio import GFF
from Bio import SeqIO


def convert(gff_path, genome_path, output_path):
    seq_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
    chrs = [_ for _ in GFF.parse(gff_path, base_dict=seq_dict)]
    for chromosome in chrs:
        chromosome.annotations = {}
        for gene in chromosome.features:
            if gene.type not in ['gene', 'tRNA', 'mRNA', 'tmRNA', 'CDS', 'pseudogene']:
                continue
            gene.qualifiers['ID'][0] = gene.qualifiers['locus_tag'][0] + "_gene"
            gene.id = gene.qualifiers['ID'][0]
            for _idx, cds in enumerate(gene.sub_features):
                if cds.type == 'CDS':
                    cds.qualifiers['ID'][0] = cds.qualifiers['locus_tag'][0]
                    cds.qualifiers['Parent'][0] = gene.id
                    cds.id = cds.qualifiers['ID'][0]
                if _idx > 0:
                    cds.qualifiers['ID'][0] = cds.qualifiers['locus_tag'][0] + "_" + str(_idx)
                    cds.id = cds.qualifiers['ID'][0]

    GFF.write(chrs, open(output_path, 'w'), include_fasta=True)


def parse_args():
    parser = argparse.ArgumentParser(
        prog="RefSeq2Prokka",
        formatter_class=argparse.RawTextHelpFormatter,
        description=''' 
        -------------------------------------------------------------------------------------------------------
        %(prog)s 
        Author:  Zhongyi Hua
        Convert RefSeq gff and seq to Prokka format. Tidy name at the same time.
        --------------------------------------------------------------------------------------------------------
        ''')
    parser.add_argument('-a', '--anno', required=True,
                        help='<file_path>  A fasta file containing query species protein sequences')
    parser.add_argument('-s', '--seq', default="",
                        help='<file_path>  A file containing seed accession')
    parser.add_argument('-o', '--output', required=True,
                        help='<file_path>  output gff file')
    args = parser.parse_args()
    return args


def main(args):
    convert(args.anno, args.seq, args.output)


if __name__ == '__main__':
    main(parse_args())
