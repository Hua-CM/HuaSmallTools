# -*- coding: utf-8 -*-
# @Time    : 2021/8/23 13:09
# @Author  : Zhongyi Hua

# @File    : validate_genome.py
# @Usage   :
# @Note    :
# @E-mail  : njbxhzy@hotmail.com

import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation
from Bio.Data.CodonTable import TranslationError
from collections import defaultdict


def validate(_gff, _genome):
    _seqs = SeqIO.to_dict(SeqIO.parse(_genome, 'fasta'))
    _gff = [_ for _ in GFF.parse(_gff, base_dict=_seqs)]
    correct_list = []
    error_dict = defaultdict(list)
    for scaffold in _gff:
        correct_scaffold = SeqRecord(seq="1234", id=scaffold.id, name=scaffold.name, description=scaffold.description)
        for gene in scaffold.features:
            try:
                cds_features = [_.location.extract(scaffold.seq)
                                for _ in gene.sub_features[0].sub_features
                                if _.type == "CDS"]
                if gene.strand == -1:
                    cds_features = cds_features[::-1]
                cds_seq_full_length = Seq('')
                for cds_seq in cds_features:
                    cds_seq_full_length += cds_seq
                cds_seq_full_length.translate(cds=True)
                correct_scaffold.features.append(gene)
            except TranslationError as e:
                if e.args[0].startswith("First codon"):
                    error_dict.setdefault('first', []).append(gene.id)
                elif e.args[0] == "Extra in frame stop codon found":
                    error_dict.setdefault('internal', []).append(gene.id)
                elif e.args[0].endswith("is not a stop codon"):
                    error_dict.setdefault('final', []).append(gene.id)
                elif e.args[0].endswith("is not a multiple of three"):
                    error_dict.setdefault('three', []).append(gene.id)
            except:
                print(gene.id)
        correct_list.append(correct_scaffold)
    return correct_list, error_dict


def write_result(correct_list: list, error_dict: dict, prefix):
    GFF.write(correct_list, open(prefix+'.gff3', 'w'))
    with open(prefix+'.log', 'w') as f_log:
        f_log.write('01. These genes codon_start values are invalid:\n' +
                    '\n'.join(error_dict['first']) + '\n')
        f_log.write('02. These genes has internal stop codon:\n' +
                    '\n'.join(error_dict['internal']) + '\n')
        f_log.write('03. These genes codon_end values are invalid:\n' +
                    '\n'.join(error_dict['final']) + '\n')
        f_log.write('04. These genes lengths are not a multiple of three:\n' +
                    '\n'.join(error_dict['three']) + '\n')


def getArgs():
    parser = argparse.ArgumentParser(description="Validate whether the features in gff is legal. Will generate a log \
                                                  file and a corresponding gff file")
    parser.add_argument('fasta', help="genome fasta file")
    parser.add_argument('gff', help="gff file")
    parser.add_argument('-p', '--prefix', type=str, default='validated', help="output file prefix. Default: genome")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    Args = getArgs()
    correct, err = validate(Args.gff, Args.fasta)
    write_result(correct, err, Args.prefix + '_validated')
