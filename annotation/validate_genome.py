# -*- coding: utf-8 -*-
# @Time    : 2021/8/23 13:09
# @Author  : Zhongyi Hua

# @File    : validate_genome.py
# @Usage   :
# @Note    :
# @E-mail  : njbxhzy@hotmail.com

import argparse
import copy
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
from Bio.Data.CodonTable import TranslationError
from collections import defaultdict

start_codon = {'ATG'}
stop_codon = {'TAG', 'TAA', 'TGA'}


def get_cds(gene, scaffold_seq):
    cds_features = [_ for _ in gene.sub_features[0].sub_features if _.type == "CDS"]
    cds_features.sort(key=lambda x: x.location.start)
    if gene.strand == -1:
        cds_features = cds_features[::-1]
    if not cds_features[0].qualifiers['phase'] == ['0']:
        raise TranslationError('The phase of first CDS is not 0')
    cds_features = [_.location.extract(scaffold_seq.seq) for _ in cds_features]
    cds_seq_full_length = Seq('')
    for cds_seq in cds_features:
        cds_seq_full_length += cds_seq
    protein_seq = cds_seq_full_length.translate(cds=True)
    return protein_seq


def to_gene_dict(_gff: list):
    gene_dict = {}
    for scaffold in _gff:
        for gene in scaffold.features:
            gene_dict[gene.id] = gene
    return gene_dict


def correct_stop_codon(gene, scaffold_seq):
    _gene = copy.deepcopy(gene)
    cds_features = [_ for _ in _gene.sub_features[0].sub_features if _.type == "CDS"]
    cds_features.sort(key=lambda x: x.location.start)
    cds_seq_full_length = Seq('')
    if _gene.strand == 1:
        for _feature in cds_features[:-1]:
            cds_seq_full_length += _feature.location.extract(scaffold_seq.seq)
        # extension
        while True:
            cds_features[-1].location = FeatureLocation(cds_features[-1].location.start,
                                                        cds_features[-1].location.end + 3,
                                                        strand=1)
            cds_seq_test = cds_seq_full_length + cds_features[-1].location.extract(scaffold_seq.seq)
            try:
                cds_seq_test.translate(cds=True)
                break
            except TranslationError:
                continue
        # modify_gene
        new_sub_features = []
        for _ in _gene.sub_features[0].sub_features:
            if _.type == "exon" and _.location.start == cds_features[-1].location.start:
                _.location = FeatureLocation(_.location.start, cds_features[-1].location.end, strand=1)
            if _.type == 'three_prime_UTR':
                continue
            new_sub_features.append(_)
        _gene.location = FeatureLocation(_gene.location.start, cds_features[-1].location.end, strand=1)
        _gene.sub_features[0].location = FeatureLocation(_gene.sub_features[0].location.start,
                                                         cds_features[-1].location.end,
                                                         strand=1)
    else:
        cds_features = cds_features[::-1]
        for _feature in cds_features[:-1]:
            cds_seq_full_length += _feature.location.extract(scaffold_seq.seq)
        while True:
            cds_features[-1].location = FeatureLocation(cds_features[-1].location.start - 3,
                                                        cds_features[-1].location.end,
                                                        strand=-1)
            cds_seq_test = cds_seq_full_length + cds_features[-1].location.extract(scaffold_seq.seq)
            try:
                cds_seq_test.translate(cds=True)
                break
            except TranslationError:
                continue
        # modify_gene
        new_sub_features = []
        for _ in _gene.sub_features[0].sub_features:
            if _.type == 'three_prime_UTR':
                continue
            if _.type == "exon" and _.location.start == cds_features[-1].location.start:
                _.location = FeatureLocation(cds_features[-1].location.start, _.location.end, strand=-1)
            new_sub_features.append(_)
        _gene.location = FeatureLocation(cds_features[-1].location.start, _gene.location.end, strand=-1)
        _gene.sub_features[0].location = FeatureLocation(cds_features[-1].location.start,
                                                         _gene.sub_features[0].location.end,
                                                         strand=-1)
    _gene.sub_features[0].sub_features = new_sub_features
    return _gene


def correct_start_codon(gene, scaffold_seq):
    _gene = copy.deepcopy(gene)
    cds_features = [_ for _ in _gene.sub_features[0].sub_features if _.type == "CDS"]
    cds_features.sort(key=lambda x: x.location.start)
    if _gene.strand == 1:
        while True:
            _start_pos = cds_features[0].location.start
            cds_features[0].location = FeatureLocation(_start_pos-3, cds_features[0].location.end, strand=1)
            _codon = scaffold_seq.seq[_start_pos-3: _start_pos]
            if _codon in start_codon:
                break
            if _codon in stop_codon:
                raise TranslationError('First codon could not be found')
        # modify_gene
        new_sub_features = []
        for _ in _gene.sub_features[0].sub_features:
            if _.type == "exon" and _.location.end == cds_features[0].location.end:
                _.location = cds_features[0].location
            if _.type == 'five_prime_UTR':
                continue
            new_sub_features.append(_)
        _gene.location = FeatureLocation(cds_features[0].location.start, _gene.location.end, strand=1)
        _gene.sub_features[0].location = FeatureLocation(cds_features[0].location.start,
                                                         _gene.sub_features[0].location.end,
                                                         strand=1)
    else:
        while True:
            _start_pos = cds_features[-1].location.end
            cds_features[-1].location = FeatureLocation(cds_features[-1].location.start, _start_pos + 3, strand=-1)
            _codon = scaffold_seq.seq[_start_pos - 3: _start_pos].reverse_complement()
            if _codon in start_codon:
                break
            if _codon in stop_codon:
                raise TranslationError('First codon could not be found')
        # modify_gene
        new_sub_features = []
        for _ in _gene.sub_features[0].sub_features:
            if _.type == "exon" and _.location.start == cds_features[-1].location.start:
                _.location = cds_features[-1].location
            if _.type == 'five_prime_UTR':
                continue
            new_sub_features.append(_)
        _gene.location = FeatureLocation(_gene.location.start, cds_features[-1].location.end, strand=-1)
        _gene.sub_features[0].location = FeatureLocation(_gene.sub_features[0].location.start,
                                                         cds_features[-1].location.end,
                                                         strand=-1)
    _gene.sub_features[0].sub_features = new_sub_features
    return _gene


def correct_phase(gene, scaffold_seq):
    _gene = copy.deepcopy(gene)
    cds_features = [_ for _ in _gene.sub_features[0].sub_features if _.type == "CDS"]
    cds_features.sort(key=lambda x: x.location.start)
    if _gene.strand == 1:
        if cds_features[0].qualifiers['phase'] == ['1']:
            cds_features[0].location = FeatureLocation(cds_features[0].location.start - 1,
                                                       cds_features[0].location.end,
                                                       strand=1)
        elif cds_features[0].qualifiers['phase'] == ['2']:
            cds_features[0].location = FeatureLocation(cds_features[0].location.start - 2,
                                                       cds_features[0].location.end,
                                                       strand=1)
        else:
            raise TranslationError('The phase was not 0/1/2')
    else:
        if cds_features[-1].qualifiers['phase'] == ['1']:
            cds_features[-1].location = FeatureLocation(cds_features[-1].location.start,
                                                        cds_features[-1].location.end + 1,
                                                        strand=-1)
        elif cds_features[-1].qualifiers['phase'] == ['2']:
            cds_features[-1].location = FeatureLocation(cds_features[-1].location.start,
                                                        cds_features[-1].location.end + 2,
                                                        strand=-1)
        else:
            raise TranslationError('The phase was not 0/1/2')
        cds_features = cds_features[::-1]
    cds_seqs = [_.location.extract(scaffold_seq.seq) for _ in cds_features]
    cds_seq_full_length = Seq('')
    for cds_seq in cds_seqs:
        cds_seq_full_length += cds_seq
    # modify_gene
    cds_features[0].qualifiers['phase'] = ['0']
    new_sub_features = []
    if _gene.strand == 1:
        for _ in _gene.sub_features[0].sub_features:
            if _.type == "exon" and _.location.end == cds_features[0].location.end:
                _.location = cds_features[0].location
            if _.type == 'five_prime_UTR':
                continue
            new_sub_features.append(_)
        _gene.location = FeatureLocation(cds_features[0].location.start, _gene.location.end, strand=1)
        _gene.sub_features[0].location = FeatureLocation(cds_features[0].location.start,
                                                         _gene.sub_features[0].location.end,
                                                         strand=1)
    else:
        for _ in _gene.sub_features[0].sub_features:
            if _.type == "exon" and _.location.start == cds_features[0].location.start:
                _.location = cds_features[0].location
            if _.type == 'five_prime_UTR':
                continue
            new_sub_features.append(_)
        _gene.location = FeatureLocation(_gene.location.start, cds_features[0].location.end, strand=-1)
        _gene.sub_features[0].location = FeatureLocation(_gene.sub_features[0].location.start,
                                                         cds_features[0].location.end,
                                                         strand=-1)
    _gene.sub_features[0].sub_features = new_sub_features
    try:
        cds_seq_full_length.translate(cds=True)
        return _gene
    except TranslationError:
        raise TranslationError('These genes need another round correction', _gene)


def correct(_gff, _genome):
    _seqs = SeqIO.to_dict(SeqIO.parse(_genome, 'fasta'))
    _gff = [_ for _ in GFF.parse(_gff, base_dict=_seqs)]
    correct_list = []
    error_dict = defaultdict(list)
    for scaffold in _gff:
        correct_scaffold = SeqRecord(seq="1234", id=scaffold.id, name=scaffold.name, description=scaffold.description)
        for gene in scaffold.features:
            try:
                get_cds(gene, scaffold)
                correct_scaffold.features.append(gene)
            except TranslationError as e:
                try:
                    if e.args[0].startswith("First codon"):
                        correct_scaffold.features.append(correct_start_codon(gene, scaffold))
                        error_dict.setdefault('corrected', []).append(gene.id)
                    elif e.args[0].startswith('The phase of first CDS is not 0'):
                        correct_scaffold.features.append(correct_phase(gene, scaffold))
                        error_dict.setdefault('corrected', []).append(gene.id)
                    elif e.args[0].endswith("is not a stop codon"):
                        correct_scaffold.features.append(correct_stop_codon(gene, scaffold))
                        error_dict.setdefault('corrected', []).append(gene.id)
                    # can not handle for now
                    elif e.args[0] == "Extra in frame stop codon found":
                        error_dict.setdefault('internal', []).append(gene.id)
                    elif e.args[0].endswith("is not a multiple of three"):
                        error_dict.setdefault('three', []).append(gene.id)
                except TranslationError as e2:
                    if e2.args[0].startswith('These genes need another round correction'):
                        correct_scaffold.features.append(e2.args[1])
                        error_dict.setdefault('phase', []).append(gene.id)
                    # for second round
                    elif e.args[0] == "Extra in frame stop codon found":
                        error_dict.setdefault('internal', []).append(gene.id)
                    elif e2.args[0].startswith('First codon could not be found'):
                        error_dict.setdefault('first2', []).append(gene.id)
                    elif e.args[0].endswith("is not a stop codon"):
                        error_dict.setdefault('final', []).append(gene.id)
                    elif e.args[0].endswith("is not a multiple of three"):
                        error_dict.setdefault('three', []).append(gene.id)
            except Exception as e:
                print(e)
                print(gene.id)
        correct_list.append(correct_scaffold)
    return correct_list, error_dict


def write_result(correct_list: list, error_dict: dict, prefix):
    GFF.write(correct_list, open(prefix+'.gff3', 'w'))
    with open(prefix+'.log', 'w') as f_log:
        f_log.write('01. These genes start codon are illegal (not ATG):\n' +
                    '\n'.join(error_dict['first']) + '\n')
        f_log.write('02. These genes start codon are illegal (not ATG) and could not be corrected:\n' +
                    '\n'.join(error_dict['first2']) + '\n')
        f_log.write('03. These genes\'s first phase is not 0 and need a second round correction:\n' +
                    '\n'.join(error_dict['phase']) + '\n')
        f_log.write('04. These genes has internal stop codon:\n' +
                    '\n'.join(error_dict['internal']) + '\n')
        f_log.write('05. These genes end codon are illegal:\n' +
                    '\n'.join(error_dict['final']) + '\n')
        f_log.write('06. These genes lengths are not a multiple of three:\n' +
                    '\n'.join(error_dict['three']) + '\n')
        f_log.write('07. These genes corrected:\n' +
                    '\n'.join(error_dict['corrected']) + '\n')


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
    correct, err = correct(Args.gff, Args.fasta)
    write_result(correct, err, Args.prefix + '_validated')
