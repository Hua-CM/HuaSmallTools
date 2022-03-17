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
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Data.CodonTable import TranslationError
from collections import defaultdict

start_codon = {'ATG'}
stop_codon = {'TAG', 'TAA', 'TGA'}


def get_cds(mRNA, scaffold_seq):
    cds_features = [_ for _ in mRNA.sub_features if _.type == "CDS"]
    cds_features.sort(key=lambda x: x.location.start)
    if mRNA.strand == -1:
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


def correct_stop_codon(mRNA, scaffold_seq):
    _mRNA = copy.deepcopy(mRNA)
    cds_features = [_ for _ in _mRNA.sub_features if _.type == "CDS"]
    cds_features.sort(key=lambda x: x.location.start)
    cds_seq_full_length = Seq('')
    extention_len = 0
    if _mRNA.strand == 1:
        for _feature in cds_features[:-1]:
            cds_seq_full_length += _feature.location.extract(scaffold_seq.seq)
        # extension
        while True:
            extention_len += 3
            cds_features[-1].location = FeatureLocation(cds_features[-1].location.start,
                                                        cds_features[-1].location.end + 3,
                                                        strand=1)
            cds_seq_test = cds_seq_full_length + cds_features[-1].location.extract(scaffold_seq.seq)
            try:
                cds_seq_test.translate(cds=True)
                break
            except TranslationError:
                if extention_len > 3000:
                    break
                else:
                    continue
        # modify_gene
        new_sub_features = []
        utr_features = []
        for _ in _mRNA.sub_features:
            if _.type == "exon" and _.location.start == cds_features[-1].location.start:
                if _.location.end > cds_features[-1].location.end:
                    # with UTR and UTR modified, but not extend exon boundary
                    continue
                else:
                    # no UTR or extend UTR boundary
                    _.location = FeatureLocation(_.location.start, cds_features[-1].location.end, strand=1)

            if _.type == 'three_prime_UTR':
                utr_features.append(_)
                continue
            new_sub_features.append(_)
        # handle 3' UTR and corresponding exon
        for utr_3 in utr_features:
            if utr_3.location.start < cds_features[-1].location.end:
                if utr_3.location.end < cds_features[-1].location.end:
                    continue
            else:
                # handle corresponding exon
                for _ in new_sub_features:
                    if _.type == "exon" and _.location == utr_3.location:
                        _.location = FeatureLocation(cds_features[-1].location.end + 1, utr_3.location.end, strand=1)
                utr_3.location = FeatureLocation(cds_features[-1].location.end + 1, utr_3.location.end, strand=1)
                new_sub_features.append(utr_3)
    else:
        cds_features = cds_features[::-1]
        for _feature in cds_features[:-1]:
            cds_seq_full_length += _feature.location.extract(scaffold_seq.seq)
        while True:
            extention_len += 3
            cds_features[-1].location = FeatureLocation(cds_features[-1].location.start - 3,
                                                        cds_features[-1].location.end,
                                                        strand=-1)
            cds_seq_test = cds_seq_full_length + cds_features[-1].location.extract(scaffold_seq.seq)
            try:
                cds_seq_test.translate(cds=True)
                break
            except TranslationError:
                if extention_len > 3000:
                    break
                else:
                    continue
        # modify_gene
        new_sub_features = []
        utr_features = []
        for _ in _mRNA.sub_features:
            if _.type == "exon" and _.location.end == cds_features[-1].location.end:
                if _.location.start < cds_features[-1].location.start:
                    # with UTR and UTR modified, but not extend exon boundary
                    continue
                else:
                    # no UTR or extend UTR boundary
                    _.location = FeatureLocation(cds_features[-1].location.start, _.location.end, strand=-1)
            if _.type == 'three_prime_UTR':
                utr_features.append(_)
                continue
            new_sub_features.append(_)
        # handle 3' UTR and corresponding exon
        for utr_3 in utr_features:
            if utr_3.location.end > cds_features[-1].location.start:
                if utr_3.location.start > cds_features[-1].location.start:
                    continue
            else:
                for _ in new_sub_features:
                    if _.type == "exon" and _.location == utr_3.location:
                        _.location = FeatureLocation(utr_3.location.start,
                                                     cds_features[-1].location.start - 1,
                                                     strand=-1)
                utr_3.location = FeatureLocation(utr_3.location.start,
                                                 cds_features[-1].location.start - 1,
                                                 strand=-1)
                new_sub_features.append(utr_3)
    _mRNA.sub_features = new_sub_features
    _mRNA.location = FeatureLocation(min(_.location.start for _ in _mRNA.sub_features),
                                     max(_.location.end for _ in _mRNA.sub_features),
                                     strand=_mRNA.strand)
    return _mRNA


def correct_start_codon(mRNA, scaffold_seq):
    _mRNA = copy.deepcopy(mRNA)
    cds_features = [_ for _ in _mRNA.sub_features if _.type == "CDS"]
    cds_features.sort(key=lambda x: x.location.start)
    if _mRNA.strand == 1:
        extention_len = 0
        while True:
            extention_len += 3
            _start_pos = cds_features[0].location.start
            cds_features[0].location = FeatureLocation(_start_pos-3, cds_features[0].location.end, strand=1)
            _codon = scaffold_seq.seq[_start_pos-3: _start_pos]
            if _codon in start_codon or extention_len > 3000:
                break
            if _codon in stop_codon:
                raise TranslationError('First codon could not be found')
        # modify_gene
        new_sub_features = []
        utr_features = []
        for _ in _mRNA.sub_features:
            if _.type == "exon" and _.location.end == cds_features[0].location.end:
                if cds_features[0].location.start > _.location.start:
                    # no UTR or extend UTR boundary
                    _.location = cds_features[0].location
            if _.type == 'five_prime_UTR':
                utr_features.append(_)
                continue
            new_sub_features.append(_)
        # handle 5' UTR
        for utr_5 in utr_features:
            if utr_5.location.end > cds_features[0].location.start:
                if utr_5.location.start > cds_features[0].location.start:
                    continue
            else:
                # handle corresponding exon
                for _ in new_sub_features:
                    if _.location == utr_5.location:
                        _.location = FeatureLocation(utr_5.location.start, cds_features[0].location.start - 1, strand=1)
                utr_5.location = FeatureLocation(utr_5.location.start, cds_features[0].location.start - 1, strand=1)
                new_sub_features.append(utr_5)
    else:
        extention_len = 0
        while True:
            extention_len += 3
            _start_pos = cds_features[-1].location.end
            cds_features[-1].location = FeatureLocation(cds_features[-1].location.start, _start_pos + 3, strand=-1)
            _codon = scaffold_seq.seq[_start_pos: _start_pos + 3].reverse_complement()
            if _codon in start_codon or extention_len > 3000:
                break
            if _codon in stop_codon:
                raise TranslationError('First codon could not be found')
        # modify_gene
        new_sub_features = []
        utr_features = []
        for _ in _mRNA.sub_features:
            if _.type == "exon" and _.location.start == cds_features[-1].location.start:
                if cds_features[-1].location.end > _.location.end:
                    # no UTR or extend UTR boundary
                    _.location = cds_features[-1].location
            if _.type == 'five_prime_UTR':
                utr_features.append(_)
                continue
            new_sub_features.append(_)
        # handle 5' UTR
        for utr_5 in utr_features:
            if utr_5.location.start < cds_features[-1].location.end:
                if utr_5.location.end < cds_features[-1].location.end:
                    continue
            else:
                # for corresponding exon
                for _ in new_sub_features:
                    if _.location == utr_5.location:
                        _.location = FeatureLocation(cds_features[0].location.end + 1, utr_5.location.end, strand=-1)
                utr_5.location = FeatureLocation(cds_features[0].location.end + 1, utr_5.location.end, strand=-1)
                new_sub_features.append(utr_5)
    _mRNA.sub_features = new_sub_features
    _mRNA.location = FeatureLocation(min(_.location.start for _ in _mRNA.sub_features),
                                     max(_.location.end for _ in _mRNA.sub_features),
                                     strand=_mRNA.strand)
    return _mRNA


def correct_phase(mRNA, scaffold_seq):
    _mRNA = copy.deepcopy(mRNA)
    cds_features = [_ for _ in _mRNA.sub_features if _.type == "CDS"]
    cds_features.sort(key=lambda x: x.location.start)
    if _mRNA.strand == 1:
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
    # since the distance from exon to UTRs are always greater than 1 bp, we do not take UTR into consideration here
    if _mRNA.strand == 1:
        for _ in _mRNA.sub_features:
            if _.type == "exon" and _.location.end == cds_features[0].location.end:
                _.location = cds_features[0].location
            new_sub_features.append(_)
        _mRNA.location = FeatureLocation(cds_features[0].location.start, _mRNA.location.end, strand=1)
        _mRNA.sub_features[0].location = FeatureLocation(cds_features[0].location.start,
                                                         _mRNA.location.end,
                                                         strand=1)
    else:
        for _ in _mRNA.sub_features[0].sub_features:
            if _.type == "exon" and _.location.start == cds_features[0].location.start:
                _.location = cds_features[0].location
            new_sub_features.append(_)
        _mRNA.location = FeatureLocation(_mRNA.location.start, cds_features[0].location.end, strand=-1)
        _mRNA.sub_features[0].location = FeatureLocation(_mRNA.location.start,
                                                         cds_features[0].location.end,
                                                         strand=-1)
    _mRNA.sub_features[0].sub_features = new_sub_features
    try:
        cds_seq_full_length.translate(cds=True)
        return _mRNA
    except TranslationError:
        raise TranslationError('These mRNAs need another round correction', _mRNA)


def correct(_gff, _genome):
    """
    :param _gff: gff path
    :param _genome: fasta path
    :return:
    """
    _seqs = SeqIO.to_dict(SeqIO.parse(_genome, 'fasta'))
    _gff = [_ for _ in GFF.parse(_gff, base_dict=_seqs)]
    correct_list = []
    gene_error_dict = defaultdict(list)
    for scaffold in _gff:
        correct_scaffold = SeqRecord(seq="", id=scaffold.id, name=scaffold.name, description=scaffold.description)
        for gene in scaffold.features:
            correct_gene = SeqFeature(location=gene.location,
                                      type='gene',
                                      strand=gene.strand,
                                      id=gene.id,
                                      qualifiers={'ID': [gene.id]})
            correct_gene.sub_features = []
            error_dict = defaultdict(list)
            for mRNA in gene.sub_features:
                try:
                    get_cds(mRNA, scaffold)
                    correct_gene.sub_features.append(mRNA)
                except TranslationError as e:
                    try:
                        if e.args[0].startswith("First codon"):
                            _tmp_mrna = correct_start_codon(mRNA, scaffold)
                            get_cds(_tmp_mrna, scaffold)
                            correct_gene.sub_features.append(_tmp_mrna)
                            error_dict.setdefault('corrected', []).append(mRNA.id)
                        elif e.args[0].startswith('The phase of first CDS is not 0'):
                            # the translation was checked in function
                            _tmp_mrna = correct_phase(mRNA, scaffold)
                            correct_gene.sub_features.append(_tmp_mrna)
                            error_dict.setdefault('corrected', []).append(mRNA.id)
                        elif e.args[0].endswith("is not a stop codon"):
                            _tmp_mrna = correct_stop_codon(mRNA, scaffold)
                            get_cds(_tmp_mrna, scaffold)
                            correct_gene.sub_features.append(_tmp_mrna)
                            error_dict.setdefault('corrected', []).append(mRNA.id)
                        # can not handle for now
                        elif e.args[0] == "Extra in frame stop codon found":
                            error_dict.setdefault('internal', []).append(mRNA.id)
                        elif e.args[0].endswith("is not a multiple of three"):
                            error_dict.setdefault('three', []).append(mRNA.id)
                    except TranslationError as e2:
                        if e2.args[0].startswith('These mRNAs need another round correction'):
                            correct_gene.sub_features.append(e2.args[1])
                            error_dict.setdefault('phase', []).append(mRNA.id)
                        # for second round
                        elif e2.args[0] == "Extra in frame stop codon found":
                            error_dict.setdefault('internal', []).append(mRNA.id)
                        elif e2.args[0].startswith('First codon'):
                            error_dict.setdefault('first2', []).append(mRNA.id)
                        elif e2.args[0].endswith("is not a stop codon"):
                            error_dict.setdefault('final', []).append(mRNA.id)
                        elif e2.args[0].endswith("is not a multiple of three"):
                            error_dict.setdefault('three', []).append(mRNA.id)
                except Exception as e:
                    print(e)
                    print(mRNA.id)
            # handle mRNA and gene relationship
            if not correct_gene.sub_features:
                # Raise error for genes whose all mRNAs have error.
                for _key, value in error_dict.items():
                    _tmp_error = [gene.id + ' ' + _ for _ in value]
                    gene_error_dict[_key] += _tmp_error
            else:
                # check boundary conflict between gene and mRNA.
                gene_start, gene_end = gene.location.start, gene.location.end
                for mRNA in correct_gene.sub_features:
                    if mRNA.location.start < correct_gene.location.start:
                        gene_start = mRNA.location.start
                    if mRNA.location.end > correct_gene.location.end:
                        gene_end = mRNA.location.end
                correct_gene.location = FeatureLocation(gene_start, gene_end, strand=correct_gene.strand)
                correct_scaffold.features.append(correct_gene)
        correct_list.append(correct_scaffold)
    # tidy correct list
    for scaffold in correct_list:
        for gene in scaffold.features:
            for mRNA in gene.sub_features:
                for ele in mRNA.sub_features:
                    if ele.sub_features:
                        ele.sub_features = []
    return correct_list, gene_error_dict


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
    parser.add_argument('-p', '--prefix', type=str, default='genome', help="output file prefix. Default: genome")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    Args = getArgs()
    correct, err = correct(Args.gff, Args.fasta)
    write_result(correct, err, Args.prefix + '_validated')
