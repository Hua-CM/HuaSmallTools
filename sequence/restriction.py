
from Bio.Restriction import Restriction as Rst
from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict
from functools import reduce
from itertools import chain
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import argparse
import os


def get_seq(seq, _pos, _flank):
    """
    :param seq: Bio.seq
    :param _pos: the natural position.
    get_seq would handle the conflict between natural position and computational postion
    :param _flank:
    :return:
    """
    if _pos <= _flank:
        start = 0
        end = _pos + _flank
    elif _pos >= len(seq.seq) - _flank:
        start = _pos - _flank - 1
        end = len(seq.seq)
    else:
        start = _pos - _flank - 1
        end = _pos + _flank
    return seq.seq[start:end]


def create_enzyme(name):
    e_types = [x for t, (x, y) in typedict.items() if name in y][0]
    enzyme_types = tuple(getattr(Rst, x) for x in e_types)

    return Rst.RestrictionType(name, enzyme_types, rest_dict[name])


def get_site_interval(name):
    enzyme = create_enzyme(name)
    rco_pos, *_, rco_site = enzyme.characteristic()
    return {name: np.array([rco_pos, len(rco_site)-rco_pos])}


def analysis_pos(_seq, _snp_pos, _sub_nuc, _rb, _rbi):
    """
    NOTE:
    1. Bio.Restriction could not handle degenerate bases
    2. Bio.Restriction did not scan complement sequence
    :param _seq: Bio.Seq.Seq
    :param _snp_pos: int: SNP position in sequence (NATURE Position !!)
    :param _sub_nuc,
    :param _rb: Bio.Restriction.RestrictionBatch
    :param _rbi:
    :return:  list: enzymes' names
    """
    def parse_search(rb_search):
        return defaultdict(list, {k.__name__: v for k, v in rb_search.items()})

    _sub_seq = list(_seq)
    _sub_seq[_snp_pos] = _sub_nuc
    _sub_seq = Seq(''.join(_sub_seq))
    res_lst = []
    for __seq in [_seq, _seq.complement(), _sub_seq, _sub_seq.complement()]:
        res_lst.append(parse_search(_rb.search(__seq)))
    rb_total = {enzyme: list(set(list(chain(*[_[enzyme] for _ in res_lst])))) for enzyme in res_lst[0].keys()}
    res = {}
    res_enzymes = []
    for _enzyme, _pos_list in rb_total.items():
        _enzyme_list = []
        if len(_pos_list) == 1:
            _enzyme_list.append(_pos_list[0] + _rbi[_enzyme])
        res[_enzyme] = _enzyme_list
    for _enzyme, _pos in res.items():
        if len(_pos) > 0:
            if _pos[0][0] < _snp_pos < _pos[0][1]:
                res_enzymes.append(_enzyme)
    return res_enzymes


def main(args):
    enzyme_list = ['ApaI', 'BamHI', 'BglII', 'DraI', 'EcoRI', 'EcoRV', 'HindIII',
                   'KpnI', 'NcoI', 'NdeI', 'NheI', 'NotI', 'PstI', 'SacI', 'SalI',
                   'SphI', 'XbaI', 'XhoI']
    rb = reduce(lambda x, y: x + y, map(create_enzyme, enzyme_list))
    rbi = defaultdict()
    for _ in map(get_site_interval, enzyme_list):
        rbi.update(_)
    results = []
    seq_results = []
    # parse genome sequence
    geo_seq = SeqIO.to_dict(SeqIO.parse(args.input, 'fasta'))
    # prepare per sequence
    with open(args.bed) as f_in:
        for snp in f_in.read().split('\n'):
            try:
                _id, _pos, _sub = snp.split()
                _seq = get_seq(geo_seq[_id], int(_pos), args.flank)
                enzymes = analysis_pos(_seq, args.flank, _sub, rb, rbi)
                if enzymes:
                    results.append(f'{_id}\t{_pos}\t{",".join(enzymes)}')
                    seq_results.append(SeqRecord(id='_'.join([_id, str(int(_pos) - args.flank - 1), str(int(_pos) + args.flank)]),
                                                 seq=_seq))
                else:
                    continue
            except:
                continue
    with open(args.output, 'w') as f_out:
        f_out.write('\n'.join(results))
    SeqIO.write(seq_results, os.path.splitext(args.output)[0] + '.fasta', 'fasta')


def getargs():
    parser = argparse.ArgumentParser(description="Detect whether SNP site is on restriction sites")
    parser.add_argument('-i', '--input', required=True,
                        help='<file_path> Genome sequence file in fasta format')
    parser.add_argument('-b', '--bed', required=True,
                        help='<file_path> a bed-like file: chromosome_id/position/substitution nucleotide')
    parser.add_argument('-l', '--length', dest='flank', default=400, type=int,
                        help='<int> output flanking sequence length. Default: 400')
    parser.add_argument('-o', '--output', required=True,
                        help='<file_path>  Output file path')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main(getargs())
