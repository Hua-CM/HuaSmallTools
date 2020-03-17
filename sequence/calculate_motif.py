# -*- coding: utf-8 -*-
# @Time : 2020/3/17 11:36
# @Author : Zhongyi Hua
# @FileName: calculate_AT123.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

from operator import itemgetter
import numpy as np


def motif123(seq, motif):
    """
    Calculate A/T/C/G content: total, for first, second and third positions.
    :param seq: a seq
    :param motif: a list of motif. eg.["A"], ["T"], ["AT"], ["GC"]
    :return:a tuple of four floats (percentages between 0 and 100) for the
    entire sequence, and the three codon positions.
    """
    d = {}
    for nt in ["A", "T", "G", "C"]:
        d[nt] = [0, 0, 0]

    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if len(codon) < 3:
            codon += "  "
        for pos in range(0, 3):
            for nt in ["A", "T", "G", "C"]:
                if codon[pos] == nt or codon[pos] == nt.lower():
                    d[nt][pos] += 1
    gc = {}
    gcall = 0
    nall = 0
    if motif.__len__() == 1:
        for i in range(0, 3):
            try:
                n = d["G"][i] + d["C"][i] + d["T"][i] + d["A"][i]
                gc[i] = np.array(itemgetter(*motif)(d))[i].sum() * 100.0 / n
            except Exception:
                gc[i] = 0

            gcall = gcall + np.array(itemgetter(*motif)(d))[i].sum()
            nall = nall + n
    else:
        for i in range(0, 3):
            try:
                n = d["G"][i] + d["C"][i] + d["T"][i] + d["A"][i]
                gc[i] = np.array(itemgetter(*motif)(d))[:, i].sum() * 100.0 / n
            except Exception:
                gc[i] = 0

            gcall = gcall + np.array(itemgetter(*motif)(d))[:, i].sum()
            nall = nall + n
    gcall = 100.0 * gcall / nall
    return gcall, gc[0], gc[1], gc[2]
