# -*- coding: utf-8 -*-
# @Time : 2020/6/4 15:52
# @Author : Zhongyi Hua
# @FileName: calculate_GC_window.py
# @Usage: Use slide windows to calculate GC and GC skew for Circos Plot
# @Note:
# @E-mail: njbxhzy@hotmail.com

from os import path
from Bio.SeqUtils import GC, GC_skew
from Bio import SeqIO
from pandas import DataFrame


def gc_window(genome, window=1000):
    interval = [x for x in range(int(len(genome.seq) / window) + 1)]
    result_dict = dict()
    for i in interval[:-1]:
        gc = GC(genome.seq[interval[i]*window:interval[i+1]*window])
        gc_cell = {i: {"chr": genome.id,
                       "start": interval[i]*window+1,
                       "end": interval[i+1]*window,
                       "GC": gc
                       }
                   }
        result_dict.update(gc_cell)
    gc = GC(genome.seq[interval[-1] * window: len(genome.seq)])
    gc_cell = {interval[-1]: {"chr": genome.id,
                              "start": interval[-1] * window + 1,
                              "end": len(genome.seq),
                              "GC": gc
                              }
               }
    result_dict.update(gc_cell)
    return DataFrame.from_dict(result_dict, "index")


def gc_skew_window(genome, window=1000):
    interval = [x for x in range(int(len(genome.seq) / window) + 1)]
    result_dict = dict()
    for i in interval[:-1]:
        gc_skew = GC_skew(genome.seq[interval[i] * window: interval[i + 1] * window])[0]
        gc_skew_cell = {i: {"chr": genome.id,
                            "start": interval[i] * window + 1,
                            "end": interval[i + 1] * window,
                            "GC_skew": gc_skew
                            }
                        }
        result_dict.update(gc_skew_cell)
    gc_skew = GC_skew(genome.seq[interval[-1] * window: len(genome.seq)])[0]
    gc_skew_cell = {interval[-1]: {"chr": genome.id,
                                   "start": interval[-1] * 100 + 1,
                                   "end": len(genome.seq),
                                   "GC_skew": gc_skew
                                   }
                    }
    result_dict.update(gc_skew_cell)
    return DataFrame.from_dict(result_dict, "index")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This is the script for calculating GC and GC_skew window which could "
                                                 "be used in Circos")
    parser.add_argument('-i', '--input_file', required=True,
                        help='<filepath>  The genome fasta file')
    parser.add_argument('-w', '--window', default=1000, type=int,
                        help='<integer> Window length')
    parser.add_argument('-o', '--output_dir', required=True)
    args = parser.parse_args()
    genomeG = SeqIO.read(args.input_file, "fasta")
    df1 = gc_window(genomeG, args.window)
    df2 = gc_skew_window(genomeG, args.window)
    df1.to_csv(path.join(args.output_dir, "GC.tsv"), sep="\t", index=False)
    df2.to_csv(path.join(args.output_dir, "GC_skew.tsv"), sep="\t", index=False)
