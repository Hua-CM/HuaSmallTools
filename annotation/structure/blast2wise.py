# -*- coding: utf-8 -*-
# @Time : 2021/6/26 17:48
# @Author : Zhongyi Hua
# @FileName: primer_utils.py
# @Usage:
# @Note:

# @E-mail: njbxhzy@hotmail.com

# version 0.1: it will run genewise parallel with Python multiprocessing library 
# version 0.2: run before testing the existed file

import os
import sys
import subprocess as sp
import pandas as pd
from multiprocessing import Pool
import os
from tempfile import mktemp


def parse_genblasta(_file):
    query = 'None'
    results = []
    with open(_file) as f_in:
        lines = f_in.read().strip().split('\n')
        module = []
        for line in lines:
            if line.startswith(">"):
                query = line.lstrip(">")
                results += module
                module = []
            else:
                eles = line.split('\t')
                module.append({'query': query,
                               'chr': eles[6],
                               'start': eles[7],
                               'end': eles[8],
                               'strand': eles[10]})
    results += module
    results = pd.DataFrame(results)
    results = results[~(results['query'] == 'None')]
    return results


def genewise(gblasta, _tmp, genome_path, protein_path):
    gffpath = mktemp(dir=_tmp) + '.gff'
    logpath = gffpath +'.log'
    old_prt_path = gffpath + '.target'
    os.system('touch {}'.format(old_prt_path))
    for _idx, _row in gblasta.iterrows():
        prt_path = os.path.join(_tmp, _row['query'] + '.fa')
        if not os.path.exists(prt_path):
            # substract protein sequence
            prtshell = "seqkit faidx {} '{}' -o {}".format(protein_path, _row['query'], prt_path)
            sp.run(prtshell, shell=True)
            os.remove(old_prt_path)
        d_start = _row['start'] if int(_row['start']) - 2000 < 0 else int(_row['start']) - 2000
        d_end = _row['end'] if int(_row['start']) - 2000 < 0 else int(_row['end']) + 2000
        region = str(d_start) + "-" + str(d_end)
        dnapath = os.path.join(_tmp, _row['chr'] + '_' + region + '.fa')
        dnashell = "seqkit faidx {} '{}:{}' -o {}".format(genome_path, _row['chr'], region, dnapath)
        sp.run(dnashell, shell=True)
        genewiseshell = "genewise {} {} -quiet -gff {} >> {} 2>> {}".\
            format(prt_path, dnapath, "-tfor" if _row['strand'] == "1" else "-trev", gffpath, logpath)
        sp.run(genewiseshell, shell=True)
        old_prt_path = prt_path
        os.remove(dnapath)
    os.remove(old_prt_path)


def main():
    genome = sys.argv[1]
    protein = sys.argv[2]
    genblasta = sys.argv[3]
    tmp = sys.argv[4]
    genblasta_res = parse_genblasta(genblasta)
    genewise(genblasta_res, tmp, genome, protein)


if __name__ == "__main__":
    main()
