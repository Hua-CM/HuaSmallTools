# -*- coding: utf-8 -*-
# @Time : 2020/11/23 10:10
# @Author : Zhongyi Hua
# @FileName: judge_path.py
# @Usage: Judge the "correct" GetOrganelle result path
# @Note:
# @E-mail: njbxhzy@hotmail.com
import shutil
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from tempfile import gettempdir
from tqdm import tqdm

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

logger = logging.getLogger('judge_path')
logger.setLevel(level=logging.INFO)
logger.addHandler(logging.StreamHandler())


def parse_ssc(gfa: Path):
    """Get SSC sequence from the gfa Path

    Args:
        gfa (Path): The *.gfa file recording segment sequences. 

    Returns:
        ssc_seq (SeqRecord): The SSC SeqRecord 
    """
    seq_dict = {}
    count_dict = defaultdict(int)
    ssc_name = ''
    ssc_length = 9999999999 # Just a big number
    for line in gfa.read_text().split('\n'):
        _tmp = line.split('\t')
        if line.startswith('L'):
            count_dict[_tmp[1]] += 1 
            count_dict[_tmp[3]] += 1
        elif line.startswith('S'):
            seq_dict[_tmp[1]] = _tmp[2]
    for seq_name, count in count_dict.items():
        if count == 2 and len(seq_dict.get(seq_name)) < ssc_length:
            ssc_name = seq_name
            ssc_length = len(seq_dict.get(seq_name))
    ssc_seq = SeqRecord(id=ssc_name, 
                        seq=Seq(seq_dict.get(ssc_name)),
                        description='')
    return ssc_seq



class ChloroGraph:
    def __init__(self, res_dir:Path, ref: Path, tmp: Path = Path(gettempdir())) -> None:
        """_summary_

        Args:
            res_dir (Path): the GetOrganelle result directory
            ref (Path): The reference fasta path
            tmp (Path): The temporary directory path (for blast)
        """
        self.res_dir  = res_dir
        self.ref = ref
        self.tmp = tmp / res_dir.name
    
    def _initiate(self) -> None:
        """
        Make temporary directory
        """
        if not self.tmp.exists():
            self.tmp.mkdir()

    def _autoclean(self) -> None:
        """
        Remove temporary directory
        """
        shutil.rmtree(self.tmp)
    
    def _detect_file(self):
        """Detect the specific files' path in result directory

        Args:
            res_dir (Path): The assemble result directory path.

        Returns:
            path1 (Path): *.complete.graph1.1.path_sequence.fasta
            path2 (Path): *.complete.graph1.2.path_sequence.fasta
            gfa (Path): *.complete.graph1.selected_graph.gfa
        """
        for file in self.res_dir.iterdir():
            file_str = str(file)
            if file_str.endswith('graph1.1.path_sequence.fasta'):
                path1 = Path.resolve(file)
            elif file_str.endswith('graph1.2.path_sequence.fasta'):
                path2 = Path.resolve(file)
            elif file_str.endswith('graph1.selected_graph.gfa'):
                gfa = Path.resolve(file)
        res_vars = locals()
        if 'path1' in res_vars and 'path2' in res_vars and 'gfa' in res_vars:
            return path1, path2, gfa
        else: 
            logger.info('This result directory is not an complete graph!')
            return None
    
    def _blast_module(self, ssc_seq_path, ref):
        """Identify SSC direction

        Args:
            ssc_seq_path (_type_): _description_
            ref (_type_): _description_
        Returns:
            strand (str): '+' or '-'
        """
        mkdb_cmd = NcbimakeblastdbCommandline(
            input_file=ref,
            dbtype='nucl',
        )
        mkdb_cmd()
        run_cmd = NcbiblastnCommandline(
            query = ssc_seq_path,
            db = ref,
            max_target_seqs = 10,
            max_hsps=1,
            outfmt='6 std sstrand',
            out=self.tmp / 'ssc_blast.out'
        )
        run_cmd()
        # Parse the output
        strand_dict = defaultdict(int)
        for line in (self.tmp / 'ssc_blast.out').read_text().strip().split('\n'):
            strand_dict[line.split('\t')[-1]] += 1
        if strand_dict.get('plus', 0) > strand_dict.get('minus', 0):
            return '+'
        else:
            return '-'
    
    @staticmethod
    def _judge_file(strand:str, ssc_name:str, file1: Path, file2: Path):
        """_summary_

        Args:
            strand (str): The correct SSC strand.
            ssc_name (str): SSC fragment name.
            file1 (Path): The graph1 fasta file path.
            file2 (Path): The graph2 fasta file path.
        Returns:
            selected_graph: (SeqRecord): The selected graph
        """

        def parse_header(file):
            with open(file) as f_in:
                header = f_in.readline()
                fragments = header.split(',')
                for fragment in fragments:
                    if fragment.startswith(ssc_name) and fragment.endswith(strand):
                        return True
            return False
        
        if parse_header(file1):
            return file1
        if parse_header(file2):
            return file2
     
    def select_graph(self):
        try:
            self._initiate()
            _tmp = self._detect_file()
            if not _tmp:
                logger.info('Skip!')
                return None
            path1, path2, gfa = _tmp
            ssc_seq = parse_ssc(gfa)
            # Write SSC seq to temporary directory
            SeqIO.write(ssc_seq, self.tmp / 'ssc.fasta', 'fasta')
            # soft link file
            if (self.tmp / 'ref.fasta').exists():
                Path.unlink((self.tmp / 'ref.fasta')) 
            (self.tmp / 'ref.fasta').symlink_to(self.ref, )
            # Do BLAST
            strand = self._blast_module(self.tmp / 'ssc.fasta', self.tmp / 'ref.fasta')
            # Judge File
            selected_graph_path = self._judge_file(strand, ssc_seq.id, path1, path2)
            # Rename and output
            selected_graph = SeqIO.read(selected_graph_path, 'fasta')
            selected_graph.id = self.res_dir.name
            selected_graph.name = ''
            selected_graph.description = ''
        except:
            selected_graph = None
        finally:
            self._autoclean()
        return selected_graph

def parse_args():
    parser = argparse.ArgumentParser(
         description='This script was for selecting "correct" graph')
    parser.add_argument('-i', '--input', type=Path, required=True,
                                 help='<File path>  The meta information. Two coloumns: result directory / reference fasta file')
    parser.add_argument('-o', '--output', type=Path, required=True,
                                 help='<File path>  The output result directory')
    args = parser.parse_args()
    return args

def main(args):
    if not args.output.exists():
        args.output.mkdir()
    lines = args.input.read_text().strip().split('\n')
    for line in tqdm(lines):
        try:
            res_dir, ref_fasta = line.split('\t')
            res_dir = Path(res_dir)
            ref_fasta = Path(ref_fasta)
            ins = ChloroGraph(res_dir, ref_fasta)
            selected_graph = ins.select_graph()
            SeqIO.write(selected_graph,
                        args.output / (res_dir.name + '.fasta') ,
                        'fasta')
            logger.info(f'{res_dir.name} done')
        except Exception as e:
            logger.info(f'{res_dir.name} is error! check it！')
        

if __name__ == '__main__':
    main(parse_args())
