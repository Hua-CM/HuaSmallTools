# -*- coding: utf-8 -*-
# @Time : 2020/3/16 19:34
# @Author : Zhongyi Hua
# @FileName: rhizobium_codon.py
# @Usage:
# @Note:
# @E-mail: njbxhzy@hotmail.com

from sequence.calculate_GC123 import calculate_codon_analysis
from scipy.stats import pearsonr
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


class CodonAnalysis(object):
    """
    :param cds_seq: cds sequence in fasta format
    :param tab: result file from calculate GC123.py
    :param out_dir: result directory
    """
    def __init__(self, cds_seq, tab, out_dir):
        self.cds_seq = cds_seq
        self.tab = tab
        self.out_dir = out_dir
        self.table = None

    def get_table(self):
        tmp_tb1 = calculate_codon_analysis(self.cds_seq)
        tmp_tb2 = pd.read_csv(self.tab, sep=r"\s+")
        tmp_tb2.drop(tmp_tb2[tmp_tb2["Nc"] == "*****"].index,
                     axis=0, inplace=True)
        columns_to_use = list(set(tmp_tb2).difference(set(tmp_tb1)))
        table = pd.merge(
            tmp_tb1,
            tmp_tb2[columns_to_use],
            left_on="seqid",
            right_on="title",
            how="right")
        table.drop(["title"], axis=1, inplace=True)
        table["Nc"] = table["Nc"].astype(float)
        self.table = table

    @staticmethod
    def enc_curve(xdata):
        ydata = [2 + x + 29 / (x ** 2 + (1 - x) ** 2) for x in xdata]
        return xdata, ydata

    def enc_plot(self):
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_xlabel("GC3s", size="large")
        ax.set_ylabel("ENC", size="large")
        plt.ylim(0, 70)
        plt.yticks(np.arange(0, 70, 10))
        plt.scatter(self.table["GC3s"], self.table["Nc"], 5, "black")
        x_pred = np.arange(0, 1, 0.01)
        x_pred, y_pred = CodonAnalysis.enc_curve(x_pred)
        plt.plot(x_pred, y_pred, "k-")
        plt.savefig(os.path.join(self.out_dir, "enc_plot.pdf"))

    def neutrality_plot(self):
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.set_xlabel("GC3(%)", size="large")
        ax.set_ylabel("GC12(%)", size="large")
        ax.scatter(self.table["GC3"]*100, self.table["GC12"]*100, 5, "black")
        x_data = np.arange(round(min(self.table["GC3"])*100) - 10,
                           round(max(self.table["GC3"])*100) + 10)
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            self.table["GC3"]*100, self.table["GC12"]*100)
        y_data = [slope * x + intercept for x in x_data]
        ax.plot(x_data, y_data, "r-")
        ax.set_xlim([30, 100])
        ax.set_ylim([30, 80])
        ax.set_xticks(np.arange(30, 110, 10))
        ax.set_yticks(np.arange(30, 90, 10))
        plt.text(
            80, 70, "y={:.4f}x+{:.2f}\n$R^2$={:.4f}".format(slope, intercept, r_value**2))
        plt.savefig(os.path.join(self.out_dir, "neutrality_plot.pdf"))

    def correlation(self):
        pearson_v = pd.DataFrame(
            index=self.table.columns.drop(
                ["seqid"]), columns=self.table.columns.drop(
                ["seqid"]))
        pearson_p = pd.DataFrame(
            index=self.table.columns.drop(
                ["seqid"]), columns=self.table.columns.drop(
                ["seqid"]))
        for var1 in self.table.columns.drop(["seqid"]):
            for var2 in self.table.columns.drop(["seqid"]):
                pearson_v.loc[var1, var2], pearson_p.loc[var1, var2] = pearsonr(
                    self.table[var1], self.table[var2])
        pearson_p.to_csv(os.path.join(self.out_dir, "p_value.tsv"), sep="\t", index=False)
        pearson_v.to_csv(os.path.join(self.out_dir, "cor_value.tsv"), sep="\t", index=False)

    def pr2plot(self):
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.set_xlabel("$\\frac{A3}{A3+T3}$", size="large")
        ax.set_ylabel("$\\frac{G3}{G3+C3}$", size="large")
        ax.scatter(self.table["A3"] / (self.table["A3"] + self.table["T3"]),
                   self.table["G3"] / (self.table["G3"] + self.table["C3"]),
                   5,
                   "black")
        ax.axhline(0.5, 0, 1)
        ax.axvline(0.5, 0, 1)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_xlim([0, 1])
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.tick_params(axis="both", labelsize=10)
        ax.set_aspect(1)
        plt.savefig(os.path.join(self.out_dir, "pr2plot.pdf"))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This is the script for codon analysis")
    parser.add_argument('-c', '--cds', required=True,
                        help='<file path> The cds sequence in fasta format')
    parser.add_argument('-t', '--table', required=True,
                        help='<file path> The indices table calculated by codonW')
    parser.add_argument('-o', '--out_dir', required=True,
                        help='<dir_path>  The result directory')
    args = parser.parse_args()
    analysis_instance = CodonAnalysis(args.cds, args.table, args.out_dir)
    analysis_instance.get_table()
    analysis_instance.correlation()
    analysis_instance.enc_plot()
    analysis_instance.neutrality_plot()
    analysis_instance.pr2plot()
    analysis_instance.table.to_csv(os.path.join(args.out_dir, "indices_table.tsv"), sep="\t", index=False)
