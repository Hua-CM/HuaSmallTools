# -*- coding: utf-8 -*-
# @Time : 2020/3/16 19:34
# @Author : Zhongyi Hua
# @FileName: rhizobium_codon.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

from sequence.calculate_GC123 import calculate123
from scipy.stats import pearsonr
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class CodonAnalysis(object):
    def __init__(self, cds_seq, tab, out_dir):
        self.cds_seq = cds_seq
        self.tab = tab
        self.out_dir = out_dir

    @property
    def table(self):
        tmp_tb1 = calculate123(self.cds_seq)
        tmp_tb2 = pd.read_csv(self.tab, sep="\s+")
        tmp_tb2.drop(tmp_tb2[tmp_tb2["Nc"] == "*****"].index, axis=0, inplace=True)
        columns_to_use = list(set(tmp_tb2).difference(set(tmp_tb1)))
        table = pd.merge(tmp_tb1, tmp_tb2[columns_to_use], left_on="seqid", right_on="title", how="right")
        table.drop(["title"], axis=1, inplace=True)
        return table

    @staticmethod
    def enc_curve(xdata):
        ydata = [2 + x + 29 / (x ** 2 + (1 - x) ** 2) for x in xdata]
        return xdata, ydata

    def enc_plot(self):
        plt.figure(figsize=(12, 8), dpi=600)
        plt.xlabel("GC3s")
        plt.ylabel("ENC")
        plt.ylim(0, 70)
        plt.yticks(np.arange(0, 70, 10))
        plt.scatter(self.table["GC3s"], self.table["Nc"], 25, "black")
        x_pred = np.arange(0, 1, 0.01)
        x_pred, y_pred = CodonAnalysis.enc_curve(x_pred)
        plt.plot(x_pred, y_pred, "k-")
        plt.show()
        plt.savefig()

    def neutrality_plot(self):
        plt.figure(figsize=(12, 8), dpi=600)
        plt.xlabel("GC3(%)")
        plt.ylabel("GC12(%)")
        plt.scatter(self.table["GC3"], self.table["GC12"], 25, "black")
        x_data = np.arange(round(min(self.table["GC3"]))-10,
                           round(max(self.table["GC3"]))+10)
        slope, intercept, r_value, p_value, std_err = stats.linregress(self.table["GC3"], self.table["GC12"])
        y_data = [slope * x + intercept for x in x_data]
        plt.plot(x_data, y_data, "r-")
        plt.text(80, 70, "y={:.4f}x+{:.2f}\nR^2={:.4f}".format(slope, intercept, r_value**2))
        plt.show()
        plt.savefig()

    def correlation(self):
        pearson_v = pd.DataFrame(index=self.table.columns.drop(["seqid"]), columns=self.table.columns.drop(["seqid"]))
        pearson_p = pd.DataFrame(index=self.table.columns.drop(["seqid"]), columns=self.table.columns.drop(["seqid"]))
        for var1 in self.table.columns.drop(["seqid"]):
            for var2 in self.table.columns.drop(["seqid"]):
                pearson_v.loc[var1, var2], pearson_p.loc[var1, var2] = pearsonr(self.table[var1], self.table[var2])
        return pearson_v, pearson_p

    def pr2plot(self):
        plt.figure(figsize=(12, 8), dpi=600)
        plt.xlabel("$\frac{A3}{A3+T3}$")
        plt.ylabel("$\frac{G3}{G3+C3}$")
        plt.scatter(self.table[""], self.table["GC12"], 25, "black")




