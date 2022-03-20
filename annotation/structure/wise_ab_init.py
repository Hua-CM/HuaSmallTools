import io
import sys
import re
import copy

filepath = sys.argv[1]
fp = open(filepath, 'r')

exon_count = 0
choose = 0
for line in fp:
    if line.startswith("#"): # or col_id in match_list:
        continue
    elif re.search(r"\tmatch\t", line):
        # construct gene feature
        choose += 1
        gene_List = line.strip().split("\t")
        mRNA_List = copy.deepcopy(gene_List)
        # Public
        gene_start = gene_List[0].split(':')[1].split('-')[0]
        gene_List[0] = gene_List[0].split(':')[0]
        gene_List[3] = str(int(gene_start) + int(gene_List[3]) - 1)
        gene_List[4] = str(int(gene_start) + int(gene_List[4]) - 1)
        if gene_List[6] == '-':
            gene_List[3], gene_List[4] = gene_List[4], gene_List[3]
        mRNA_List[0], mRNA_List[3], mRNA_List[4], mRNA_List[5] = gene_List[0], gene_List[3], gene_List[4], gene_List[5]
        # gene
        gene_List[2] = "gene"
        gene_id = '_'.join([gene_List[0], gene_List[3], gene_List[4]]) + "_gene_" + str(choose)
        gene_List[8] = "ID={}".format(gene_id)
        # mRNA
        mRNA_List[2] = 'mRNA'
        mRNA_id = '_'.join([mRNA_List[0], mRNA_List[3], mRNA_List[4]]) + "_mRNA_" + str(choose)
        mRNA_List[8] = "ID={};Parent={}".format(mRNA_id, gene_id)
        print("\t".join(gene_List))
        print("\t".join(mRNA_List))
        exon_count = 0
    elif re.search(r"\tcds\t", line):
        exon_count += 1
        CDS_List = line.strip().split('\t')
        exon_List = copy.deepcopy(CDS_List)
        # Public
        exon_start = exon_List[0].split(':')[1].split('-')[0]
        exon_List[0] = exon_List[0].split(':')[0]
        exon_List[3] = str(int(exon_start) + int(exon_List[3]) - 1)
        exon_List[4] = str(int(exon_start) + int(exon_List[4]) - 1)
        if exon_List[6] == '-':
            exon_List[3], exon_List[4] = exon_List[4], exon_List[3]     
        CDS_List[0], CDS_List[3], CDS_List[4], CDS_List[5] = exon_List[0], exon_List[3], exon_List[4], exon_List[5]
        # exon
        exon_List[7] = '.'
        exon_List[2] = 'exon'
        exon_id = mRNA_id + '_exon_' + str(exon_count)
        exon_List[8] = 'ID={};Parent={}'.format(exon_id, mRNA_id)
        # CDS 
        CDS_List[2] = 'CDS' #upper
        cds_id = mRNA_id + '_cds_' + str(exon_count)
        CDS_List[8] = "ID={};Parent={}".format(cds_id, mRNA_id)
        print("\t".join(exon_List))
        print("\t".join(CDS_List))
    else:
        continue
fp.close()
