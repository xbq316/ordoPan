#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: geneindex2old.py
# Created time: 2025/11/24 21:17:23
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python geneindex2old.py -h



import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-i","--input",required=True,type=str,help="输入文件")
parser.add_argument("-s","--sample",required=True,type=str,help="样本列表文件")
parser.add_argument("-o","--out",required=True,type=str,help="输出文件")
args=parser.parse_args()

infile=args.input
outfile=args.out
samplefile=args.sample

sample_list=[]
for i in open(samplefile):
    line=i.strip()
    if not line:
        continue
    sample_list.append(line)
sample_name_len=[len(i) for i in sample_list]

geneindex=[]
for i in open(infile):
    line=i.strip()
    if not line:
        continue
    line=line.split("\t")
    geneindex.append(line)

if len(sample_list)!=len(geneindex[0]):
    raise ValueError("样本数量与基因组数量不匹配，请检查输入文件和样本列表文件")

with open(outfile,"w") as out:
    for line in geneindex:
        this_new_line=[]
        for sample_index,sample in enumerate(sample_list):
            genes=line[sample_index]
            if genes=="-":
                this_new_line.append("-")
            elif "," in genes:
                a=[]
                genes_list=genes.split(",")
                for gene in genes_list:
                    new_gene=gene[sample_name_len[sample_index]+1:]
                    a.append(new_gene)
                this_new_line.append(",".join(a))
            else:
                new_gene=genes[sample_name_len[sample_index]+1:]
                this_new_line.append(new_gene)
        out.write("\t".join(this_new_line)+"\n")