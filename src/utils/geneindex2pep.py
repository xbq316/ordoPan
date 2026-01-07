#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: geneindex2pep.py
# Created time: 2025/12/26 17:25:55
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python geneindex2pep.py -h



import os
from Bio import SeqIO 
import argparse



NONE_GENE_SET=set([".","-",""])
SEP="\t"

helptext=r"""
This script extracts peptide sequences for genes from multiple samples based on a gene index file.
"""

def core(line: list[str],core_rate: float=1) -> bool:
    if core_rate<0 or core_rate>1:
        raise ValueError("core_rate must be between 0 and 1")
    if core_rate==1:
        for j in line:
            if j in NONE_GENE_SET:
                return False
        return True
    else:
        present_num=0
        total_num=len(line)
        for j in line:
            if j not in NONE_GENE_SET:
                present_num+=1
        if present_num/total_num>=core_rate:
            return True
        else:
            return False

parser = argparse.ArgumentParser(description=helptext, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-p", "--pep", type=str, required=True, help="Directory containing BED files for each sample")
parser.add_argument("-s", "--sample", type=str, required=True, help="File containing the list of sample names(The order of the samples must match the order of samples in the geneindex file.)")
parser.add_argument("--pep-suffix", type=str, default=".pep", help="Suffix of the peptide files (default: .pep)")
parser.add_argument("-g", "--geneindex", type=str, required=True, help="Gene index file")
parser.add_argument("-o", "--output", type=str, required=True, help="Output file")
parser.add_argument("--core", action="store_true", help="Only extract peptides for core genes,default: extract peptides for all genes")
parser.add_argument("-r","--rate",type=float,default=1,help="Core gene definition: the minimum presence rate across samples (default: 1, meaning only genes present in all samples are considered core genes)")
args = parser.parse_args()

pep_suffix=args.pep_suffix

sample_list=[]
for i in open(args.sample):
    sample_list.append(i.strip())
# 检查文件是否存在
for sample in sample_list:
    pep_file=os.path.join(args.pep, f"{sample}{pep_suffix}")
    if not os.path.exists(pep_file):
        raise FileNotFoundError(f"PEP file for sample {sample} not found: {pep_file}")
        print(f"Maybe use the --pep-suffix option")
        exit(1)
# (gene,sample_index)
gene_list=[] 
sampleindex2num={} # 记录 gene_list 中每个样本需要提取的基因数量
for i in open(args.geneindex):
    i=i.strip()
    if not i:
        continue
    line=i.split(SEP)
    
    if args.core:
        # 只提取核心基因
        this_line_is_core=core(line,core_rate=args.rate)
        if not this_line_is_core:
            continue
        else:
            pass
    else:
        # 提取所有基因
        pass
    for sample_index,j in enumerate(line):
        if j in NONE_GENE_SET:
            continue
        elif "," in j:
            gene_sublist=j.split(",")
            gene=gene_sublist[0]
            gene_list.append((gene,sample_index))
            sampleindex2num[sample_index]=sampleindex2num.get(sample_index,0)+1
            break
        else:
            gene_list.append((j,sample_index))
            sampleindex2num[sample_index]=sampleindex2num.get(sample_index,0)+1
            break

gene_sample2pep={}
for gene,sample_index in gene_list:
    gene_sample2pep[(gene,sample_index)]=""


for sample_index in sampleindex2num.keys():
    sample_name=sample_list[sample_index]
    pep_file=os.path.join(args.pep, f"{sample_name}{pep_suffix}")
    this_sample_gene_num=sampleindex2num[sample_index]
    a=0
    for record in SeqIO.parse(pep_file, "fasta"):
        gene=str(record.id)
        if (gene,sample_index) in gene_sample2pep:
            gene_sample2pep[(gene,sample_index)]=str(record.seq)
            a+=1
            if a>=this_sample_gene_num:
                break

with open(args.output, "w") as out_f:
    for (gene,sample_index) in gene_list:
        pep=gene_sample2pep[(gene,sample_index)]
        if not pep:
            Warning(f"Warning: No peptide found for gene {gene} in sample({pep_suffix}) {sample_list[sample_index]}")
        out_f.write(f">{gene}\n{pep}\n")
