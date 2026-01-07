#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: geneindex2genelist.py
# Created time: 2025/05/27 10:48:37
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python geneindex2genelist.py -h



import sys
import argparse

helptext=r"""
输入文件格式：
-  - spe3-gene1
spe1-gene1  spe2-gene1 spe3-gene2
spe1-gene2  spe2-gene2 spe3-gene3

输出文件格式：
spe1-gene1
spe2-gene1
spe3-gene2

"""

parser=argparse.ArgumentParser()
parser.add_argument("-i","--input",required=True,type=str,help="输入文件")
parser.add_argument("-o","--out",required=True,type=str,help="输出文件")
args=parser.parse_args()

infile=args.input
outfile=args.out


with open(outfile,"w") as f:
    geneindex=open(infile).readlines()
    for i in geneindex:
        line=i.strip().split("\t")
        for genes in line:
            if "," in genes:
                geneList=genes.split(",")
                gene=geneList[0]
                break
            elif genes=="-":continue
            else:
                gene=genes
                break
        f.write(f"{gene}\n")
