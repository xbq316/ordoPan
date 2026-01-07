#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: geneindex2pair.py
# Created time: 2025/11/24 20:39:00
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python geneindex2pair.py -h


import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-i","--input",required=True,type=str,help="输入文件")
parser.add_argument("-o","--out",required=True,type=str,help="输出文件")
args=parser.parse_args()

infile=args.input
outfile=args.out


geneindex=[]
for i in open(infile):
    line=i.strip()
    if not line:
        continue
    line=line.split("\t")
    geneindex.append(line)

sample_num=len(geneindex[0])

with open(outfile,"w") as out:
    for i in geneindex:
        
        a=[]
        for j in i:
            if j=="-":
                continue
            else:
                a.append(j)
        a_size=len(a)
        for m in range(a_size-1):
            genes1=a[m].split(",")
            genes2=a[m+1].split(",")
            for g1 in genes1:
                for g2 in genes2:
                    out.write(f"{g1}\t{g2}\n")
            