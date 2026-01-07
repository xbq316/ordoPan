#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: split.py
# Created time: 2025/06/27 15:36:17
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python split.py -h


from collections import defaultdict
import os
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-i","--input","--geneindex",required=True,help="geneindex 文件(行数与 vote 一致即可)")
parser.add_argument("-o","--out",required=True,help="输出目录")
parser.add_argument("-v","--vote",required=True,help="vote文件")
parser.add_argument("-s","--suffix",default="txt",help="文件后缀")
parser.add_argument("-c","--chro",required=True,help="常染色体文件")
args=parser.parse_args()


votefile=args.vote
infile=args.input
outdir=args.out
suffix=args.suffix

if not os.path.exists(outdir):
    os.mkdir(outdir)


chros=set()
for chro in open(args.chro):
    chro=chro.strip()
    chros.add(chro)

voteList=[]
for i in open(votefile):
    voteList.append(i.strip())

chro2genes = defaultdict(list)

for index,i in enumerate(open(infile)):
    chro=voteList[index]
    if chro not in chros:
        chro="Scaffold"
    chro2genes[chro].append(i)

for chro,genes in chro2genes.items():
    with open(os.path.join(outdir,f"{chro}{suffix}"),"w") as f:
        f.write("".join(genes))
