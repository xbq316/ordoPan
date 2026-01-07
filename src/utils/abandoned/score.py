#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: score.py
# Created time: 2025/07/03 17:24:11
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python score.py -h


import networkx as nx
import argparse
from collections import defaultdict

parser=argparse.ArgumentParser()
parser.add_argument("input",type=str,help="输入 num 文件")
args=parser.parse_args()

numFile=args.input

a=[]
for i in open(numFile):
    i=i.strip()
    if not i :
        continue
    line=i.split("\t")
    newline=[int(x) if x!="-" else -9 for x in line]
    a.append(newline)

a=list(zip(*a))
U=set()
for i in a:
    i=[-float("inf"),*i,float("inf")]
    b=[]
    for gene_index,j in enumerate(i):
        if j==-9:
            continue
        b.append((gene_index,j))

    b.sort(key=lambda x: x[1])
    len_b = len(b)
    for i in range(len_b-1):
        this_index=b[i][0]
        next_index=b[i+1][0]
        U.add((this_index,next_index))


U_=set()
for i in a:
    i=[-float("inf"),*i,float("inf")]
    b=[]
    for gene_index,j in enumerate(i):
        if j==-9:
            continue
        b.append((gene_index,j))

    len_b = len(b)
    for i in range(len_b-1):
        this_index=b[i][0]
        next_index=b[i+1][0]
        U_.add((this_index,next_index))

score=(1-len(U&U_)/len(U|U_))*100
print(f"error rate: {score:.2f}%")

adjacency_accuracy=len(set([i for i in U_ if i in U]))/len(U)
print(f"adjacency accuracy: {adjacency_accuracy:.2f}")
