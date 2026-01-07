#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: reindex.py
# Created time: 2025/06/26 21:48:39
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python reindex.py -h


import argparse
parser=argparse.ArgumentParser()
infile=parser.add_argument("-i","--input",required=True,type=str,help="输入文件")
outfile=parser.add_argument("-o","--out",required=True,type=str,help="输出文件")
args=parser.parse_args()

infile=args.input
outfile=args.out

data=[]
for i in open(infile):
    line=i.strip().split("\t")
    data.append([int(j) if j!="-" else -9 for j in line])

data=list(zip(*data))

newdata=[]
for i in data:
    a=[]
    line=[]
    for j in i:
        if j==-9:
            continue
        else:
            line.append(j)
    line.sort()
    start2times={}
    start2index={}
    for index,start in enumerate(line):
        if start not in start2times:
            start2index[start]=index
        start2times[start]=start2times.get(start,0)+1
        
    for j in i:
        if j==-9:
            a.append("-")
        else:
            if start2times[j]==1:
                a.append(start2index[j])
            else:
                a.append(start2index[j])
                start2index[j]=start2index[j]+1
                start2times[j]=start2times[j]-1
    newdata.append(a)

newdata=list(zip(*newdata))

with open(outfile,"w") as f:
    for i in newdata:
        f.write("\t".join([str(j) for j in i])+"\n")