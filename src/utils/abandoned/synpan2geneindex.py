#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: synpan2geneindex.py
# Created time: 2025/06/26 20:30:34
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python synpan2geneindex.py -h

import os
import argparse
import re
import warnings

parser=argparse.ArgumentParser("Convert Synpan results to GeneIndex format")
parser.add_argument("-l","--sample_list",help="Sample list file",required=True)
parser.add_argument("-s","--synpan",help="Synpan output file: *.SG.pan",required=True)
parser.add_argument("-o","--output",help="Output file name",required=True)
args=parser.parse_args()

sample_list_file=args.sample_list
sample_list=[]
for i in open(sample_list_file):
    i=i.strip()
    if not i:
        continue
    sample_list.append(i)

synpan_file=args.synpan

synpan=[]
for i in open(synpan_file):
    line=i.strip("\n").split("\t")[4:]
    for j in range(len(line)):
        if line[j] =="":
            line[j]="-"
        elif "," in line[j]:
            a=[]
            for k in line[j].split(","):
                prifix=f"{sample_list[j]}_"
                indexStart=len(prifix)
                genename=k[indexStart:]
                a.append(genename)
            line[j]=",".join(a)
        else:
            prifix=f"{sample_list[j]}_"
            indexStart=len(prifix)
            genename=line[j][indexStart:]
            line[j]=genename
    synpan.append(line)


with open(args.output,"w") as f:
    for i in synpan:
        f.write("\t".join(i)+"\n")

