#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: bed2wgdi.py
# Created time: 2025/07/14 17:24:28
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python bed2wgdi.py -h




import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Convert bed(format 6) file to wgdi(gff)')
parser.add_argument("-i",'--bed',required=True, help='bed file')
parser.add_argument("-s",'--species',required=True, help='gff文件新基因名称前缀')
parser.add_argument("-o",'--out',required=True, help='prefix of output file')
args = parser.parse_args()

infile=args.bed
out_prefix=args.out
species=args.species


chro2gff = defaultdict(list)
for i in open(infile):
    i=i.strip()
    if not i: continue
    line=i.split('\t')
    chro=line[0]
    start=int(line[1])
    end=int(line[2])
    name=line[3]
    chro2gff[chro].append([start,end,name])

for i in chro2gff:
    chro2gff[i].sort(key=lambda x:x[0])

with open(out_prefix+".lens","w") as f:
    for chro in chro2gff:
        chromosomeLength=chro2gff[chro][-1][1]
        geneNum=len(chro2gff[chro])
        f.write(f"{chro}\t{chromosomeLength}\t{geneNum}\n")

with open(out_prefix+".gff","w") as f:
    for chro in chro2gff:
        for index,(start,end,name) in enumerate(chro2gff[chro]):
            chrome=f"{chro}"
            newGeneName=f"{species}{chro}g{index+1:07}"
            start=f"{start}"
            end=f"{end}"
            strand=f"+"
            f.write(f"{chrome}\t{newGeneName}\t{start}\t{end}\t{strand}\t{index+1}\t{name}\n")
