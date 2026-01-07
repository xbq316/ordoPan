#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: pep2bed.py
# Created time: 2025/09/23 16:14:24
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python pep2bed.py -h



from Bio import SeqIO
import argparse
import os


parser=argparse.ArgumentParser()
parser.add_argument("-i","--pep",required=True,help="pep file")
parser.add_argument("-o","--out",required=True,help="out file:bed file")
args=parser.parse_args()

pepfile=args.pep
outfile=args.out

basename = os.path.basename(pepfile)
chro = os.path.splitext(basename)[0]

start=1
gap=2000
with open(outfile,"w")as f:
    for index,record in enumerate(SeqIO.parse(pepfile, "fasta")):
        i=index+1
        genename=record.id
        end=start+gap
        f.write(f"{chro}\t{start}\t{end}\t{genename}\t.\t+\n")
        start=end+1
