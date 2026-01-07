#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: start.py
# Created time: 2025/06/26 21:29:16
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python start.py -h


helptext=r"""
bed:
A01	358	1313	A01p00010.1_BnaDAR
A01	1351	2505	A01p00020.1_BnaDAR
A01	3679	5038	A01p00030.1_BnaDAR

input:
BnaA01T0000100ZS	ZS11A01G000110	A01p00080.1_BnaDAR	A01p000100.1_BnaEXP	BnaA01T0001200GG	-	A01g000013	Bnascaffold2421T0001200NO	Bnascaffold2787T0075200QU	BnaA01T0001300SL,Bnascaffold3315T0007000SL	Bnascaffold1809T0000600TA,Bnascaffold2359T0000700TA	Bnascaffold1716T0018400WE	BnaA01T0001600ZY,Bnascaffold2414T0000800ZY	ZY821000048
BnaA01T0000200ZS	ZS11A01G000100	A01p00060.1_BnaDAR	A01p000090.1_BnaEXP	BnaA01T0001100GG	-	-	Bnascaffold2421T0001000NO	Bnascaffold2787T0075300QU,Bnascaffold520T0001000QU	BnaA01T0001200SL	Bnascaffold1809T0000500TA	Bnascaffold1716T0018500WE,Bnascaffold286T0000600WE	BnaA01T0001500ZY	ZY821000047
BnaA01T0000300ZS	-	-	-	-	-	-	-	-	-	-	-	-	-

vote:
A01
A01
A01

out:
124814	71724	45758	65827	64819	-	110119	-	-	229359	-	-	97535	227663
152729	48475	14808	35498	45464	-	-	-	-	198380	-	-	78224	208551
156327	-	-	-	-	-	-	-	-	-	-	-	-	-

"""


import os
import argparse
parser=argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-b","--beddir",required=True,help="bed文件夹路径")
parser.add_argument("-i","--input","--geneindex",required=True,help="geneindex 文件(全是gene名称)")
parser.add_argument("-o","--out",required=True,help="输出文件")
parser.add_argument("-v","--vote",required=True,help="vote文件")
parser.add_argument("-a","--all",help="默认展示所有 start",action="store_true")
args=parser.parse_args()

votefile=args.vote
SG_file=args.input
beddir=args.beddir
outfile=args.out
allstart=args.all


class Gene:
    def __init__(self,name,chro,start,end):
        self.name=name
        self.chro=chro
        self.start=start
        self.end=end

    def __str__(self):
        return f"{self.chro};{self.name};range:{self.start}-{self.end}"


index2votechro={}
for index,chro in enumerate(open(votefile)):
    chro=chro.strip()
    index2votechro[index]=chro

gene2annotation={}
for bedfile in os.listdir(beddir):
    if bedfile.endswith(".bed"):
        bedpath=os.path.join(beddir,bedfile)
        with open(bedpath) as f:
            for line in f:
                line=line.strip().split("\t")
                name=line[3]
                chro=line[0]
                start=line[1]
                end=line[2]
                gene2annotation[name]=Gene(name,chro,start,end)


out=open(outfile,"w")
for index,i in enumerate(open(SG_file)):
    i=i.strip().split("\t")
    c=[]
    for j in i:
        if j=="-":
            b="-"
        elif "," in j:
            a=[]
            anogap=[]
            for k in j.split(","):
                annotation=gene2annotation[k]
                thischro=str(annotation.chro)
                votechro=index2votechro[index]
                if thischro==votechro:
                    a.append(str(annotation.start))
                    anogap.append(annotation.start)
                else:
                    a.append("-")
            if allstart:
                b=",".join(a)
            else:
                if len(anogap)==0:
                    b="-"
                else:
                    b=str(min(anogap))
        else:
            annotation=gene2annotation[j]
            thischro=str(annotation.chro)
            votechro=index2votechro[index]
            if thischro==votechro:
                b=str(annotation.start)
            else:
                b="-"
        c.append(b)
    out.write("\t".join(c)+"\n")

out.close()
