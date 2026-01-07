#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: annotation.py
# Created time: 2025/06/27 16:29:50
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python annotation.py -h


import argparse
import os



helptext=r"""
输入文件格式:

geneindex:
BnaA01T0000100ZS	ZS11A01G000110	A01p00080.1_BnaDAR	A01p000100.1_BnaEXP	BnaA01T0001200GG	-	A01g000013	Bnascaffold2421T0001200NO	Bnascaffold2787T0075200QU	BnaA01T0001300SL,Bnascaffold3315T0007000SL	Bnascaffold1809T0000600TA,Bnascaffold2359T0000700TA	Bnascaffold1716T0018400WE	BnaA01T0001600ZY,Bnascaffold2414T0000800ZY	ZY821000048
BnaA01T0000200ZS	ZS11A01G000100	A01p00060.1_BnaDAR	A01p000090.1_BnaEXP	BnaA01T0001100GG	-	-	Bnascaffold2421T0001000NO	Bnascaffold2787T0075300QU,Bnascaffold520T0001000QU	BnaA01T0001200SL	Bnascaffold1809T0000500TA	Bnascaffold1716T0018500WE,Bnascaffold286T0000600WE	BnaA01T0001500ZY	ZY821000047
BnaA01T0000300ZS	-	-	-	-	-	-	-	-	-	-	-	-	-
BnaA01T0000400ZS	ZS11A01G000090	A01p00050.1_BnaDAR	A01p000080.1_BnaEXP	BnaA01T0001000GG	-	A01g000011	Bnascaffold2421T0000900NO	Bnascaffold2787T0075400QU,Bnascaffold520T0000900QU	BnaA01T0001100SL	Bnascaffold1809T0000400TA	Bnascaffold1716T0018700WE,Bnascaffold286T0000500WE	BnaA01T0001400ZY	ZY821000046

bed:
A01	358	1313	A01p00010.1_BnaDAR
A01	1351	2505	A01p00020.1_BnaDAR
A01	3679	5038	A01p00030.1_BnaDAR
A01	9617	13540	A01p00040.1_BnaDAR

out 输出两个文件
annotation:
A01;BnaA01T0000100ZS;range:124814-125347;strand:+	A01;ZS11A01G000110;range:71724-72409;strand:-	A01;A01p00080.1_BnaDAR;range:45758-46339;strand:-	A01;A01p000100.1_BnaEXP;range:65827-66408;strand:-	A01;BnaA01T0001200GG;range:64819-65423;strand:-	-	A01;A01g000013;range:110119-110700;strand:-	scaffold2421;Bnascaffold2421T0001200NO;range:71789-72412;strand:-	scaffold2787;Bnascaffold2787T0075200QU;range:4524356-4524889;strand:+	A01;BnaA01T0001300SL;range:229359-229892;strand:-,scaffold3315;Bnascaffold3315T0007000SL;range:366565-367098;strand:+	scaffold1809;Bnascaffold1809T0000600TA;range:51425-52005;strand:-,scaffold2359;Bnascaffold2359T0000700TA;range:71052-71632;strand:-	scaffold1716;Bnascaffold1716T0018400WE;range:1078408-1078941;strand:+	A01;BnaA01T0001600ZY;range:97535-98144;strand:-,scaffold2414;Bnascaffold2414T0000800ZY;range:109372-109981;strand:-	A01;ZY821000048;range:227663-228228;strand:-
A01;BnaA01T0000200ZS;range:152729-153872;strand:+	A01;ZS11A01G000100;range:48475-50775;strand:-	A01;A01p00060.1_BnaDAR;range:14808-16914;strand:-	A01;A01p000090.1_BnaEXP;range:35498-36899;strand:-	A01;BnaA01T0001100GG;range:45464-46933;strand:-	-	-	scaffold2421;Bnascaffold2421T0001000NO;range:42916-44390;strand:-	scaffold2787;Bnascaffold2787T0075300QU;range:4550986-4552262;strand:+,scaffold520;Bnascaffold520T0001000QU;range:27719-29063;strand:-	A01;BnaA01T0001200SL;range:198380-200427;strand:-	scaffold1809;Bnascaffold1809T0000500TA;range:14582-16050;strand:-	scaffold1716;Bnascaffold1716T0018500WE;range:1183601-1184848;strand:+,scaffold286;Bnascaffold286T0000600WE;range:22727-28558;strand:-	A01;BnaA01T0001500ZY;range:78224-79694;strand:-	A01;ZY821000047;range:208551-209734;strand:-
A01;BnaA01T0000300ZS;range:156327-156924;strand:-	-	-	-	-	-	-	-	-	-	-	-	-	-

chro:
A01	A01	A01	A01	A01	-	A01	scaffold2421	scaffold2787	A01,scaffold3315	scaffold1809,scaffold2359	scaffold1716	A01,scaffold2414	A01
A01	A01	A01	A01	A01	-	-	scaffold2421	scaffold2787,scaffold520	A01	scaffold1809	scaffold1716,scaffold286	A01	A01
A01	-	-	-	-	-	-	-	-	-	-	-	-	-

"""
parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-b","--beddir",required=True,help="bed文件夹路径")
parser.add_argument("-i","--input","--geneindex",required=True,help="geneindex 文件(全是gene名称)")
parser.add_argument("-a","--annotation",default="all",choices=["gene","all"],help="输出的注释使用gene名称还是所有信息")
parser.add_argument("-o","--out",required=True,help="输出目录")
args=parser.parse_args()

beddir=args.beddir
geneindex=args.input
outdir=args.out
annotation=args.annotation

class Gene:
    def __init__(self,name,chro,start,end):
        self.name=name
        self.chro=chro
        self.start=start
        self.end=end
    def __str__(self):
        return f"{self.chro};{self.name};range:{self.start}-{self.end}"



gene2annotation={}
for bedfile in os.listdir(beddir):
    if bedfile.endswith(".bed"):
        bedpath=os.path.join(beddir,bedfile)
        with open(bedpath) as f:
            if annotation=="all":
                for line in f:
                    line=line.strip().split("\t")
                    name=line[3]
                    chro=line[0]
                    start=line[1]
                    end=line[2]
                    gene2annotation[name]=Gene(name,chro,start,end)
            elif annotation=="gene":
                for line in f:
                    line=line.strip().split("\t")
                    name=line[3]
                    chro=line[0]
                    start=line[1]
                    end=line[2]
                    gene_name=name.split("_")[0]
                    gene2annotation[gene_name]=gene_name


geneindex_chro_file=os.path.join(outdir,"geneindex.chro.txt")
geneindex_annotation_file=os.path.join(outdir,"geneindex.annotation.txt")

geneindex_chro=open(geneindex_chro_file,"w")
geneindex_annotation=open(geneindex_annotation_file,"w")


for i in open(geneindex):
    i=i.strip().split("\t")
    c1=[]
    c2=[]
    for j in i:
        if j=="-":
            b1=b2="-"
        elif "," in j:
            a1=[]
            a2=[]
            for k in j.split(","):
                a1.append(str(gene2annotation[k].chro))
                a2.append(str(gene2annotation[k]))
            b1=",".join(a1)
            b2=",".join(a2)
        else:
            b1=str(gene2annotation[j].chro)
            b2=str(gene2annotation[j])
        c1.append(b1)
        c2.append(b2)
    geneindex_chro.write("\t".join(c1)+"\n")
    geneindex_annotation.write("\t".join(c2)+"\n")

geneindex_chro.close()
geneindex_annotation.close()
