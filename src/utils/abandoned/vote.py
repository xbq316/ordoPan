#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: vote.py
# Created time: 2025/06/26 21:17:10
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python vote.py -h



helptext=r"""
输入文件格式:
A01	A01	A01	A01	A01	-	A01	scaffold2421	scaffold2787	A01,scaffold3315	scaffold1809,scaffold2359	scaffold1716	A01,scaffold2414	A01
A01	A01	A01	A01	A01	-	-	scaffold2421	scaffold2787,scaffold520	A01	scaffold1809	scaffold1716,scaffold286	A01	A01
A01	-	-	-	-	-	-	-	-	-	-	-	-	-

输出文件格式：
A01
A01
A01
"""

import argparse
parser=argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i","--input",required=True,type=str,help="输入文件")
parser.add_argument("-o","--out",required=True,type=str,help="输出文件")
args=parser.parse_args()

infile=args.input
outfile=args.out


with open(outfile,"w") as f:
    for i in open(infile):
        a={}
        i=i.strip().split("\t")
        for j in i:
            if j =="-":
                continue
            elif "," in j:
                for x in j.split(","):
                    a[x]=a.get(x,0)+1
            else:
                a[j]=a.get(j,0)+1
        bestchro=max(a.items(),key=lambda x:x[1])[0]
        f.write(f"{bestchro}\n")


