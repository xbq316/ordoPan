#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: num2accuracy.py
# Created time: 2025/07/29 15:29:25
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python num2accuracy.py -h



import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-n", "--num", required=True,help="排序后的num文件")
args = parser.parse_args()
num_file=args.num


a=[]
for i in open(num_file):
    i=i.strip()
    if not i:
        continue
    line=i.split("\t")
    b=[]
    for j in line:
        if j=="-":
            b.append(-9)
        else:
            b.append(int(j))
    a.append(b)
a=list(zip(*a))

x=[]
for species_index,i in enumerate(a):
    species_name=f"spe_{species_index}"
    b=[]
    for j in i:
        if j!=-9:
            b.append(j)
    if len(b)<=100:
        continue
    T=set()
    for j in range(len(b)-1):
        T.add(frozenset([b[j],b[j+1]]))

    b.sort()
    T_=set()
    for j in range(len(b)-1):
        T_.add(frozenset([b[j],b[j+1]]))

    # 邻接相似度
    S=len(T&T_)/len(T|T_)
    x.append(S)


z=sum(x)/len(x)
print(S)
#print(f"输入文件{num_file} 邻接相似度:{S*100}%")
