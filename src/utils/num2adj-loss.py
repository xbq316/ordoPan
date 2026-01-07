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
parser.add_argument("-n", "--num", required=True,help="未排序的num文件")
parser.add_argument("-s", "--sortnum", required=True,help="已排序的num文件")
args = parser.parse_args()

num_file=args.num
sort_num_file=args.sortnum


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

c=[]
for i in open(sort_num_file):
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
    c.append(b)
c=list(zip(*c))


T=set();T_=set() # 用于计算邻接相似度
G=set();G_=set() # 用于计算基因缺失率
for species_index,i in enumerate(a):
    species_name=f"spe_{species_index}"
    b=[]
    for j in i:
        if j!=-9:
            b.append(j)
    for j in b:
        G.add(frozenset([species_name,j]))
    
    b.sort()
    for j in range(len(b)-1):
        T_.add(frozenset([species_name,b[j],b[j+1]]))

for species_index,i in enumerate(c):
    species_name=f"spe_{species_index}"
    b=[]
    for j in i:
        if j!=-9:
            b.append(j)
    for j in range(len(b)-1):
        T.add(frozenset([species_name,b[j],b[j+1]]))
    
    for j in b:
        G_.add(frozenset([species_name,j]))

    # b.sort()
    # for j in range(len(b)-1):
    #     T_.add(frozenset([species_name,b[j],b[j+1]]))



S=len(T&T_)/len(T|T_)

in_G_=set()
not_in_G_=set()
for i in G:
    if i in G_:
        in_G_.add(i)
    else:
        not_in_G_.add(i)
L=len(not_in_G_)/len(G)
print(S,L)
