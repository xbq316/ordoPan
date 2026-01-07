#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: kendallTau.py
# Created time: 2025/12/25 16:19:07
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python kendallTau.py -h



import sys
from pathlib import Path

parent_dir = str(Path(__file__).resolve().parent.parent)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

import cluster
if len(sys.argv)!=2:
    print("Usage: python kendallTau.py <geneindex num file>")
    sys.exit(1)
geneindexFile=sys.argv[1]


geneindex_num=cluster.readNum(geneindexFile)
repo=cluster.geneindex_permutation_report(geneindex_num,alpha=0.05)
print(f"{'*'*20} Processing $filename {'*'*20}")
print(repo)

