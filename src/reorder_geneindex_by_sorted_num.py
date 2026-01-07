#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: resort.py
# Created time: 2025/06/27 17:25:06
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python resort.py -h

import cluster

def reSortGeneIndex(numFile,sortedNumFile,geneindexFile)->list[list[str]]:
    """
    Reorder the geneindex file based on the unsorted and sorted num files
    numFile: Unsorted num file
    sortedNumFile: Sorted num file
    geneindexFile: Unsorted geneindex file
    outFile: Output file containing the sorted geneindex
    """

    geneNum2index={}
    for index,i in enumerate(cluster.readNum(numFile)):
        geneNum2index[tuple(i)]=index
    path=[]
    for i in cluster.readNum(sortedNumFile):
        try:
            index=geneNum2index[tuple(i)]
        except KeyError:
            raise Exception(f"sorted num file {sortedNumFile} line {i} not in num {numFile}file")
        path.append(index)

    geneindex=cluster.read_geneindex(geneindexFile)

    result=[geneindex[index] for index in path]
    return result


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("-n","--num",required=True,help="Unsorted num files")
    parser.add_argument("-s","--sortednum",required=True,help="Sorted num file")
    parser.add_argument("-g","--geneindex",required=True,help="Unsorted geneindex file")
    parser.add_argument("-o","--out",required=True,help="Output file: Sorted geneindex file")
    args=parser.parse_args()

    originalNumFile=args.num
    sortedNumFile=args.sortednum
    geneindexFile=args.geneindex
    outFile=args.out

    with open(outFile,"w") as f:
        for line in reSortGeneIndex(originalNumFile,sortedNumFile,geneindexFile):
            f.write("\t".join(line)+"\n")

