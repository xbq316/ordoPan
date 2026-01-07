#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: sort2cluster3.py
# Created time: 2025/08/17 16:58:07
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python sort2cluster.py -h


from cluster import *
import sys

helptext="""
Convert the clusterSort-sorted geneindex file into a cluster order file.
"""

# Global setting: stdout line buffer (writes to the file immediately upon encountering a newline character)
sys.stdout.reconfigure(line_buffering=True)


def sort2orderClusters(clusters: list[Cluster],geneindex: list[list[int]])->list[Cluster]:
    line2index={} #tuple(line):int(index)
    for index,cluster in enumerate(clusters):
        for line in cluster.content:
            line2index[tuple(line)]=index

    for line in geneindex:
        if tuple(line) not in line2index:
            raise Exception(f"sorted geneindex line {line} not in clusters")
    
    firstline=geneindex[0]
    indexs=[line2index[tuple(firstline)]]
    for line in geneindex[1:]:
        thisindex=line2index[tuple(line)]
        if thisindex == indexs[-1]:
            continue
        else:
            indexs.append(thisindex)
    
    if len(indexs)!=len(clusters):
        raise Exception(f"sorted geneindex clusters num {len(geneindex)} not equal to original clusters num {len(clusters)}")
    
    order_clusters=[]
    for index in indexs:
        order_clusters.append(clusters[index])
    return order_clusters

def reCluster(clusters: list[Cluster], rate:float=0.5)-> list[Cluster]:
    result=[]
    lastCluster=clusters[0]

    for thisCluster in clusters[1:]:
        link_rate=lastCluster.link_rate(thisCluster)
        if link_rate>=rate:
            lastCluster=lastCluster+thisCluster
        else:
            result.append(lastCluster)
            lastCluster=thisCluster
    result.append(lastCluster)

    return result


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cluster", type=str, help="Input (raw) cluster file", required=True)
    parser.add_argument("-g", "--geneindex", type=str, help="The geneindex file output by clusterSort", required=True)
    parser.add_argument("-o", "--out", type=str, help="The output cluster file (re-clustered cluster)", required=True)
    parser.add_argument("-r", "--rate", type=float, help="cut off threshold",default=0.5)
    args = parser.parse_args()

    originalClusterFile=args.cluster
    sortedClusterFile=args.geneindex
    outClusterFile=args.out
    rate=args.rate

    # Output parameters
    print("original cluster file:",originalClusterFile)
    print("sorted cluster file:",sortedClusterFile)
    print("output cluster file:",outClusterFile)
    print("...")

    print("read original clusters...")
    originalClusters=read_clusters(originalClusterFile)

    print("read sorted clusters...")
    geneindex=readNum(sortedClusterFile)

    print("recluster...")
    order_clusters=sort2orderClusters(originalClusters,geneindex)
    
    re_clusters=reCluster(order_clusters,rate)
    writeClusterList(re_clusters,outClusterFile)

    print("done")
