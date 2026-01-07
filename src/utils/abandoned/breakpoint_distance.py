#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: 03.breakpoint_distance.py
# Created time: 2025/06/03 16:02:40
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python 03.breakpoint_distance.py -h


from typing import List

def breakpoint_distance(reference_order:List[str], pred_order:List[str])->float:
    refset=set(reference_order)
    pred_order=[i for i in pred_order if i in refset]
    pred_order=["IMITATE_STRAT_GENE",*pred_order,"IMITATE_END_GENE"]
    reference_order=["IMITATE_STRAT_GENE",*reference_order,"IMITATE_END_GENE"]

    adjacant_set=set()
    for i in range(len(reference_order)-1):
        adjacant_set.add((reference_order[i],reference_order[i+1]))
    
    pred_adjacant_num=0
    true_adjacant_num=len(adjacant_set)
    for i in range(len(pred_order)-1):
        if (pred_order[i],pred_order[i+1]) in adjacant_set:
            pred_adjacant_num+=1
        elif (pred_order[i+1],pred_order[i]) in adjacant_set:
            pred_adjacant_num+=1
    
    TP=pred_adjacant_num#预测正确的相邻对数目
    TN=len(reference_order) - 1 - true_adjacant_num
    FN=true_adjacant_num - pred_adjacant_num#在参考中存在但预测中缺失的相邻对
    FP=len(pred_order) - 1 - pred_adjacant_num
    print(TP,FN,FP,TN)

    return {"TP":TP,"FN":FN,"FP":FP,"TN":TN}

def Accuracy_Recall_F1(dict_data:dict) -> str:
    precision=dict_data["TP"]/(dict_data["TP"]+dict_data["FP"])
    recall=dict_data["TP"]/(dict_data["TP"]+dict_data["FN"])
    f1=2*precision*recall/(precision+recall)
    return " ".join([f"{precision:.4f}",f"{recall:.4f}",f"{f1:.4f}"])

if __name__=="__main__":
    import argparse

    parser=argparse.ArgumentParser("计算断点距离")
    parser.add_argument("-q","--query",type=str,help="query file(gene list)",required=True)
    parser.add_argument("-a","--ancestor",type=str,help="ancestor file(gene list)",required=True)
    args=parser.parse_args()

    queryfile=args.query
    ancestorfile=args.ancestor


    true_order = []
    pred_order =[]

    with open(queryfile,"r") as f:
        for line in f:
            pred_order.append(line.strip())
            
    with open(ancestorfile,"r") as f:
        for line in f:
            true_order.append(line.strip())

    #print(Accuracy_Recall_F1(breakpoint_distance(true_order, pred_order)))
    breakpoint_distance(true_order, pred_order)
