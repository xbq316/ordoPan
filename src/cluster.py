#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: cluster.py
# Created time: 2025/06/27 17:27:41
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python cluster.py -h


import networkx as nx
import tqdm
import math
from scipy.stats import kendalltau
import time
import os
from typing import List
import random
from collections import defaultdict
import subprocess
from pathlib import Path

def format_duration(seconds: float) -> str:
    """格式化耗时：秒 / 分钟:秒 / 小时:分钟:秒"""
    if seconds < 60:
        return f"{seconds:.2f} 秒"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        sec = seconds % 60
        return f"{minutes} 分 {sec:.2f} 秒"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        sec = seconds % 60
        return f"{hours} 小时 {minutes} 分 {sec:.2f} 秒"

def readNum(file:str)-> list[list[int]]:
    result=[]
    for i in open(file):
        i=i.strip()
        if i:
            line=[int(x) if x!="-" else -9 for x in i.split("\t")]
            result.append(line)
    return result

def read_bed4(file:str)-> tuple[dict[str,int],dict[str,str]]:
    '''
    A01	831	1437	BnaA01g00010D
    '''
    gene2start={}
    gene2chromosome={}
    for i in open(file):
        line=i.strip().split()
        gene2start[line[3]]=int(line[1])
        gene2chromosome[line[3]]=line[0]
    return gene2start,gene2chromosome

def read_geneindex(file:str)-> list[list[str]]:
    result=[]
    for i in open(file):
        i=i.strip()
        if i:
            line=["-" if not x else x for x in i.split("\t") ]
            result.append(line)
    return result

class Gene:
    def __init__(self,chromosome: str,start: int,end: int,name: str,score:str,strand: str):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand

    def __str__(self):
        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.name}\t{self.score}\t{self.strand}"

def read_bed6(file:str)->list[Gene]:
    """
    chromosome\tstart\tend\name\tscore\tstrand
    """
    result=[]
    with open(file) as bed6:
        for i in bed6:
            i=i.strip()
            if not i:
                continue
            line=i.split("\t")
            try:
                chromosome,start,end,name,score,strand=line
                start=int(start)
                end=int(end)
            except :
                raise ValueError(f"{line} format error")
            result.append(Gene(chromosome,start,end,name,score,strand))
    return result

def bed6_resultTobed4_result(bed6: list[Gene])->tuple[dict[str,int],dict[str,str]]:
    gene2start={}
    gene2chromosome={}
    for gene in bed6:
        gene2start[gene.name]=gene.start
        gene2chromosome[gene.name]=gene.chromosome
    return gene2start,gene2chromosome

def read_sample(file:str)->list[str]:
    result=[]
    with open(file) as sample:
        for i in sample:
            i=i.strip()
            if not i:
                continue
            result.append(i)
    return result
    
def pan2index(geneindex: list[list[str]],beddir:str,sample_list:list[str],bestchromosome:str,bed_type:str)-> tuple[list[list[str]],list[list[str]]]:
    """
    将当前染色体构建的 geneindex 的基因转为 gene在该染色体上的序号
    return:
        1: 这条染色体上(geneindex)的每一个基因顺序
        2: 这条染色体上(geneindex)的每一个单元格顺序
    """

    if bed_type not in ["bed6","bed4"]:
        raise Exception("bed_type must be bed6 or bed4")
        
    # 检查输入文件
    for this_sample in sample_list:
        bed_file=os.path.join(beddir,f"{this_sample}.bed")
        if not os.path.exists(bed_file):
            raise FileNotFoundError(f"{bed_file} does not exist")

    geneindex_list=list(zip(*geneindex))
    if geneindex_list.__len__()!=sample_list.__len__():
        raise Exception(f"number of samples {sample_list.__len__()} not equal to number of geneindex columns {geneindex_list.__len__()}")


    # 先找到最佳染色体
    if bestchromosome:
        bestChromosome=bestchromosome
    else:
        chromosome_list=[]
        for sample_index,sample in enumerate(sample_list):
            bed_file=os.path.join(beddir,f"{sample}.bed")
            if bed_type=="bed6":
                bed6=read_bed6(bed_file)
                gene2start,gene2chromosome=bed6_resultTobed4_result(bed6)

            elif bed_type=="bed4":
                gene2start,gene2chromosome=read_bed4(bed_file)
            
            for genes in geneindex_list[sample_index]:
                try:
                    if genes == "-":
                        continue
                    elif "," in genes:
                        for gene in genes.split(","):
                            chromosome_list.append(gene2chromosome[gene])
                    else:
                        chromosome_list.append(gene2chromosome[genes])
                except KeyError:
                    raise KeyError(f"gene {genes} not found in {bed_file}")
        bestChromosome=chromosome_list.max(set(chromosome_list),key=chromosome_list.count)

    result1=[]
    result2=[]
    for sample_index,sample in enumerate(sample_list):
        bed_file=os.path.join(beddir,f"{sample}.bed")
        
        if bed_type=="bed6":
            bed6=read_bed6(bed_file)
            gene2start,gene2chromosome=bed6_resultTobed4_result(bed6)
        elif bed_type=="bed4":
            gene2start,gene2chromosome=read_bed4(bed_file)

        startIndexGene_list=[]
        minstartIndex_list=[]
        for thisindex,genes in enumerate(geneindex_list[sample_index]):
            try:
                if genes == "-":
                    continue
                elif "," in genes:
                    starts=[]
                    for gene in genes.split(","):
                        this_chro=gene2chromosome[gene]
                        if this_chro!=bestChromosome:
                            continue
                        startIndexGene_list.append([gene2start[gene],gene])
                        starts.append(gene2start[gene])
                    if not starts:
                        continue
                    minstart=min(starts)
                    minstartIndex_list.append([minstart,thisindex])
                else:
                    this_chro=gene2chromosome[genes]
                    if this_chro!=bestChromosome:
                        continue
                    startIndexGene_list.append([gene2start[genes],genes])
                    minstartIndex_list.append([gene2start[genes],thisindex])
            except KeyError:
                raise KeyError(f"gene {genes} not found in {bed_file}")
        startIndexGene_list.sort(key=lambda x: x[0])
        gene2order={i[1]:order for order,i in enumerate(startIndexGene_list)}
        minstartIndex_list.sort(key=lambda x: x[0])
        index2order={i[1]:order for order,i in enumerate(minstartIndex_list)}

        a1=[]
        a2=[]
        for thisindex,genes in enumerate(geneindex_list[sample_index]):
            try:
                if genes == "-":
                    a1.append("-")
                    a2.append("-")
                elif "," in genes:
                    if thisindex in index2order:
                        a2.append(index2order[thisindex])
                    else:
                        a2.append("-")
                    orders=[]
                    for gene in genes.split(","):
                        this_chro=gene2chromosome[gene]
                        if this_chro!=bestChromosome:
                            continue
                        orders.append(str(gene2order[gene]))
                    a1.append(",".join(orders))
                else:
                    this_chro=gene2chromosome[genes]
                    if this_chro!=bestChromosome:
                        continue
                    a2.append(index2order[thisindex])
                    a1.append(gene2order[genes])
            except KeyError:
                raise KeyError(f"gene {genes} not found in {bed_file}")
        result1.append(a1)
        result2.append(a2)
    return result1,result2
    
    

class Cluster:
    def __init__(self,line:list[int]):
        self.start=line
        self.end=line
        self.speciesNum=len(line)
        self.content=[line]

    @classmethod
    def frommultiline(cls, lines:list[list[int]]) -> "Cluster":    
        if not lines:
            raise ValueError("Input lines cannot be empty")
        newCluster = cls([])
        newCluster.content = lines
        newCluster.start=newCluster.updateStart()
        newCluster.end=newCluster.updateEnd()
        newCluster.speciesNum = len(lines[0]) if lines else 0
        return newCluster

    def updateStart(self):
        newLine=self.content[0]
        for line in self.content:
            newLine=[i if i!=-9 else j for i,j in zip(newLine,line)]   
        return newLine
    
    def updateEnd(self):
        newLine=self.content[-1]
        for line in self.content[::-1]:
            newLine=[i if i!=-9 else j for i,j in zip(newLine,line)]   
        return newLine

    def __add__(self,nextCluster:"Cluster")-> "Cluster":
        if not isinstance(nextCluster, Cluster):
            raise TypeError("Can only add another Cluster instance")
        newCluster=Cluster.frommultiline(self.content+nextCluster.content)
        return newCluster

    def cut_off(self,other:"Cluster")-> float:
        """
        计算cutoff值，基于三个因素：
        1. 邻接密度：t/min(m,n)
        2. 基因平衡：(2*sqrt(m*n))/(m+n)
        3. 归一化邻接数：log(t+1)/log(global_max_t+1)
        """
        if self.speciesNum != other.speciesNum:
            raise ValueError("Clusters must have the same number of species")
        adjacency=0
        for i,j in zip(self.end, other.start):
            if i==-9 or j==-9:
                continue
            adjacency+=1 if i+1==j else 0
        m=self.suportedSpecies()
        n=other.suportedSpecies()
        # 使用min(m,n)作为分母
        adjacency_density = adjacency / min(m, n)
        # 基因平衡因子，范围[0,1]
        gene_balance = (2 * math.sqrt(m * n)) / (m + n) if (m + n) > 0 else 0.0
        # 归一化邻接数，确保global_max_t至少为1
        norm_t = math.log(adjacency + 1) / math.log(max(2, self.speciesNum))
        # 返回三个因子的乘积
        return adjacency_density * gene_balance * norm_t

    def __str__(self):
        return f"content len:{len(self.content)}"

    def isLinked(self, nextCluster:"Cluster")-> bool:
        if not isinstance(nextCluster, Cluster):
            raise TypeError("Can only check linkage with next Cluster instance")
        for lastEnd,thisStart in zip(self.end, nextCluster.start):
            if lastEnd==-9 or thisStart==-9:
                continue
            if lastEnd+1==thisStart:
                return True #相邻
        return False #不相邻
    def suportedSpecies(self) -> int:
        """返回支持的物种数量"""
        return self.speciesNum-self.start.count(-9)

    def link_num(self,nextCluster:"Cluster")-> int:
        result=0
        if not isinstance(nextCluster, Cluster):
            raise TypeError("Can only check linkage with next Cluster instance")
        for lastEnd,thisStart in zip(self.end, nextCluster.start):
            if lastEnd==-9 or thisStart==-9:
                continue
            if lastEnd+1==thisStart:
                result+=1
        return result

    def unlink_num(self,nextCluster:"Cluster")-> int:
        result=0
        if not isinstance(nextCluster, Cluster):
            raise TypeError("Can only check linkage with next Cluster instance")
        for lastEnd,thisStart in zip(self.end, nextCluster.start):
            if lastEnd==-9 or thisStart==-9:
                continue
            if lastEnd+1!=thisStart:
                result+=1
        return result
        

    def link_rate(self,nextCluster:"Cluster")-> float:
        """
        计算相邻率:
        邻接数/(能进行比较的物种数)
        能进行比较的物种数：都不是 缺失 就+1
        """
        if not isinstance(nextCluster, Cluster):
            raise TypeError("Can only check linkage with next Cluster instance")
        norm_num=0
        adjacency=0
        for i,j in zip(self.end, nextCluster.start):
            if i==-9 or j==-9:
                continue
            norm_num+=1
            if abs(i-j)==1:
                adjacency+=1
        return 0 if norm_num==0 else adjacency/norm_num

def read_clusters(cluster_file:str)->list[Cluster]:
    """
    读取 Cluster 文件
    
    :param cluster_file: Cluster 文件名
    :return: 一个 list，包含了所有 Cluster 对象
    """
    result=[]
    cluster=[]

    for i in open(cluster_file):
        i=i.strip()
        if i.startswith("#Cluster-"):
            if cluster:
                result.append(Cluster.frommultiline(cluster))
            cluster=[]
        else:
            line=[int(x) if x!="-" else -9 for x in i.split()]
            cluster.append(line)

    if cluster:
        result.append(Cluster.frommultiline(cluster))
    return result

def geneindex2clusters(geneindex:list[list[int]])-> list[Cluster]:
    """
    将基因索引转换为 Cluster 列表
    :param geneindex: 基因索引列表
    :return: Cluster 列表
    """
    clusters=[]
    for line in geneindex:
        clusters.append(Cluster(line))
    return clusters
    
def geneindex2rankScore(geneindex:list[list[int]],Minimum_genome_length=100)-> float:
    newgeneindex=list(zip(*geneindex))
    ranks=[]
    for genome in newgeneindex:
        a=[i for i in genome if i!=-9]
        if len(a) < Minimum_genome_length:
            continue
        predict_ancestor=list(range(len(a)))

        index_gene=list(enumerate(a))
        index_gene.sort(key=lambda x: x[1])
        current_genome=[x[0] for x in index_gene]
        this_rank,_=kendalltau(predict_ancestor,current_genome)
        ranks.append(abs(this_rank))
    return sum(ranks)/len(ranks) if ranks else 0.0


def randome_combinations(this_list: List[int], n: int, pair_size: int = 2) -> dict:
    """
    从列表中随机抽取 n 次组合，每个组合大小为 pair_size，
    使用 frozenset 作为 key 返回组合出现次数。
    
    参数:
        this_list: 原始列表
        n: 抽取次数
        pair_size: 每次组合大小，默认2
        
    返回:
        dict: {frozenset(组合): 出现次数}
    """
    pair_counts = defaultdict(int)
    
    for _ in range(n):
        pair = random.sample(this_list, pair_size)  # 无放回选择 pair_size 个元素
        key = frozenset(pair)                       # 转为 frozenset 作为 key
        pair_counts[key] += 1
    return dict(pair_counts)


from typing import List
from collections import defaultdict
from scipy.stats import kendalltau
import random
import numpy as np

def randome_combinations(this_list: List[int], n: int, pair_size: int = 2) -> dict:
    pair_counts = defaultdict(int)
    for _ in range(n):
        pair = random.sample(this_list, pair_size)
        key = frozenset(pair)
        pair_counts[key] += 1
    return dict(pair_counts)


def geneindex_permutation_report(
    geneindex: List[List[int]],
    Minimum_genome_length: int = 100,
    random_combinations: int = 100,
    n_permutations: int = 1000,
    alpha: float = 0.05
) -> str:
    """
    对 geneindex 的两个群体进行置换检验，返回可直接报告的字符串。
    单尾检验方向: mean(ancestor) > mean(random)
    """
    # 转置基因矩阵
    newgeneindex = list(zip(*geneindex))
    n_genomes = len(newgeneindex)
    
    # 随机组合
    pairs2count = randome_combinations(list(range(n_genomes)), random_combinations, pair_size=2)
    
    pair2legal = defaultdict(bool)
    pair2kendalltau = defaultdict(float)
    
    # 计算随机组合的 Kendall tau
    for pair, count in pairs2count.items():
        x = []
        genome1 = newgeneindex[list(pair)[0]]
        genome2 = newgeneindex[list(pair)[1]]
        for index, (gene1, gene2) in enumerate(zip(genome1, genome2)):
            if gene1 == -9 or gene2 == -9:
                continue
            x.append((gene1, gene2, index))
        if len(x) < Minimum_genome_length:
            pair2legal[pair] = False
            continue
        x.sort(key=lambda item: item[0])
        x_1 = [sd[2] for sd in x]
        x.sort(key=lambda item: item[1])
        x_2 = [sd[2] for sd in x]
        score, _ = kendalltau(x_1, x_2)
        pair2kendalltau[pair] = abs(score)
        pair2legal[pair] = True
    
    original_kendalltau = []
    for pair in pairs2count:
        if pair2legal[pair]:
            count = pairs2count[pair]
            score = pair2kendalltau[pair]
            original_kendalltau.extend([score] * count)
    
    # 祖先预测群体
    ancestor_kendalltau = []
    for genome in newgeneindex:
        a = [i for i in genome if i != -9]
        if len(a) < Minimum_genome_length:
            continue
        predict_ancestor = list(range(len(a)))
        index_gene = list(enumerate(a))
        index_gene.sort(key=lambda x: x[1])
        current_genome = [x[0] for x in index_gene]
        this_kendalltau, _ = kendalltau(predict_ancestor, current_genome)
        ancestor_kendalltau.append(abs(this_kendalltau))
    
    # 计算均值和差异
    mean_random = np.mean(original_kendalltau)
    mean_ancestor = np.mean(ancestor_kendalltau)
    diff = mean_ancestor - mean_random
    
    # 合并两群体
    combined = original_kendalltau + ancestor_kendalltau
    n_orig = len(original_kendalltau)
    
    # 置换检验
    count_extreme_two_tail = 0
    count_extreme_one_tail = 0
    
    for _ in range(n_permutations):
        np.random.shuffle(combined)
        perm_orig = combined[:n_orig]
        perm_anc = combined[n_orig:]
        perm_diff = np.mean(perm_anc) - np.mean(perm_orig)
        
        if abs(perm_diff) >= abs(diff):
            count_extreme_two_tail += 1
        if perm_diff >= diff:
            count_extreme_one_tail += 1
    
    p_two_tail = count_extreme_two_tail / n_permutations
    p_one_tail = count_extreme_one_tail / n_permutations
    
    # 是否拒绝原假设
    # 判定
    reject_two_tail = p_two_tail < alpha
    reject_one_tail = p_one_tail < alpha

    # 检查效应方向
    direction_note = ""
    if diff < 0:
        direction_note = (
            "NOTE: Observed effect direction is ancestor < random, "
            "which is opposite to the one-tailed alternative hypothesis (ancestor > random). "
            "Thus the one-tailed p-value is expected to be large.\n"
        )
    elif diff > 0:
        direction_note = (
            "NOTE: Observed effect direction is ancestor > random, "
            "which matches the one-tailed alternative hypothesis.\n"
        )
    else:
        direction_note = (
            "NOTE: Observed effect shows no difference between groups.\n"
        )

# 生成报告字符串（包含明确定义的 H0 与 H1）
    report = (
        f"Permutation Test Report (α = {alpha})\n"
        f"-------------------------------------------\n"
        f"Randomized Group Mean Kendall tau: {mean_random:.4f}\n"
        f"Ancestor Predicted Group Mean Kendall tau: {mean_ancestor:.4f}\n"
        f"Observed Mean Difference (ancestor - random): {diff:.4f}\n"
        f"\n"
        f"Two-tailed Hypothesis Test:\n"
        f"  H0: mean_ancestor = mean_random\n"
        f"  H1: mean_ancestor ≠ mean_random\n"
        f"  p-value = {p_two_tail:.4f} → {'Reject H0' if reject_two_tail else 'Do not reject H0'}\n"
        f"\n"
        f"One-tailed Hypothesis Test (Specified Direction: ancestor > random):\n"
        f"  H0: mean_ancestor ≤ mean_random\n"
        f"  H1: mean_ancestor > mean_random\n"
        f"  p-value = {p_one_tail:.4f} → {'Reject H0' if reject_one_tail else 'Do not reject H0'}\n"
        f"\n"
        f"{direction_note}"
    )
    return report

    

def geneindex2lcsScore(geneindex:list[list[int]],Minimum_genome_length=100)-> float:
    newgeneindex=list(zip(*geneindex))
    LCS_scores=[]
    for genome in newgeneindex:
        a=[i for i in genome if i!=-9]
        if len(a) < Minimum_genome_length:
            continue
        n=len(a)
        dp=[[0]*(n+1) for _ in range(n+1)]
        for i in range(1,n+1):
            for j in range(1,n+1):
                if a[i-1]==a[j-1] and i!=j:
                    dp[i][j]=dp[i-1][j-1]+1
                else:
                    dp[i][j]=max(dp[i-1][j],dp[i][j-1])
        LCS_scores.append(dp[n][n]/n)
    return sum(LCS_scores)/len(LCS_scores) if LCS_scores else 0.0



def clusters2linkNum(clusters:list[Cluster])-> int:
    linkNum=0

    newClusters=[]
    for thisCluster in clusters:
        for line in thisCluster.content:
            newClusters.append(Cluster(line))
    
    lastCluster=newClusters[0]
    for cluster in newClusters[1:]:
        linkNum += lastCluster.link_num(cluster)
        lastCluster = lastCluster + cluster
    
    return linkNum

def mergeClusters(clusters:list[Cluster])-> list[Cluster]:
    #建图
    clustersNum= len(clusters)
    G = nx.DiGraph()
    for node in range(clustersNum):
        G.add_node(node)
    for i in range(clustersNum):
        for j in range(clustersNum):
            if i != j and clusters[i].isLinked(clusters[j]):
                G.add_edge(i, j)

    #对出度或入度为1的节点进行合并
    mergedClusters = []
    visited = set()

    for node in G.nodes():
        if node not in visited:
            in_degree = G.in_degree(node)
            out_degree = G.out_degree(node)
            if in_degree==1:
                pred_node = list(G.predecessors(node))[0]
                if node in visited or pred_node in visited:
                    continue
                mergedCluster=clusters[pred_node]+clusters[node]
                mergedClusters.append(mergedCluster)
                visited.add(node)
                visited.add(pred_node)
            elif out_degree==1:
                next_node=list(G.successors(node))[0]
                if node in visited or next_node in visited:
                    continue
                mergedCluster=clusters[node]+clusters[next_node]
                mergedClusters.append(mergedCluster)
                visited.add(node)
                visited.add(next_node)
    for node in G.nodes():
        if node not in visited:
            mergedClusters.append(clusters[node])
            visited.add(node)
    return mergedClusters

def writeClusterList(clusters:list[Cluster],outFile:str):
    with open(outFile, "w") as f:
        for index,cluster in enumerate(clusters):
            f.write(f"#Cluster-{index}:\n")
            for line in cluster.content:
                line=[x if x!=-9 else "-" for x in line]
                f.write("\t".join(str(x) for x in line) + "\n")

def writeGeneIndex(geneindex:list[list[str]],outFile:str):
    with open(outFile, "w") as f:
        for line in geneindex:
            line=[str(x) for x in line]
            f.write("\t".join(str(x) for x in line) + "\n")

def writeGeneIndexNum(geneindex:list[list[int]],outFile:str):
    with open(outFile, "w") as f:
        for line in geneindex:
            line=[str(x) if x!=-9 else "-" for x in line]
            f.write("\t".join(str(x) for x in line) + "\n")

def sort2cluster(originalClusterFile: str,
                 sortedClusterFile: str,
                 outClusterFile: str,
                 rate: float = 0.5):
    """
    Convert the clusterSort-sorted geneindex file into a cluster order file.
    cluster_file: str
        Input (raw) cluster file
    geneindex_file: str
        The geneindex file output by clusterSort
    output_cluster_file: str
        The output cluster file (re-clustered cluster)
    rate: float
        The linkage rate threshold for merging clusters
    """
    # Output parameters
    # print("original cluster file:",originalClusterFile)
    # print("sorted cluster file:",sortedClusterFile)
    # print("output cluster file:",outClusterFile)
    # print("...")
    # print("read original clusters...")
    originalClusters=read_clusters(originalClusterFile)
    # print("read sorted clusters...")
    geneindex=readNum(sortedClusterFile)
    # print("recluster...")
    order_clusters=sort2orderClusters(originalClusters,geneindex)
    
    # 根据阈值将 cluster 进行合并
    re_clusters=reCluster(order_clusters,rate)
    writeClusterList(re_clusters,outClusterFile)
    # print("done")
    

def imitateStart(geneindex:list[list[int]])-> list[int]:
    speciesNum=len(geneindex[0])
    start=[-9]*speciesNum
    for speciesIndex in range(speciesNum):
        a=[]
        for geneIndex in range(len(geneindex)):
            if geneindex[geneIndex][speciesIndex]!=-9:
                a.append(geneindex[geneIndex][speciesIndex])
        if a:
            start[speciesIndex] = min(a)-1
    return start

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

def geneindex_SampleName(geneindex:list[list[str]],sample_list: list[str],model="remove")-> list[list[str]]:
    models=["remove","add"]
    if model not in models:
        raise Exception(f"model must be in {models}")
    # 检查数据
    if len(geneindex[0])!=len(sample_list):
        raise Exception(f"geneindex columns {len(geneindex[0])} not equal to sample num {len(sample_list)}")
    
    result=[]
    if model=="remove":
        for line in geneindex:
            this_new_line=[]
            for sample_index,sample in enumerate(sample_list):
                genes=line[sample_index]
                if genes=="-":
                    this_new_line.append("-")
                elif "," in genes:
                    a=[]
                    genes_list=genes.split(",")
                    for gene in genes_list:
                        new_gene=gene[sample_name_len[sample_index]+1:]
                        a.append(new_gene)
                    this_new_line.append(",".join(a))
                else:
                    new_gene=genes[sample_name_len[sample_index]+1:]
                    this_new_line.append(new_gene)
            result.append(this_new_line)
    elif model=="add":
        for line in geneindex:
            this_new_line=[]
            for sample_index,sample in enumerate(sample_list):
                genes=line[sample_index]
                if genes=="-":
                    this_new_line.append("-")
                elif "," in genes:
                    a=[]
                    genes_list=genes.split(",")
                    for gene in genes_list:
                        new_gene=f"{sample}_{gene}"
                        a.append(new_gene)
                    this_new_line.append(",".join(a))
                else:
                    new_gene=f"{sample}_{genes}"
                    this_new_line.append(new_gene)
            result.append(this_new_line)
    return result
    

def sort2index(clusters: list[Cluster],geneindex: list[list[int]])->list[int]:
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
    
    return indexs

def removeImitateStartFromClusters(clusters:list[Cluster], start:list[int])-> bool:

    startIsIsolated=False
    isolatedIndex= -1

    removeImitateStart=False
    for index,cluster in enumerate(clusters):
        if cluster.start == start:
            cluster.content = cluster.content[1:]  # Remove the first line
            if cluster.content:  # Only add non-empty clusters
                cluster.start = cluster.updateStart()
                cluster.end = cluster.updateEnd()
                cluster.speciesNum = len(cluster.start)
                removeImitateStart = True
            elif not cluster.content:
                startIsIsolated = True
                isolatedIndex = index
        
    if startIsIsolated:
        clusters.pop(isolatedIndex)
        removeImitateStart = True
    
    return removeImitateStart

def run_cluster_python(infile:str,outfile:str)->None:
    """
    使用 python 进行 cluster
    """
    geneNumFile=infile
    outputFile=outfile

    start_Time = time.time()
    print("read gene index file...")
    geneindex=readNum(geneNumFile)

    print("Get geneindex start...")
    start=imitateStart(geneindex)
    
    print("Adding imitate start to geneindex...")
    geneindex.insert(0, start)

    # geneindex(List[List[int]]) 转为 List[Cluster]
    print("Clustering started...")
    lastClusterList=geneindex2clusters(geneindex)
    lastClusterNum=len(lastClusterList)
    
    #开始迭代
    print(f"initial clustering... {lastClusterNum} clusters (containing imitate start).")
    newClusterList=mergeClusters(lastClusterList)
    newClusterNum=len(newClusterList)

    iteration=0
    while(lastClusterNum!=newClusterNum):
        print(f"Iteration {iteration}: Cluster number {lastClusterNum} -> {newClusterNum}. Merging clusters ...")
        iteration += 1
        lastClusterNum=newClusterNum
        newClusterList=mergeClusters(newClusterList)
        newClusterNum=len(newClusterList)
    print(f"Clustering completed. {len(newClusterList)}")


    print("Removing imitate start from clusters...")
    isRemoved=removeImitateStartFromClusters(newClusterList, start)
    if isRemoved:
        print("Imitate start removed from clusters.")
    else:
        print("error:No imitate start found in clusters!!!")

    writeClusterList(newClusterList, outputFile)
    print(f"Clusters written to {outputFile}.")
    end_Time = time.time()
    duration = end_Time - start_Time
    print(f"Total time taken: {format_duration(duration)}")
    print("Done.")

def run_cluster_cpp(infile:str,outfile:str)->None:
    """
    使用 C++ 进行 cluster
    """


    BASE_DIR = Path(__file__).resolve().parent

    cluster_cpp_executable = BASE_DIR / "cluster"
    if not cluster_cpp_executable.exists():
        raise FileNotFoundError(f"C++ executable not found at {cluster_cpp_executable}")
    run_success = subprocess.run([cluster_cpp_executable, "-i", infile, "-o", outfile]).returncode
    if run_success != 0:
        raise RuntimeError("C++ clustering process failed.")
    else:
        print(f"C++ clustering completed. Clusters written to {outfile}.")

if __name__=="__main__":

    import argparse

    parser=argparse.ArgumentParser("Cluster gene index")
    parser.add_argument("-i","--input",type=str,help="gene index file",required=True)
    parser.add_argument("-o","--output",type=str,help="output file",required=True)
    args=parser.parse_args()

    geneNumFile=args.input
    outputFile=args.output
    run_cluster_cpp(geneNumFile,outputFile)
    
