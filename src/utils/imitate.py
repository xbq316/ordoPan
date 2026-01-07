#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: imitate.py
# Created time: 2025/05/26 19:56:54
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python imitate.py -h


from __future__ import annotations
from itertools import combinations
from collections import defaultdict
import random
import numpy as np
import string
from typing import List, Tuple, Dict,Optional
import networkx as nx
import argparse
import os
from tqdm import tqdm
from scipy.cluster.hierarchy import linkage, dendrogram,to_tree
from scipy.stats import kendalltau
import editdistance
import matplotlib.pyplot as plt
import logging


# DEL        DUP        INS        INV        TRA 
# 0.31018622 0.06562710 0.28303792 0.04285394 0.29829482


# type  DEL     DUP     INS         INV         TRA
# mean  1.32    4.47    0.000485    49.1        5.43
# sd    3.55    13.7    0           154.        33.9 

# 5101.105
# 1368.689

# Insertion    Deletion     Inversion   Translocation
# 107151       76915        363         18455
# 0.52813923   0.37910826   0.0017892   0.09096331


# Insertion Deletion Inversion Duplication
# 125611    124744   6146      2364
# 0.485     0.4818   0.0237    0.00913

# 推导出sv频率
# Insertion   Deletion      Inversion     Duplication    Translocation
# 125611      124744        6146          2364           29931
# 0.43494716  0.43194504    0.02128146    0.00818571     0.10364063
# 由 一条染色体 5000 个基因为例的数据
EVENT_P={"DEL":0.43194504,"DUP":0.00818571,"INS":0.43494716,"INV":0.02128146,"TRA":0.10364063}
EVENT_MEAN={"DEL":1.32,"DUP":4.47,"INS":1.32,"INV":49.1,"TRA":5.43}
EVENT_STD={"DEL":3.55,"DUP":13.7,"INS":3.55,"INV":154,"TRA":33.9}

class Chromosome:
    def __init__(self,name:str|int, genes:list[str]):
        self.genes = genes
        self.name=name
        self.length=len(genes)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __str__(self):
        if len(self.genes)>10:
            return f"Chromosome(name={self.name}, genes={self.genes[:10]},length={self.length}...)"
        else:
            return f"Chromosome(name={self.name}, genes={self.genes},length={self.length})"
    
    def __len__(self):
        return self.length

    def __repr__(self):
        return self.__str__()
    
    def copy(self):
        return Chromosome(self.name, self.genes.copy())

    #倒位
    def inversion(self, start:int, end:int):
        start=start-1
        self.genes[start:end] = reversed(self.genes[start:end])
        #self.length=len(self.genes)

    #重复
    def duplication(self, start:int, end:int):
        start=start-1
        segment = self.genes[start:end]
        self.genes[end:end] = segment
        self.length=len(self.genes)

    #缺失
    def deletion(self, start:int, end:int):
        start=start-1
        del self.genes[start:end]
        self.length=len(self.genes)

    #插入
    def insertion(self, position:int, chromosome:Chromosome):
        self.genes[position:position]=chromosome.genes
        self.length=len(self.genes)

    #易位-染色体间
    def translocation(self, other_chromosome: Chromosome, start1: int, end1: int, start2: int, end2: int):
        """
        对两个染色体对象进行易位操作。若为同一染色体，则要求区间不重叠。
        所有索引为1-based，内部转换为0-based处理。
        """
        # 转为 0-based 索引
        start1 -= 1
        start2 -= 1

        # 获取片段副本
        seg1 = self.genes[start1:end1]
        seg2 = other_chromosome.genes[start2:end2]

        # 如果是同一条染色体，需要特别处理（避免切片赋值冲突）
        if self is other_chromosome:
            # 确保 start1 < start2（方便后续插入操作）
            if start1 > start2:
                start1, end1, seg1, start2, end2, seg2 = start2, end2, seg2, start1, end1, seg1

            # 删除靠后的段，插入交换内容
            del self.genes[start2:end2]
            self.genes[start2:start2] = seg1

            del self.genes[start1:end1]
            self.genes[start1:start1] = seg2

            self.length = len(self.genes)
        else:
            # 跨染色体，直接交换片段
            self.genes[start1:end1], other_chromosome.genes[start2:end2] = seg2, seg1
            self.length = len(self.genes)
            other_chromosome.length = len(other_chromosome.genes)

    @classmethod
    def generate(cls, name:str|None,geneNum:int)->Chromosome:
        genes=[]
        if name is None:
            name="".join(random.choices(string.ascii_uppercase + string.digits, k=5))
        for i in range(1,geneNum+1):
            genes.append(f"{name}g{i:07d}")
        return cls(name,genes)


    def write(self,path:str)->None:
        with open(path,"w") as f:
            for gene in self.genes:
                f.write(f"{gene}\n")

    def similarity_lcs(self,other_chromosome: Chromosome) -> float:
        seq1 = self.genes
        seq2 = other_chromosome.genes
        m, n = len(seq1), len(seq2)
        dp = [[0]*(n+1) for _ in range(m+1)]
        for i in range(m):
            for j in range(n):
                if seq1[i] == seq2[j]:
                    dp[i+1][j+1] = dp[i][j] + 1
                else:
                    dp[i+1][j+1] = max(dp[i+1][j], dp[i][j+1])
        lcs_len = dp[m][n]
        return (2 * lcs_len) / (len(seq1) + len(seq2)) if (seq1 or seq2) else 0.0

    def similarity_lcs_fast(self, other_chromosome: Chromosome) -> float:
        seq1 = self.genes
        seq2 = other_chromosome.genes
        dist = editdistance.eval(seq1, seq2)
        lcs_len = (len(seq1) + len(seq2) - dist) / 2
        return (2 * lcs_len) / (len(seq1) + len(seq2)) if (seq1 or seq2) else 0.0


    def similarity_kendall_tau(self,other_chromosome: Chromosome) -> float:
        """
        返回两个品系基因顺序的 Kendall tau 相似度（0~1，越大越相似）
        """
        seq1=self.genes
        seq2=other_chromosome.genes
        common_genes = set(seq1) & set(seq2)
        if len(common_genes) < 2:
            return 0.0  # 没法比，相似度最低
        s1 = [gene for gene in seq1 if gene in common_genes]
        s2 = [gene for gene in seq2 if gene in common_genes]
        gene_to_index = {gene: i for i, gene in enumerate(s1)}
        s2_mapped = [gene_to_index[gene] for gene in s2]
        tau, _ = kendalltau(range(len(s2_mapped)), s2_mapped)
        return tau if tau is not None else 0.0

    def similarity_breakpoint(self, other_chromosome) -> float:
        """
        值越大，越相似
        """
        seq1 = self.genes
        seq2 = other_chromosome.genes
        common_genes = set(seq1) & set(seq2)
        if len(common_genes) < 2:
            return 0.0
        s1 = [gene for gene in seq1 if gene in common_genes]
        s2 = [gene for gene in seq2 if gene in common_genes]
        s1 = ['START'] + s1 + ['END']
        s2 = ['START'] + s2 + ['END']
        pos2 = {gene: i for i, gene in enumerate(s2) if gene in s2}
        breakpoints = 0
        for i in range(len(s1) - 1):
            a, b = s1[i], s1[i + 1]
            if a not in pos2 or b not in pos2:
                breakpoints += 1  # treat missing as a breakpoint
            elif abs(pos2[a] - pos2[b]) != 1:
                breakpoints += 1
        max_breakpoints = len(common_genes) + 1
        return 1 - breakpoints / max_breakpoints
    
    def similarity_global_alignment(
        self, other_chromosome: Chromosome,
        match_score: int = 4, 
        mismatch_penalty: int = -1, 
        gap_penalty: int = -1
    ) -> float:
        """
        全局序列比对（Needleman-Wunsch算法），支持 List[str] 输入
        
        参数:
            seq1: 第一个序列（字符串列表）
            seq2: 第二个序列（字符串列表）
            match_score: 匹配得分（默认1）
            mismatch_penalty: 不匹配惩罚（默认-1）
            gap_penalty: 空位惩罚（默认-1）
        
        返回:
            比对得分和比对结果（均为 List[str] 格式）
        """
        seq1=self.genes
        seq2=other_chromosome.genes
        m, n = len(seq1), len(seq2)
        # 初始化得分矩阵
        score = [[0] * (n + 1) for _ in range(m + 1)]
        # 初始化回溯指针矩阵
        traceback = [[0] * (n + 1) for _ in range(m + 1)]
        
        # 初始化第一行和第一列
        for i in range(1, m + 1):
            score[i][0] = score[i-1][0] + gap_penalty
            traceback[i][0] = 2  # 向上
        for j in range(1, n + 1):
            score[0][j] = score[0][j-1] + gap_penalty
            traceback[0][j] = 3  # 向左
        
        # 填充得分矩阵
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = score[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
                delete = score[i-1][j] + gap_penalty
                insert = score[i][j-1] + gap_penalty
                
                score[i][j] = max(match, delete, insert)
                
                # 记录回溯方向
                if score[i][j] == match:
                    traceback[i][j] = 1  # 对角线
                elif score[i][j] == delete:
                    traceback[i][j] = 2  # 向上
                else:
                    traceback[i][j] = 3  # 向左
        
        # 回溯构建比对结果（使用append，最后反转）
        align1, align2 = [], []
        i, j = m, n
        
        while i > 0 or j > 0:
            if i > 0 and j > 0 and traceback[i][j] == 1:  # 对角线
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif i > 0 and traceback[i][j] == 2:  # 向上
                align1.append(seq1[i-1])
                align2.append("-")
                i -= 1
            else:  # 向左
                align1.append("-")
                align2.append(seq2[j-1])
                j -= 1
        # 反转结果
        align1.reverse()
        align2.reverse()
        matches = 0
        total = len(align1)
        for a, b in zip(align1, align2):
            if a == b:
                matches += 1
        return matches / total if total > 0 else 0.0

class Species:
    def __init__(self,name: str,chroList:List[Chromosome]):
        self.name=name
        self.chromosomes={}        
        for i in chroList:
            self.chromosomes[i.name]=i
    
    @classmethod
    def generate(cls,name:str,chroNum:int,geneNumMean:int=5101,geneNumStd:int=1368)->Species:
        """
        生成一个品系实例
        
        参数:
        name: 品系名称 (默认 "Bna")
        chroNum: 染色体数量 (默认 19)
        geneNumMean: 每个染色体基因数量的均值 (默认 5101)
        geneNumStd: 每个染色体基因数量的标准差 (默认 1368)
        """
        chroList=[]
        for i in range(1,chroNum+1):
            geneNum=int(random.gauss(geneNumMean,geneNumStd))
            if geneNum<1:
                geneNum=1
            chroList.append(Chromosome.generate(i,geneNum))
        return cls(name,chroList)

    @classmethod
    def ZS11(cls)->Species:
        chroIDs=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        geneNums=[4296,4203,5946,3014,5002,4533,3904,3222,7218,3021,5168,5574,8178,6183,5761,4597,5507,5349,6245]
        return Species("ZS11", [Chromosome.generate(f"{i}",j) for i,j in zip(chroIDs,geneNums)])

    def __str__(self):
        a=[]
        for chroname,genelist in self.chromosomes.items():
            if len(genelist.genes)<10:
                a.append(f"{chroname}:{','.join(genelist.genes)}")
            else:
                a.append(f"{chroname}:{','.join(genelist.genes[:10])}...")
        b="\n".join(a)
        return f"{str(self.name)}\n{b}"

    def __len__(self):
        return len(self.chromosomes)

    def __repr__(self):
        return self.__str__()

    def copy(self,name: str)->Species:
        return Species(name,[i.copy() for i in self.chromosomes.values()])

    def chromosomeNum(self)->int:
        return len(self.chromosomes)


    def random_subinterval(self,length:int)-> Tuple[str, int, int]:
        """
        返回一个随机的子区间
        分别是染色体编号,开始位置,结束位置
        """
        chromosome = random.choice(list(self.chromosomes.keys()))
        end=self.chromosomes[chromosome].length
        if end-length<=0:
            raise ValueError(f"变异过程中,染色体长度已缩减至0,请尝试改变(变异次数,基因数量,随机种子)")
        sub_start = random.randint(1, end-length)
        return (chromosome,sub_start,sub_start+length)
 
    def random_position(self)-> Tuple[str, int]:
        """
        返回一个随机的位置
        """
        chromosome = random.choice(list(self.chromosomes.keys()))
        position = random.randint(1, self.chromosomes[chromosome].length)
        return (chromosome,position)

    def random_event(self)-> str:
        """
        返回一个随机的事件和长度
        """
        event = random.choices(list(EVENT_P.keys()), weights=list(EVENT_P.values()), k=1)[0]
        return event

    def random_length(self,event:str)->int|Tuple[int, int]:
        """
        返回一个随机的事件的长度
        根据基因组大小，更新变异长度
        EVENT_MEAN
        """
        if event == "TRA":
            # 小变异不管
            # 大变异保留
            length1 = int(random.gauss(EVENT_MEAN[event], EVENT_STD[event])*GENOME_ZOOM)
            length1=max(1,length1)
            length1=int(min((EVENT_STD[event]*3+EVENT_MEAN[event])*GENOME_ZOOM,length1))
            
            length2 = int(random.gauss(EVENT_MEAN[event], EVENT_STD[event])*GENOME_ZOOM)
            length2=max(1,length2)
            length2=int(min((EVENT_STD[event]*3+EVENT_MEAN[event])*GENOME_ZOOM,length2))
            return (length1,length2)
        else:
            # 指数分布
            length=int(np.random.exponential(scale=1/EVENT_MEAN[event], size=1)[0]*GENOME_ZOOM)
            length=max(1,length)
            length=int(min((EVENT_STD[event]*3+EVENT_MEAN[event])*GENOME_ZOOM,length))
            return length

    def new(self,name:str,n:int)->Species:
        newSpecies=self.copy(name)
        events={}
        n=max(1,int(n))
        for i in range(n):
            event=newSpecies.random_event()
            if event=="DEL":
                length=newSpecies.random_length(event)
                chro,start,end=newSpecies.random_subinterval(length)
                newSpecies.chromosomes[chro].deletion(start,end)
            elif event=="DUP":
                length=newSpecies.random_length(event)
                chro,start,end=newSpecies.random_subinterval(length)
                newSpecies.chromosomes[chro].duplication(start,end)
            elif event=="INS":
                length=newSpecies.random_length(event)
                chro,position=newSpecies.random_position()
                random_chromosome=Chromosome.generate(None,length)
                newSpecies.chromosomes[chro].insertion(position,random_chromosome)
            elif event=="INV":
                length=newSpecies.random_length(event)
                chro,start,end=newSpecies.random_subinterval(length)
                newSpecies.chromosomes[chro].inversion(start,end)
            elif event=="TRA":
                max_attempts = 100
                for attempt in range(max_attempts):
                    length1, length2 = newSpecies.random_length(event)
                    chro1, start1, end1 = newSpecies.random_subinterval(length1)
                    chro2, start2, end2 = newSpecies.random_subinterval(length2)
                    if chro1 != chro2:
                        break  # 染色体不同，不可能重叠，跳出循环
                    # 同一染色体，判断是否不重叠
                    if end1 <= start2 or end2 <= start1:
                        break  # 区间不重叠，合法，跳出循环
                else:
                    raise RuntimeError("Failed to find non-overlapping intervals after 100 attempts")
                # 最终执行易位操作
                newSpecies.chromosomes[chro1].translocation(newSpecies.chromosomes[chro2], start1, end1, start2, end2)
            else:
                pass
            events[event]=events.get(event,0)+1
        logging.info(f"品系 {self.name} 变异 {n} 次，事件分布为:{str(events)} 形成新品系 {newSpecies.name}")
        return newSpecies
    
    def dotplot(self, other: Optional['Species'] = None):
        if other is None:
            other = self

        # 获取基因线性位置信息，每个染色体独立编号
        def get_gene_positions(species):
            positions = {}
            chro_offsets = {}
            chro_names = sorted(species.chromosomes.keys())  # 染色体按名字排序
            offset = 0
            for idx, chro_name in enumerate(chro_names):
                chro = species.chromosomes[chro_name]
                chro_offsets[chro_name] = offset
                for i, gene in enumerate(chro.genes):
                    positions[gene] = (idx, i)  # 染色体编号+染色体内位置
                offset += chro.length
            return positions, chro_names

        pos1, chro_order1 = get_gene_positions(self)
        pos2, chro_order2 = get_gene_positions(other)

        shared_genes = set(pos1.keys()) & set(pos2.keys())
        x = []
        y = []
        for gene in shared_genes:
            chr_idx1, gene_pos1 = pos1[gene]
            chr_idx2, gene_pos2 = pos2[gene]
            x.append((chr_idx1, gene_pos1))
            y.append((chr_idx2, gene_pos2))

        # 映射为全局位置：每个染色体独立编号，位置作为染色体内的位置
        def convert_to_axis(coords, chro_sizes):
            axis_coords = []
            chro_boundaries = []
            for chr_idx, pos in coords:
                offset = sum(chro_sizes[i] for i in range(chr_idx))
                axis_coords.append(offset + pos)
            # 计算染色体边界位置
            boundary = 0
            for size in chro_sizes:
                boundary += size
                chro_boundaries.append(boundary)
            return axis_coords, chro_boundaries

        chro_sizes1 = [self.chromosomes[chro].length for chro in chro_order1]
        chro_sizes2 = [other.chromosomes[chro].length for chro in chro_order2]


        x_vals, x_lines = convert_to_axis(x, chro_sizes1)
        y_vals, y_lines = convert_to_axis(y, chro_sizes2)

        # 绘图
        plt.figure(figsize=(12, 10))
        plt.scatter(x_vals, y_vals, s=0.05,color='blue',)
        plt.title(f"Dotplot: {self.name} vs {other.name}")
        plt.xlabel(f"{self.name} chromosomes")
        plt.ylabel(f"{other.name} chromosomes")
        plt.xlim(0, sum(chro_sizes1))
        plt.ylim(0, sum(chro_sizes2))

        # 添加染色体分界线
        for xl in x_lines[:-1]:
            plt.axvline(x=xl, color='gray', linestyle='--', linewidth=0.5)
        for yl in y_lines[:-1]:
            plt.axhline(y=yl, color='gray', linestyle='--', linewidth=0.5)

        # 设置坐标刻度（显示染色体编号）
        def label_positions(chro_sizes, chro_order):
            ticks = []
            labels = []
            acc = 0
            for size, name in zip(chro_sizes, chro_order):
                ticks.append(acc + size // 2)
                labels.append(str(name))
                acc += size
            return ticks, labels

        xticks, xlabels = label_positions(chro_sizes1, chro_order1)
        yticks, ylabels = label_positions(chro_sizes2, chro_order2)

        plt.xticks(xticks, xlabels, rotation=90)
        plt.yticks(yticks, ylabels)
        plt.grid(False)
        plt.tight_layout()
        plt.show()



    
class Population:
    def __init__(self,ancestor:Species,G:nx.DiGraph=None) -> None:
        self.ancestor=ancestor
        self.species={ancestor.name:ancestor}
        
        if G is not None:
            self.G=G
            if not self.is_directed_tree():
                raise ValueError("G 不是有向树")
            # 判断祖先节点是否为根节点
            root = next(nx.topological_sort(G))
            if root != self.ancestor.name:
                raise ValueError(f"祖先节点 {self.ancestor.name} 不是根节点")
        else:
            self.G=nx.DiGraph()
            self.G.add_node(self.ancestor.name)
    def __str__(self):
        
        leaf_species = self.leaf_species()
        leaf_names = [species.name for species in leaf_species]
        # 展示一定数量的叶子品系名称
        if len(leaf_names)<5:
            return f"Population(ancestor: {self.ancestor.name}, species count: {len(self.species)},leaf species: {', '.join(leaf_names)})"
        else:
            return f"Population(ancestor: {self.ancestor.name}, species count: {len(self.species)},leaf species: {', '.join(leaf_names[:5])}...)"
    
    def __repr__(self):
        return self.__str__()
    
    def is_directed_tree(self):
        # 1. 必须是 nx.DiGraph()
        if not isinstance(self.G, nx.DiGraph):
            return False
        # 2. 检查是否弱连通（忽略方向后是否连通）
        if not nx.is_weakly_connected(self.G):
            return False
        # 3. 检查是否有环
        if not nx.is_directed_acyclic_graph(self.G):
            return False
        # 4. 检查根节点数量（入度为0的节点）
        roots = [n for n in self.G.nodes() if self.G.in_degree(n) == 0]
        if len(roots) != 1:
            return False
        # 5. 检查其他节点入度是否为1
        for node in self.G.nodes():
            if node != roots[0] and self.G.in_degree(node) != 1:
                return False
        return True

    def run(self) -> None:
        """
        从有向无环图中生成品系。
        DAG 的边可以包含 "weight" 属性，表示变异次数。
        节点必须包含起始品系 self.ancestor.name。
        """
        if not self.is_directed_tree():
            raise ValueError("G 不是有向树")

        # 判断祖先节点是否为根节点
        root = next(nx.topological_sort(self.G))
        if root != self.ancestor.name:
            raise ValueError(f"祖先节点 {self.ancestor.name} 不是根节点")

        # 按拓扑排序顺序遍历节点，跳过祖先节点（因为祖先已存在）
        topo_order = list(nx.topological_sort(self.G))
        for node in topo_order:
            if node == self.ancestor.name:
                continue
            # 找到所有指向当前节点的前驱节点（可能多个）
            preds = list(self.G.predecessors(node))
            if not preds:
                # 如果没有前驱节点，说明是孤立点，跳过或报错
                raise ValueError(f"节点 {node} 没有前驱节点，不能生成品系")
            
            # 这里简单处理，选第一个前驱节点作为祖先
            parent_name = preds[0]
            parent_species = self.species.get(parent_name)
            if parent_species is None:
                raise ValueError(f"前驱品系 {parent_name} 尚未生成")
            
            # 读取边的变异次数weight，默认为1
            edge_data = self.G.get_edge_data(parent_name, node)
            n = edge_data.get("weight", 1)
            
            # 基于父品系生成新品系
            new_species = parent_species.new(node, n)
            self.species[node] = new_species
        
        logging.info(f"演化物种个数: {len(self.species)}")
        logging.info(f"现有物种(叶子节点)个数: {len(self.leaf_species())}")

    def get_species(self, name: str) -> Species | None:
        """返回指定名称的品系实例，找不到返回None"""
        return self.species.get(name)

    def all_species(self) -> List[Species]:
        """返回所有品系列表"""
        return list(self.species.values())

    def chromosome(self) -> Dict[str,List[Chromosome]]:
        """
        返回所有染色体列表
        """
        result=defaultdict(list)
        for species in self.species.values():
            for chro in species.chromosomes.values():
                result[chro.name].append(chro)
        return result

    def to_newick(self) -> str:
        def build_newick(node):
            children = list(self.G.successors(node))
            if not children:
                return node
            else:
                subtrees = []
                for child in children:
                    weight = self.G[node][child].get("weight", 1)  # 默认权重为1
                    subtrees.append(build_newick(child) + f":{weight}")
                return "(" + ",".join(subtrees) + ")" + node

        root = [n for n in self.G.nodes if self.G.in_degree(n) == 0][0]
        return build_newick(root) + ";"


    def leaf_species(self) -> List[Species]:
        """返回所有叶子节点的品系实例列表"""
        leaves = [n for n in self.G.nodes if self.G.out_degree(n) == 0]
        return [self.species[n] for n in leaves if n in self.species]


    def generate_binary_tree(self, N=100, weight=1):
        """
        生成一个正好有N个叶子的二叉树（有向树），并更新self.G。
        根节点名称使用 self.ancestor.name。
        """
        G = nx.DiGraph()
        leaves = [f"spe{i}" for i in range(1, N + 1)]
        for leaf in leaves:
            G.add_node(leaf)

        current_node_id = N + 1
        nodes = leaves.copy()

        while len(nodes) > 1:
            # 顺序合并前两个节点
            a = nodes.pop(0)
            b = nodes.pop(0)
            parent = f"spe{current_node_id}"
            current_node_id += 1
            G.add_node(parent)
            G.add_edge(parent, a, weight=weight)
            G.add_edge(parent, b, weight=weight)
            nodes.append(parent)

        old_root = nodes[0]  # 最终根节点的旧名字
        new_root = self.ancestor.name

        if old_root != new_root:
            # 重命名根节点
            G = nx.relabel_nodes(G, {old_root: new_root})

        self.G = G
        self.ancestor.name = new_root

        if not self.is_directed_tree():
            raise ValueError("生成的图不是有向树")

    def reset_weight(self,weight:int)->None:
        for u, v in self.G.edges:
            self.G[u][v]['weight'] = weight

    def geneindex(self,outdir,chroList:List[int])->None:
        SpeciesList=self.leaf_species()
        speList=[i.name for i in SpeciesList]
        if chroList is None:
            chroList=list(self.ancestor.chromosomes.keys())
               
        for chro in chroList:
            thischroname=chro
            chromosomeList=[ self.species[species_name].chromosomes[chro] for species_name in speList ]
            G=nx.Graph()
            indexs=list(range(len(chromosomeList)))

            geneID_index2num={}
            for index in indexs:
                chro=chromosomeList[index]
                for geneID in chro.genes:
                    geneID_index2num[(geneID,index)]=geneID_index2num.get((geneID,index),0)+1

            thiscombinations=list(combinations(indexs, 2))
            for c in tqdm(thiscombinations,desc="构建图..."):  # 选取2个元素的组合
                index1,index2=c
                chro1=chromosomeList[index1]
                chro2=chromosomeList[index2]
                for node1 in chro1.genes:
                    G.add_node((node1,index1))
                for node2 in chro2.genes:
                    G.add_node((node2,index2))
                genesset1=set(chro1.genes)
                genesset2=set(chro2.genes)
                Intersection=genesset1&genesset2
                for node in Intersection:
                    G.add_edge((node,index1),(node,index2))# 模拟序列比对，同源基因的发现

            gene2index={}
            for chro,chromosome in zip(indexs,chromosomeList):
                for index,gene in enumerate(chromosome.genes):
                    gene2index[(chro,gene)]=index+1
       
            visited=set()
            geneindex_txt=os.path.join(outdir,f"{thischroname}.txt")
            logging.info(f"生成 geneindex txt: {geneindex_txt}")
            geneindex_txt=open(geneindex_txt,"w")# 写入祖先基因的 geneindex 文件

            geneindex_num=os.path.join(outdir,f"{thischroname}.num")
            logging.info(f"生成 geneindex num : {geneindex_num}")
            geneindex_num=open(geneindex_num,"w")
            
            allnodes=list(G.nodes)
            for node in tqdm(allnodes,desc="生成 geneindex ..."):
                if node in visited:continue
                # 深度优先搜索
                path=list(nx.dfs_tree(G, node))
                for i in path:
                    visited.add(i)
                path.sort(key=lambda x: x[1])
                a=["" for i in indexs]
                for i in path:
                    a[i[1]]=i[0]
                # 写入 {chro}.txt
                b=[]
                for index,gene in enumerate(a):
                    if not gene :
                        b.append("-")
                    else:
                        geneNum=geneID_index2num[(gene,index)]
                        c=[gene]*geneNum
                        b.append(",".join(c))
                geneindex_txt.write("\t".join(b)+"\n")
                
                # 写入 {chro}.num
                axsd=[]
                for index,genes in enumerate(b):
                    if genes=="-":axsd.append("-")
                    elif "," in genes:
                        geneList=genes.split(",")
                        gene=geneList[0]
                        axsd.append(gene2index[index,gene])
                    else:
                        axsd.append(gene2index[(index,genes)])
                geneindex_num.write("\t".join(map(str,axsd))+"\n")
            geneindex_txt.close()
            geneindex_num.close()
            # 写入祖先基因顺序 list
            ancestor_gene_list=os.path.join(outdir,f"{thischroname}.{self.ancestor.name}.genelist")
            self.ancestor.chromosomes[thischroname].write(ancestor_gene_list)
            logging.info(f"生成祖先基因列表 {ancestor_gene_list}")
            print(f"完成染色体 {thischroname}")

    def orthogroup(self,outdir):
        logging.info("*"*50)
        logging.info("orthogroup...")
        SpeciesList = self.leaf_species()
        speList = [i.name for i in SpeciesList]
        chroList = list(self.ancestor.chromosomes.keys())

        # 生成 map file
        mapfile=os.path.join(outdir,"map.txt")
        logging.info(f"map file(newgene -> oldgene) {mapfile}")

        # Step 1: 生成 GFF 文件并重命名基因
        oldgene2newgene = {}
        for species in tqdm(SpeciesList, desc="生成GFF文件..."):
            species_name = species.name
            oldgene2newgene[species_name] = {}
            thisgff=os.path.join(outdir,f"{species_name}.gff")
            logging.info(f"生成GFF文件 {thisgff}")
            with open(thisgff, "w") as fw , open(mapfile,"a") as fm:
                for chroname in chroList:
                    for index, gene in enumerate(species.chromosomes[chroname].genes):
                        new_genename = f"{species_name}_{chroname}g{index+1:07d}"
                        fm.write(f"{new_genename}\t{gene}\n")
                        if gene not in oldgene2newgene[species_name]:
                            oldgene2newgene[species_name][gene] = [new_genename]
                        else:
                            oldgene2newgene[species_name][gene].append(new_genename)
                        start = (index + 1)*1500
                        end = (index + 1 + 1)*1500 - 1
                        fw.write(f"{species_name}_{chroname}\t{new_genename}\t{start}\t{end}\n")

        # Step 2: 每个品系所有基因（不分染色体）
        species_allgenes = {}
        for species in SpeciesList:
            all_genes = []
            for chro in species.chromosomes.values():
                all_genes.extend(chro.genes)
            species_allgenes[species.name] = all_genes

        # Step 3: 构建同源图
        G = nx.Graph()
        for idx, sp_name in enumerate(speList):
            for gene in species_allgenes[sp_name]:
                G.add_node((gene, idx))

        for i, j in tqdm(list(combinations(range(SPECIES_NUM), 2)), desc="构建同源基因图..."):
            genes_i = set(species_allgenes[speList[i]])
            genes_j = set(species_allgenes[speList[j]])
            shared = genes_i & genes_j
            for g in shared:
                G.add_edge((g, i), (g, j))

        # Step 4: 构建 Orthogroups.tsv
        orthogroup_file = os.path.join(outdir, "Orthogroups.tsv")
        with open(orthogroup_file, "w") as fw:
            fw.write("Orthogroup\t" + "\t".join(speList) + "\n")
            visited = set()
            og_index = 0
            for node in tqdm(G.nodes, desc="生成 Orthogroups.tsv ..."):
                if node in visited:
                    continue
                group_nodes = list(nx.dfs_tree(G, node))
                for n in group_nodes:
                    visited.add(n)

                line = [[] for _ in range(SPECIES_NUM)]
                for gene, sp_idx in group_nodes:
                    renamed_genes = oldgene2newgene[speList[sp_idx]][gene]
                    line[sp_idx].extend(renamed_genes)

                # 至少两个品系非空才保留
                nonempty = sum(1 for genes in line if genes)
                if nonempty < 2:
                    continue

                fw.write(f"OG{og_index:07d}")
                for genes in line:
                    fw.write("\t" + ", ".join(genes) if genes else "\t")
                fw.write("\n")
                og_index += 1
        print("Orthogroups.tsv 构建完成（使用新基因名，品系内同源已处理）")
        logging.info(f"Orthogroups.tsv : {orthogroup_file}")
        logging.info("*"*50)

    def leaf_distance(self,leaf1,leaf2):
        # 找到两者的最近公共祖先（LCA）
        ancestors1 = nx.ancestors(self.G, leaf1) | {leaf1} #获得这个叶子节点全部祖先节点 -> set() 包括本节点
        ancestors2 = nx.ancestors(self.G, leaf2) | {leaf2}

        # 找到所有共有祖先
        common_ancestors = ancestors1 & ancestors2
        # 计算从 原始祖先 到 每个共有祖先 的距离
        # 距离越长，则表示 共有祖先 距离两个叶子节点越近
        lca = max(common_ancestors, key=lambda x: nx.shortest_path_length(self.G, source=self.ancestor.name, target=x))

        # 计算从LCA到两个叶子的路径长度
        path1 = nx.shortest_path(self.G, source=lca, target=leaf1)
        path2 = nx.shortest_path(self.G, source=lca, target=leaf2)

        def path_length(path):
            return sum(self.G[path[i]][path[i+1]]['weight'] for i in range(len(path)-1))

        total_distance = path_length(path1) + path_length(path2)

        return total_distance

    def leaf_distance_mean(self,iterations=100):
        # 找出所有叶子节点（没有子节点）
        leaves = [n for n in self.G.nodes if self.G.out_degree(n) == 0]
        if len(leaves) < 2:
            raise ValueError("树中叶子节点不足两个")
        
        distance_list=[]
        for i in range(iterations):
            # 随机抽两个不同的叶子节点
            leaf1, leaf2 = random.sample(leaves, 2)
            d=self.leaf_distance(leaf1,leaf2)
            distance_list.append(d)
        return sum(distance_list)/iterations
    
    def togenelist(self,outputdir:str,chromosome_name:str):
        if chromosome_name not in self.ancestor.chromosomes:
            raise ValueError(f"不存在名为 {chromosome_name} 的染色体")
        SpeciesList = self.leaf_species()
        for species in SpeciesList:
            with open(os.path.join(outputdir,f"{chromosome_name}.{species.name}.genelist.txt"), "w") as f:
                for gene in species.chromosomes[chromosome_name].genes:
                    f.write(f"{gene}\n")

    def predict_tree(self,outputdir:str,chro:int,similarity_method:str="similarity_breakpoint",linkage_method: str = 'ward',show_all:bool=False):
        """
        similarity_method 支持的方法:similarity_lcs/similarity_kendall_tau/similarity_breakpoint/similarity_global_alignment/similarity_lcs_fast
        linkage_method 支持的方法:single/complete/average/weighted/centroid/median/ward
        """
        similarity_methods=["similarity_lcs","similarity_kendall_tau","similarity_breakpoint","similarity_global_alignment","similarity_lcs_fast"]
        valid_linkage_methods = ["single", "complete", "average", "weighted", "centroid", "median", "ward"]

        species_List = self.leaf_species()

        def build_distance_matrix(data: List[Species], chro: int, similarity_method: str) -> np.ndarray:
            n = len(data)
            distance_matrix = np.zeros((n, n))
            for i in tqdm(range(n), desc=f"{similarity_method} 距离矩阵计算中..."):
                chr_i = data[i].chromosomes[chro]
                func_i = getattr(chr_i, similarity_method)
                for j in range(i + 1, n):
                    chr_j = data[j].chromosomes[chro]
                    score = func_i(chr_j)
                    distance_matrix[i, j] = 1.0 - score
                    distance_matrix[j, i] = distance_matrix[i, j]
            return distance_matrix
        
        def tree_to_newick(linked: np.ndarray, labels: List[str]) -> str:
            """
            将 linkage 矩阵和标签转换为 Newick 字符串
            
            参数:
                linked: linkage 矩阵
                labels: 品系名称列表，顺序需与距离矩阵一致
            
            返回:
                Newick 格式的字符串 ,如 "(A:0.1,B:0.2):0.3;"
            """
            tree = to_tree(linked, rd=False)
            def _build_newick(node, parent_dist):
                if node.is_leaf():
                    return f"{labels[node.id]}:{parent_dist - node.dist:.6f}"
                else:
                    left = _build_newick(node.left, node.dist)
                    right = _build_newick(node.right, node.dist)
                    return f"({left},{right}):{parent_dist - node.dist:.6f}"
            return _build_newick(tree, tree.dist) + ";"
        
        if show_all:
            for similarity_method in similarity_methods:
                distance_matrix = build_distance_matrix(data=species_List,chro=chro,similarity_method=similarity_method)
                for linkage_method in valid_linkage_methods:
                    linked = linkage(distance_matrix[np.triu_indices(len(distance_matrix), k=1)], method=linkage_method)
                    lables=[species.name for species in species_List]
                    thisnwk=tree_to_newick(linked,lables)
                    with open(os.path.join(outputdir,f"{chro}.{similarity_method}.{linkage_method}.nwk"), "w") as f:
                        f.write(thisnwk+"\n")
        else:
            if similarity_method not in similarity_methods:
                raise ValueError(f"无效的 similarity_method: {similarity_method}. 有效方法为 {' '.join(similarity_methods)}")
            if linkage_method not in valid_linkage_methods:
                raise ValueError(f"无效的 linkage_method: {linkage_method}. 有效方法为 {' '.join(valid_linkage_methods)}")
            distance_matrix = build_distance_matrix(data=species_List,chro=chro,similarity_method=similarity_method)
            linked = linkage(distance_matrix[np.triu_indices(len(distance_matrix), k=1)], method=linkage_method)
            lables=[species.name for species in species_List]

            thisnwk=tree_to_newick(linked,lables)
            with open(os.path.join(outputdir,f"{chro}.{similarity_method}.{linkage_method}.nwk"), "w") as f:
                f.write(thisnwk+"\n")


if __name__ == "__main__":
    
    helptext=rf"""
    Usage: python imitate.py -sn 20 -vr 0.08 -s 42 -gn 5000 -c chr1 chr2 chr3
    EVENT_P={EVENT_P}
    EVENT_MEAN={EVENT_MEAN}
    EVENT_STD={EVENT_STD}
    """
    
    parser=argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-sn","--species_num",type=int,default=20,help="species number")
    #parser.add_argument("-vn","--variation_num",type=int,default=400,help="每个品系发生变异的次数")
    parser.add_argument("-vr","--variation_rate",type=float,default=0.08,help="每个品系发生变异的比例")
    parser.add_argument("--leaf",action='store_false',help="变异次数是发生在 叶子节点之间还是演化品系之间 , 默认为叶子节点之间")
    parser.add_argument("-s","--seed",type=int,default=42,help="随机种子")
    #parser.add_argument("-cn","--chromosome_num",type=int,default=1,help="染色体数量")
    parser.add_argument("-gn","--gene_num",type=int,default=1500,help="染色体基因数量")
    parser.add_argument("-o","--output_dir",type=str,required=True,help="输出目录")
    parser.add_argument("--Orthogroup",action="store_true",help="是否输出 Orthogroup 以及对应的GFF(MCScanX格式)文件(包含所有染色体)")
    
    parser.add_argument("-c","--chromosomes",type=int,nargs="+",default=[1],help="需要展示的染色体编号,默认展示1号染色体,可以指定多个,如 --chromosomes 1 2 3")

    # 生成 genelist 文件
    parser.add_argument("-gl","--genelist",action="store_true",help="是否输出基因列表,默认输出祖先基因列表")
    parser.add_argument("-agl","--ancestor_genelist",action="store_false",help="是否输出祖先基因列表,默认输出基因列表")

    # 生成 geneindex 文件
    parser.add_argument("--geneindex",action="store_false",help="生成 geneindex 文件,默认生成")

    # 正确进化树
    parser.add_argument("--true_tree",action="store_false",help="是否输出真实树,默认输出")

    # 预测树相关参数
    parser.add_argument("-pt","--predict_tree",action="store_true",help="是否预测树")
    parser.add_argument("--linkage_method",type=str,default="ward",
                        choices=["single","complete","average","weighted","centroid","median","ward"],
                        help="若是预测树，选择聚类方法")
    parser.add_argument("--similarity_method",type=str,default="similarity_breakpoint",
                        choices=["similarity_lcs","similarity_kendall_tau","similarity_breakpoint","similarity_global_alignment","similarity_lcs_fast"],
                        help="若是预测树，选择品系相似度计算方法")
    parser.add_argument("--show_all_methods",action="store_true",help="若是预测树，是否输出所有方法的品系树")

    args=parser.parse_args()
    SEED=args.seed
    GENE_NUM=args.gene_num
    VARIATION_NUM=max(1,int(args.variation_rate*GENE_NUM))
    SPECIES_NUM=args.species_num
    CHROMOSOME_NUM=len(args.chromosomes)
    OUTPUT_DIR=args.output_dir
    LEAF_FLAG=args.leaf
    SIMILARITY_METHOD=args.similarity_method
    LINKAGE_METHOD=args.linkage_method
    SHOW_TREE=args.predict_tree
    SHOW_ALL_METHODS_TREE=args.show_all_methods
    SHOW_TRUE_TREE=args.true_tree
    SHOW_GENEINDEX=args.geneindex
    SHOW_CHROMOSOMES_LSIT=args.chromosomes
    SHOW_GENELIST=args.genelist
    SHOW_ORTHOGROUP=args.Orthogroup
    GENOME_ZOOM=GENE_NUM/5000
    random.seed(SEED)
    np.random.seed(SEED)

    ancestor=Species.generate("ancestor",CHROMOSOME_NUM,GENE_NUM,GENE_NUM/5)

    mypopulation=Population(ancestor)
    mypopulation.generate_binary_tree(SPECIES_NUM)
    
    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    logfile=os.path.join(OUTPUT_DIR,'variation.log')
    with open(logfile, 'w'):
        pass  # 这会清空文件内容

    # 配置日志输出到文件
    logging.basicConfig(
        filename=logfile,    # 日志文件名
        level=logging.INFO,          # 日志级别
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logging.info(f"脚本名称 {__file__}")
    logging.info("开始演化")

    # 设置变异次数
    if LEAF_FLAG:
        #叶子节点之间变异次数
        b=mypopulation.leaf_distance_mean()
        weight=VARIATION_NUM/b
        mypopulation.reset_weight(weight=weight)
        logging.info(f"叶子节点之间变异次数为 {VARIATION_NUM},演化品系之间变异次数为 {weight}")
    else:
        #演化品系之间变异次数
        mypopulation.reset_weight(VARIATION_NUM)
        logging.info(f"演化品系之间变异次数为 {VARIATION_NUM}")

    # 根据生成的树，来演化种群
    mypopulation.run()

    
    if SHOW_TRUE_TREE:
        nwk=mypopulation.to_newick()
        with open(os.path.join(OUTPUT_DIR,f"true_tree.nwk"),"w") as f:
            f.write(nwk)
        logging.info(f"输出真实树: {os.path.join(OUTPUT_DIR,f'true_tree.nwk')}")

    if SHOW_GENEINDEX:
        mypopulation.geneindex(f"{OUTPUT_DIR}",SHOW_CHROMOSOMES_LSIT)

    if SHOW_GENELIST:
        for chroname in SHOW_CHROMOSOMES_LSIT:
            mypopulation.togenelist(f"{OUTPUT_DIR}",chroname)

    if SHOW_ORTHOGROUP:
        mypopulation.orthogroup(f"{OUTPUT_DIR}")
    
    #similarity_lcs","similarity_kendall_tau","similarity_breakpoint,"similarity_global_alignment"
    if SHOW_TREE:
        if SHOW_ALL_METHODS_TREE:
            for chroname in SHOW_CHROMOSOMES_LSIT:
                mypopulation.predict_tree(f"{OUTPUT_DIR}",chro=chroname,show_all=True)
        else:
            for chroname in SHOW_CHROMOSOMES_LSIT:
                mypopulation.predict_tree(f"{OUTPUT_DIR}",chro=chroname,similarity_method=SIMILARITY_METHOD,linkage_method=LINKAGE_METHOD,show_all=False)
    logging.shutdown()
