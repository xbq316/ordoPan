

import math
from collections import defaultdict
import argparse

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

    def __str__(self):
        return f"content len:{len(self.content)}"

    def isLinked(self, other:"Cluster")-> bool:
        if not isinstance(other, Cluster):
            raise TypeError("Can only check linkage with another Cluster instance")
        for lastEnd,thisStart in zip(self.end, other.start):
            if lastEnd==-9 or thisStart==-9:
                continue
            if lastEnd+1==thisStart:
                return True #相邻
        return False #不相邻
    
    def suportedSpecies(self) -> int:
        """返回支持的物种数量"""
        return self.speciesNum-self.start.count(-9)

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

def read_cluster(infile:str)-> list[Cluster]:
    """
    读取 cluster 文件，返回 Cluster 对象列表
    """
    result_multiline = []
    with open(infile, 'r') as file:
        for line in file:
            line= line.strip()
            if not line:
                continue
            if line.startswith("#Cluster-"):
                result_multiline.append([])
                continue
            i=[int(x) if x !="-" else -9 for x in line.split("\t")]
            result_multiline[-1].append(i)
    result= []
    for lines in result_multiline:
        if lines:
            result.append(Cluster.frommultiline(lines))
    return result

def build_adjacency_lists(clusters:list[Cluster]) -> tuple[dict[int, list[tuple[int, float]]], dict[int, list[tuple[int, float]]]]:
    """构建邻接表"""
    pairs = []
    for i_index,i in enumerate(clusters):
        for j_index,j in enumerate(clusters[i_index+1:], start=i_index+1):
            if i_index==j_index:
                continue
            cutoff = i.cut_off(j)
            if cutoff > 0:
                pairs.append((i_index, j_index, cutoff))
    
    # 正向邻接表：key是起点，value是(终点, cutoff)列表
    forward_adj = defaultdict(list)
    # 反向邻接表：key是终点，value是(起点, cutoff)列表
    backward_adj = defaultdict(list)
    
    for start, end, cutoff in pairs:
        forward_adj[start].append((end, cutoff))
        backward_adj[end].append((start, cutoff))
    
    # 对邻接表按cutoff降序排序
    for node in forward_adj:
        forward_adj[node].sort(key=lambda x: x[1], reverse=True)
    for node in backward_adj:
        backward_adj[node].sort(key=lambda x: x[1], reverse=True)
    
    return forward_adj, backward_adj

def find_chains(forward_adj, backward_adj)->list[list[int]]:
    """寻找最长基因链"""
    # 所有节点的集合
    all_nodes = set(forward_adj.keys()).union(set(backward_adj.keys()))
    
    # 标记已访问的节点
    visited = set()
    chains = []
    
    # 优先处理没有入边的节点（链的起始点）
    start_candidates = [node for node in all_nodes if len(backward_adj.get(node, [])) == 0]
    
    # 然后处理剩余节点
    remaining_nodes = [node for node in all_nodes if node not in start_candidates]
    
    # 按节点连接数排序，优先处理连接数多的节点
    start_candidates.sort(key=lambda x: len(forward_adj.get(x, [])), reverse=True)
    remaining_nodes.sort(key=lambda x: len(forward_adj.get(x, [])), reverse=True)
    
    for node in start_candidates + remaining_nodes:
        if node in visited:
            continue
        
        # 尝试从这个节点开始构建链
        chain = [node]
        visited.add(node)
        
        # 向前扩展链
        current = node
        while True:
            # 找到最强的出边
            next_nodes = forward_adj.get(current, [])
            if not next_nodes:
                break
            
            # 找到第一个未访问的节点
            next_node = None
            for candidate, cutoff in next_nodes:
                if candidate not in visited:
                    next_node = candidate
                    break
            
            if next_node is None:
                break
            
            chain.append(next_node)
            visited.add(next_node)
            current = next_node
        
        # 向后扩展链（反向）
        current = node
        while True:
            # 找到最强的入边
            prev_nodes = backward_adj.get(current, [])
            if not prev_nodes:
                break
            
            # 找到第一个未访问的节点
            prev_node = None
            for candidate, cutoff in prev_nodes:
                if candidate not in visited:
                    prev_node = candidate
                    break
            
            if prev_node is None:
                break
            
            # 插入到链的前面
            chain.insert(0, prev_node)
            visited.add(prev_node)
            current = prev_node
        
        # 添加到链列表中
        if len(chain) > 0:  # 只保留长度大于0的链
            chains.append(chain)
    
    # 按链的长度排序
    chains.sort(key=lambda x: len(x), reverse=True)
    return chains

def write_chains_to_file(output_file:str,clusters:list[Cluster],chains:list[list[int]]):
    """将链写入文件"""
    with open(output_file, 'w') as f:
        for cluster_ID,chain in enumerate(chains):
            f.write(f"#Cluster-{cluster_ID}\n")
            for node in chain:
                for line in clusters[node].content:
                    f.write("\t".join(str(x) if x != -9 else "-" for x in line) + "\n")
    print(f"形成 {len(chains)} 个cluster,保存到 {output_file}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Cluster gene index based on cutoff values")
    parser.add_argument("-i","--input", type=str, required=True, help="Input cluster file")
    parser.add_argument("-o","--output", type=str, required=True, help="Output file for chains")
    args = parser.parse_args()

    infile= args.input
    outfile= args.output

    print("读取 cluster 数据")
    clusters = read_cluster(infile)
    print(f"读取到 {len(clusters)} 个 cluster")

    print("构建邻接表...")
    forward_adj, backward_adj = build_adjacency_lists(clusters)

    # 寻找链
    print("寻找基因链...")
    chains = find_chains(forward_adj, backward_adj)

    # 写入文件
    write_chains_to_file(outfile,clusters,chains)
