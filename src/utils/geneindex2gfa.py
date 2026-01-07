#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: geneindex2gfa.py
# Created time: 2025/11/10 22:46:30
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python geneindex2gfa.py -h


from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Iterable
from collections import Counter

# ============================================================
# ===============  内置 cluster 模块的必要部分  ===============
# ============================================================

@dataclass
class Gene:
    chromosome: str
    start: int
    end: int
    name: str
    score: str
    strand: str

    def __str__(self):
        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.name}\t{self.score}\t{self.strand}"


def read_sample(file: str) -> List[str]:
    """读取样本名列表"""
    with open(file) as fh:
        return [line.strip() for line in fh if line.strip()]


def read_geneindex(file: str) -> List[List[str]]:
    """读取 gene index 文件（行代表节点，列代表样本）"""
    res: List[List[str]] = []
    with open(file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            row = [cell if cell else "-" for cell in line.split("\t")]
            res.append(row)
    return res


def read_bed6(file: str) -> List[Gene]:
    """读取 BED6 格式"""
    res: List[Gene] = []
    with open(file) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) == 4:
                chrom, start, end, name = parts[:4]
                res.append(Gene(chrom, int(start), int(end), name, ".", "+"))
            elif len(parts) == 6:
                chrom, start, end, name, score, strand = parts[:6]
                res.append(Gene(chrom, int(start), int(end), name, score, strand))
    return res


# ============================================================
# ================  主逻辑：pan_to_gfa 功能  =================
# ============================================================

@dataclass(frozen=True)
class Edge:
    name1: int | str
    strand1: str
    name2: int | str
    strand2: str
    overlap: int = 0

    def __str__(self) -> str:
    # 将方向与节点名拼接为一个 token（如 "1+" 或 "SG0000000-"）
        return f"L\t{self.name1}\t{self.strand1}\t{self.name2}\t{self.strand2}\t{self.overlap}M"

    def core(self) -> Tuple:
        return (self.name1, self.strand1, self.name2, self.strand2)


class Edges:
    """保持唯一 Edge 的集合"""
    __slots__ = ("_edges", "_set")
    def __init__(self, edges: Iterable[Edge] = ()):
        self._edges: List[Edge] = []
        self._set: set = set()
        for e in edges:
            self.add(e)

    def add(self, e: Edge) -> None:
        key = e.core()
        if key not in self._set:
            self._set.add(key)
            self._edges.append(e)

    def __iter__(self):
        return iter(self._edges)

    def __len__(self):
        return len(self._edges)


# ======================= 辅助函数 ========================

def _transpose_geneindex_once(geneindex: List[List[str]]) -> List[Tuple[str, ...]]:
    return list(map(tuple, zip(*geneindex)))


def _cache_bed_files(beddir: Path, sample_list: List[str]) -> Dict[str, List[Gene]]:
    cache: Dict[str, List[Gene]] = {}
    for s in sample_list:
        p = beddir / f"{s}.bed"
        if not p.exists():
            raise FileNotFoundError(f"{p} does not exist")
        cache[s] = read_bed6(str(p))
    return cache


# ===================== 核心算法函数 =====================

def get_best_chromosome(
    geneindex_transposed: List[Tuple[str, ...]],
    bed_cache: Dict[str, List[Gene]],
    sample_list: List[str]
) -> str:
    """找出在各样本中出现次数最多的染色体"""
    chrom_counter: Counter = Counter()
    for sample_idx, sample_name in enumerate(sample_list):
        bed_list = bed_cache[sample_name]
        gene2chrom = {g.name: g.chromosome for g in bed_list}
        for cell in geneindex_transposed[sample_idx]:
            if not cell or cell == "-":
                continue
            for gene in cell.split(","):
                if gene in gene2chrom:
                    chrom_counter[gene2chrom[gene]] += 1
    if not chrom_counter:
        raise ValueError("No chromosome found from geneindex and bed files.")
    return chrom_counter.most_common(1)[0][0]


def get_strand_maxcount(
    geneindex_transposed: List[Tuple[str, ...]],
    bed_cache: Dict[str, List[Gene]],
    sample_list: List[str],
    bestchromosome: str
) -> Optional[List[Optional[str]]]:
    """计算每个节点的祖先链方向"""
    per_species_strand: List[List[Optional[str]]] = []
    for sample_idx, sample_name in enumerate(sample_list):
        gene_map = {g.name: g for g in bed_cache[sample_name]}
        col = geneindex_transposed[sample_idx]
        species_strands: List[Optional[str]] = []
        for cell in col:
            if not cell or cell == "-":
                species_strands.append(None)
                continue
            strands_for_cell: List[str] = []
            for gene in cell.split(","):
                if gene in gene_map:
                    g = gene_map[gene]
                    if g.chromosome == bestchromosome:
                        strands_for_cell.append(g.strand)
            if strands_for_cell:
                species_strands.append(max(set(strands_for_cell), key=strands_for_cell.count))
            else:
                species_strands.append(None)
        per_species_strand.append(species_strands)

    all_nodes = list(map(list, zip(*per_species_strand)))
    result: List[Optional[str]] = []
    for node_strands in all_nodes:
        votes = [s for s in node_strands if s is not None]
        result.append(max(set(votes), key=votes.count) if votes else None)

    if any(s is None for s in result):
        return None
    return result


def get_edges_paths_genome(
    geneindex_transposed: List[Tuple[str, ...]],
    bed_cache: Dict[str, List[Gene]],
    sample_list: List[str],
    bestchromosome: str
) -> Tuple[Edges, Dict[str, List[str]], Dict[str, List[str]]]:
    edges = Edges()
    paths: Dict[str, List[str]] = {}
    genome: Dict[str, List[str]] = {}

    for sample_idx, sample_name in enumerate(sample_list):
        gene_map = {g.name: g for g in bed_cache[sample_name]}
        idx_genes: List[Tuple[str, int, int]] = []
        col = geneindex_transposed[sample_idx]
        for cell_idx, cell in enumerate(col, start=1):
            if not cell or cell == "-":
                continue
            for gene in cell.split(","):
                if gene in gene_map:
                    g = gene_map[gene]
                    if g.chromosome == bestchromosome:
                        idx_genes.append((gene, cell_idx, g.start))
        idx_genes.sort(key=lambda x: x[2])
        sample_path, sample_genome = [], []
        if not idx_genes:
            paths[sample_name] = []
            genome[sample_name] = []
            continue
        for a, b in zip(idx_genes, idx_genes[1:]):
            g1, id1 = a[0], a[1]
            g2, id2 = b[0], b[1]
            edges.add(Edge(id1, gene_map[g1].strand, id2, gene_map[g2].strand))
            sample_path.append(f"{id1}{gene_map[g1].strand}")
            sample_genome.append(g1)
        last_g, last_id = idx_genes[-1][0], idx_genes[-1][1]
        sample_path.append(f"{last_id}{gene_map[last_g].strand}")
        sample_genome.append(last_g)
        paths[sample_name] = sample_path
        genome[sample_name] = sample_genome

    return edges, paths, genome


def write_gfa(
    out_path: Path,
    node_num: int,
    node_name_prefix: str,
    edges: Edges,
    paths: Dict[str, List[str]],
    genome: Dict[str, List[str]],
    ancestor_strand: List[str]
) -> None:
    """输出 GFA 文件"""
    with out_path.open("w") as out:
        out.write("H\tVN:Z:1.1\n")
        for node_id in range(1, node_num + 1):
            out.write(f"S\t{node_id}\t{node_name_prefix}{node_id-1:07d}\n")
        for i in range(node_num - 1):
            out.write(f"L\t{i+1}{ancestor_strand[i]}\t{i+2}{ancestor_strand[i+1]}\t*\n")
        for e in edges:
            out.write(f"{str(e)}\n")
        for sample_name, path in paths.items():
            out.write(f"P\t{sample_name}\t{','.join(path)}\t*\n")
        for sample_name, g in genome.items():
            out.write(f"G\t{sample_name}\t{','.join(g)}\t*\n")
        out.write("P\tANCESTOR\t")
        out.write(",".join(f"{i+1}{s}" for i, s in enumerate(ancestor_strand)))
        out.write("\t*\n")


# ======================= 主入口 ========================

def main():
    parser = argparse.ArgumentParser(description="Convert gene index + BED files to GFA format.")
    parser.add_argument("-i", "--input", required=True, help="geneindex file")
    parser.add_argument("-s", "--sample", required=True, help="sample list file")
    parser.add_argument("-b", "--bed", required=True, help="BED directory (bed6 format)")
    parser.add_argument("-o", "--output", required=True, help="output GFA file")
    args = parser.parse_args()

    sample_list = read_sample(args.sample)
    geneindex = read_geneindex(args.input)
    geneindex_transposed = _transpose_geneindex_once(geneindex)
    bed_cache = _cache_bed_files(Path(args.bed), sample_list)

    bestchromosome = get_best_chromosome(geneindex_transposed, bed_cache, sample_list)
    ancestor_strand = get_strand_maxcount(geneindex_transposed, bed_cache, sample_list, bestchromosome)
    if ancestor_strand is None:
        raise ValueError("ancestor_strand is None, please check geneindex file")

    edges, paths, genome = get_edges_paths_genome(geneindex_transposed, bed_cache, sample_list, bestchromosome)

    node_num = len(geneindex)  # geneindex 行数就是节点数
    node_id_list = list(range(1, node_num + 1))
    node_name_list = [f"SG{i:07d}" for i in range(node_num)]


    with open(args.output, "w") as out:
        out.write("H\tVN:Z:1.1\n")
        # 写 S 行
        for this_node_id, this_node_name in zip(node_id_list, node_name_list):
            out.write(f"S\t{this_node_id}\t{this_node_name}\n")
        # 写祖先边（保持多列）
        for i in range(node_num - 1):
            this_node_id = node_id_list[i]
            this_strand = ancestor_strand[i]
            next_node_id = node_id_list[i + 1]
            next_strand = ancestor_strand[i + 1]
            edges.add(Edge(this_node_id, this_strand, next_node_id, next_strand, 0))
        # 写所有 edges
        for edge in edges:
            out.write(f"{edge}\n")
        # 写 paths
        for sample_name, this_path in paths.items():
            thispath = ",".join(this_path)
            out.write(f"P\t{sample_name}\t{thispath}\t*\n")
        # 写基因组信息
        for sample_name, this_genome in genome.items():
            thisgenome = ",".join(this_genome)
            out.write(f"G\t{sample_name}\t{thisgenome}\t*\n")
        # 写祖先 path
        out.write("\t".join([
            "P",
            "ANCESTOR",
            ",".join([f"{index+1}{this_strand}" for index, this_strand in enumerate(ancestor_strand)]),
            "*\n"
        ]))

    print(f"[✔] GFA successfully written to {args.output}")


if __name__ == "__main__":
    main()
