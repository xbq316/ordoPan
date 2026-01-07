#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: synpan.py
# Created time: 2025/11/18 14:53:12
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python synpan.py -h


import os
import subprocess
from rich.console import Console
from intervaltree import IntervalTree
import pandas as pd
import numpy as np
import math
from typing import Annotated 
from Bio import SeqIO
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor, as_completed
from threading import Lock
import shutil
from itertools import combinations
import sys
import tempfile
import networkx as nx
from pathlib import Path
from typing import List
from multiprocessing import Pool,Queue,Process
import argparse
from dataclasses import dataclass, field


config={
    "diamond_blastp":{"-p":1,
                      "-e":1e-20,
                      "-f":["6","qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen"],
                      "--id":40,
                       "-k":3
                      },
    "diamond_makedb":{"--threads":2,
                      },
    "diamond_path":"",
    "run_DAG_chainer_path":"",
    "perl_path":"",
    "mcl_path":"",
    "work_config":{"working_dir_name":"working",
                   "blastp_dir_name":"blastp",
                   "blastp_db_suffix":".dmnd",
                   "blastp_suffix":".blastp",
                   "chromosome_pep_dir_name":"pep",
                   "chromosome_bed_dir_name":"bed",
                   "dagchainer_input_dir_name":"dagchainer",
                   "dagchainer_output_dir_name":"dagchainer",
                   "dagchainer_input_suffix":".dagchainer",
                   "dagchainer_output_suffix":".aligncoords",
                   "geneindex_output_file_name":"geneindex.txt",
                   "geneindex_output_file_final_name":"final.geneindex.txt",
                   "RBNH_output_dir_name":"RBNH",
                   "RBNH_suffix":".RBNH",
                   "all_RBNH_file_name":"all_sample_chromosome.RBNH",
                   "mcl_output_dir_name":"mcl",
                   "mcl_output_merge_file_name":"all.cluster",
                   "mcl_input_dir_name":"mcl",
                   "mcl_input_suffix":".mcl",
                   "mcl_output_suffix":".cluster",
                   "mcl_split_geneindex_dir_name":"final_geneindex"
                   },
    "gene_filter_method":{
        "methods":["filter_synteny_blocks",
                   "filter_synteny_blocks_normal"],
        "default":"filter_synteny_blocks",
    },
    "bed_format":{"bed4":False,"bed6":True},
    "overwrite":{
        "pep_bed_split":False,
        "makedb":False,
        "blastp":False,
        "mcl":False,
        "dagchainer_input":False,
        "dagchainer_output":False,
        "geneindex_output":True,
        "RBNH_compute":False,
        "mcl_output":False,
        "merge_mcl_output":False
        },
    "mcl_config":{"-I":1.5,"--abc":""},
    "input_file_type":{"pep_suffix":".pep",
                       "bed_suffix":".bed",
                       "chromosome_file_name":"chro.list",
                       "sample_file_name":"sample.list"
                       },
    "DAGchainer_config":{
                        "-D":10,
                         "-g":1,
                         "-A":5
                         },
    "synteny_config":{"fotmat":"dagchainer"},
}
if sum([config["bed_format"]["bed4"], config["bed_format"]["bed6"]]) != 1:
    raise ValueError("Config error: choose exactly one type: bed4 or bed6.")


print_lock = Lock()   # 避免多线程输出混乱
db_lock = Lock()      # 避免多个线程重复创建相同 makedb
console=Console(soft_wrap=False)

class Tree(IntervalTree):
    """
    闭区间的 区间树
    """
    def __init__(self, intervals=None):
        # 继承 IntervalTree 的初始化
        super().__init__(intervals)

    def add_closed(self, start: int, end: int):
        """
        Add a closed interval [start, end] to the tree.

        If start is greater than end, a ValueError is raised.

        The interval is added to the tree and overlapping intervals
        are merged.

        Parameters:
        start (int): The start of the interval
        end (int): The end of the interval

        Raises:
        ValueError: If start is greater than end
        """
        if isinstance(start, int) and isinstance(end, int):
            pass
        else:
            raise ValueError("start and end must be int")            
        if start > end:
            raise ValueError("start must <= end")

        self.addi(start, end + 1)
        self.merge_overlaps()

    def remove_closed(self, start: int, end: int):
        """
        Remove a closed interval [start, end] from the tree.

        If there are existing intervals that overlap with the
        removed interval, they will be split into two intervals:

        * If the existing interval's start is less than the removed
          interval's start, the existing interval will be split into
          two intervals: [existing start, removed start] and
          [removed start, existing end].
        * If the existing interval's end is greater than the removed
          interval's end, the existing interval will be split into
          two intervals: [existing start, removed end] and
          [removed end, existing end].

        :param start: The start of the closed interval to be removed.
        :type start: int
        :param end: The end of the closed interval to be removed.
        :type end: int
        :raises ValueError: If start or end is not an integer, or
            if start is greater than end.
        """
        if isinstance(start, int) and isinstance(end, int):
            pass
        else:
            raise ValueError("start and end must be int")
        if start > end:
            raise ValueError("start must <= end")

        real_start = start
        real_end = end + 1

        to_process = sorted(self.overlap(real_start, real_end))

        for iv in to_process:
            self.remove(iv)
            if iv.begin < real_start:
                self.addi(iv.begin, real_start)
            if iv.end > real_end:
                self.addi(real_end, iv.end)

        self.merge_overlaps()

    def get_uncovered(self, start: int, end: int)-> list[tuple[int, int]]:

        """
        Return uncovered intervals in the given range.

        An uncovered interval is an interval that is not covered by
        any closed interval in the tree.

        The uncovered intervals are returned as a list of tuples,
        where each tuple represents an uncovered interval as
        (start, end).

        :param start: The start of the range.
        :type start: int
        :param end: The end of the range.
        :type end: int
        :raises ValueError: If start or end is not an integer, or
            if start is greater than end.
        :return: A list of uncovered intervals.
        :rtype: List[Tuple[int, int]]
        """
        if isinstance(start, int) and isinstance(end, int):
            pass
        else:
            raise ValueError("start and end must be int")
        if start > end:
            raise ValueError("start must <= end")

        real_start = start
        real_end = end + 1

        uncovered = []
        current = real_start

        intervals = sorted(self.overlap(real_start, real_end), key=lambda x: x.begin)

        for iv in intervals:
            if iv.begin > current:
                uncovered.append((current, iv.begin - 1))
            current = max(current, iv.end)

        if current < real_end:
            uncovered.append((current, real_end - 1))

        return uncovered

    def get_all_covered(self):
        
        """
        Get all covered intervals in the tree.

        The covered intervals are returned as a list of tuples,
        where each tuple represents a covered interval as
        (start, end).

        :return: A list of covered intervals.
        :rtype: List[Tuple[int, int]]
        """
        return [(iv.begin, iv.end - 1) for iv in sorted(self)]

class blast:
    def __init__(self, q: str, s: str, bits: float, qlen: int, slen: int, evalue: float, L_qh=None, Bprime=None):
        self.q = q
        self.s = s
        self.bits = bits
        self.qlen = qlen
        self.slen = slen
        self.evalue = evalue
        self.L_qh = L_qh
        self.Bprime = Bprime

    def compute_L_qh(self):
        # 计算并保证非零（避免 log(0)）
        val = (self.qlen or 0) * (self.slen or 0)
        if val is None or not np.isfinite(val) or val <= 0:
            val = 1e-8
        self.L_qh = val

    def compute_Bprime(self, a: float, b: float):
        # 在计算前确保 L_qh 与 bits 非零
        L = self.L_qh if (self.L_qh is not None and np.isfinite(self.L_qh) and self.L_qh > 0) else 1e-8
        bits = self.bits if (self.bits is not None and np.isfinite(self.bits) and self.bits > 0) else 1e-8
        self.Bprime = bits / (10 ** b * (L ** a))

    def __str__(self):
        return f"{self.q}\t{self.s}\t{self.Bprime}"


def read_blast_file(path: str) -> dict[frozenset[str], blast]:
    """
    读取 BLAST/DIAMOND 文件到 dict: frozenset({q,s}) -> blast
    如果在同一文件中有多条记录对应同一对 (q,s)，保留该文件内 bits 最大的那条。
    对列解析较为鲁棒：至少需要 12 列。
    """
    result: dict[frozenset[str], blast] = {}
    with open(path, 'r') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12:
                # 跳过格式不对的行（也可选择报错）
                continue
            # 标准 BLAST/DIAMOND 列索引（0-based）
            try:
                q = parts[0]
                s = parts[1]
                # bits 在索引 11，evalue 在 10
                bits = float(parts[11]) if parts[11] != '' else float('nan')
                evalue = float(parts[10]) if parts[10] != '' else float('nan')
                # qlen/slen 在 12/13（如果存在），否则用第4列 length 做备用（索引3)
                if len(parts) >= 14:
                    qlen = int(float(parts[12])) if parts[12] != '' else 0
                    slen = int(float(parts[13])) if parts[13] != '' else 0
                else:
                    # 有些 BLAST 输出没有 qlen/slen，使用 alignment length 作为保守替代
                    qlen = int(float(parts[3])) if parts[3] != '' else 0
                    slen = int(float(parts[3])) if parts[3] != '' else 0
            except Exception:
                # 跳过解析失败的行
                continue

            key = frozenset([q, s])
            obj = blast(q, s, bits, qlen, slen, evalue)
            # 如果文件内已有这对比对，保留 bits 更大的那条（同一文件去重复）
            existing = result.get(key)
            if existing is None or (np.isfinite(obj.bits) and np.isfinite(existing.bits) and obj.bits > existing.bits):
                result[key] = obj
            # 若 existing.bits 为 NaN 而 obj.bits 有效，也替换
            elif existing is not None and (not np.isfinite(existing.bits)) and np.isfinite(obj.bits):
                result[key] = obj
    return result


def merge_max_bits(blast_dict1: dict[frozenset[str], blast], blast_dict2: dict[frozenset[str], blast]) -> list[blast]:
    """
    合并两个方向（两个文件）的比对，返回 list[blast]：
    对每个无向对取两个字典中 bits 更大的那个记录（或存在的那条）。
    """
    merge_result: list[blast] = []
    keys1 = set(blast_dict1.keys())
    keys2 = set(blast_dict2.keys())
    inter = keys1 & keys2
    for k in inter:
        b1 = blast_dict1[k]
        b2 = blast_dict2[k]
        # 优先选择 bits 更大的；若 bits NaN 则选择非 NaN
        b1_bits = b1.bits if np.isfinite(b1.bits) else -np.inf
        b2_bits = b2.bits if np.isfinite(b2.bits) else -np.inf
        merge_result.append(b1 if b1_bits >= b2_bits else b2)
    for k in keys1 - keys2:
        merge_result.append(blast_dict1[k])
    for k in keys2 - keys1:
        merge_result.append(blast_dict2[k])
    return merge_result


def compute_L_qh(blastp_list: list[blast]) -> None:
    # 为每个对象计算 L_qh，并用小正数替换 0/负/NaN 保证后续 log 安全
    for obj in blastp_list:
        obj.compute_L_qh()


def bin_top_hits(blastp_list: list[blast], bin_size=1000, top_pct=0.05):
    """
    为拟合挑选 top hits：
    将 blastp_list 转为 DataFrame（显式构造），按 L_qh 分箱并在每箱取 bits 最大的 top_pct。
    返回用于拟合的 DataFrame。
    """
    # 把对象列表转成字典列表
    rows = []
    for obj in blastp_list:
        # 确保 L_qh 存在
        L = obj.L_qh if (obj.L_qh is not None and np.isfinite(obj.L_qh) and obj.L_qh > 0) else 1e-8
        bits = obj.bits if (obj.bits is not None and np.isfinite(obj.bits) and obj.bits > 0) else 1e-8
        rows.append({'q': obj.q, 's': obj.s, 'bits': bits, 'L_qh': L})
    df = pd.DataFrame(rows)
    if df.empty:
        return df  # 空 DataFrame

    df_sorted = df.sort_values('L_qh').reset_index(drop=True)
    N = len(df_sorted)
    if N < bin_size:
        chunk = max(1, int(N // 5))
    else:
        chunk = bin_size
    top_hits = []
    for i in range(0, N, chunk):
        chunk_df = df_sorted.iloc[i:i + chunk]
        k = max(1, int(math.ceil(len(chunk_df) * top_pct)))
        top_hits.append(chunk_df.nlargest(k, 'bits'))
    return pd.concat(top_hits, ignore_index=True)


def fit_log_model(sel_df: pd.DataFrame) -> tuple[float, float]:
    if sel_df is None or len(sel_df) == 0:
        return 0.0, 0.0
    x = np.log10(sel_df['L_qh'])
    y = np.log10(sel_df['bits'])
    if len(x) < 2:
        return 0.0, 0.0
    a, b = np.polyfit(x, y, 1)
    return a, b


def blastp2RBNH(blast_file1: str, blast_file2: str, out_path: str, bin_size=1000, top_pct=0.05, min_Bprime=0.0) -> str:
    """
    保持函数签名不变：读取两个 BLAST 文件，合并并计算 Bprime，写出 mcl 格式文件到 out_path
    """
    try:
        # ===== 输出文件已存在，是否跳过 =====
        if (not config["overwrite"]["RBNH_compute"]) and os.path.exists(out_path):
            return "skip"
        # ===== 输入文件合法性检查 =====
        if not (os.path.exists(blast_file1) and os.path.exists(blast_file2)):
            raise FileNotFoundError(f"BLAST file missing: {blast_file1} or {blast_file2}")
        # ===== 读取并合并 =====
        blastp_dict1 = read_blast_file(blast_file1)
        blastp_dict2 = read_blast_file(blast_file2)
        blastp_merge_list = merge_max_bits(blastp_dict1, blastp_dict2)

        # 计算 L_qh（内部做了 0->1e-8 的处理）
        compute_L_qh(blastp_merge_list)

        # 从 binned top hits 拟合 a,b
        df_binned = bin_top_hits(blastp_merge_list, bin_size, top_pct)
        a_coef, b_coef = fit_log_model(df_binned)

        # 计算 Bprime（使用 a_coef, b_coef）
        for obj in blastp_merge_list:
            obj.compute_Bprime(a_coef, b_coef)

        # 写出 mcl edge 格式：去掉 self hits，按 min_Bprime 筛选，并使用 frozenset 去重无向边
        seen = set()
        with open(out_path, 'w') as out_f:
            for obj in blastp_merge_list:
                if obj.q == obj.s:
                    continue
                key = frozenset([obj.q, obj.s])
                if key in seen:
                    continue
                # Bprime 可能为 NaN/None -> 视为 0（或你可以选择跳过）
                bprime_val = obj.Bprime if (obj.Bprime is not None and np.isfinite(obj.Bprime)) else 0.0
                if bprime_val <= min_Bprime:
                    continue
                # 为输出选择确定的 order（这里按字符串排序，保证稳定输出）
                u, v = sorted([obj.q, obj.s])
                out_f.write(f"{u} {v} {bprime_val}\n")
                seen.add(key)

        return "done"
    except Exception as e:
        raise RuntimeError(f"blastp2mcl failed ({os.path.basename(out_path)}) : {e}")
    

# def compute_RBNH(sample_list: list[str],blastp_dir: str,output_dir: str,chromosome_list: list[str],thread: int) -> None:
#     """
#     生成 MCL 输入文件(多线程) + 汇总 all.RNBH
#     """
#     os.makedirs(output_dir, exist_ok=True)
#     # ===== 准备任务列表 =====
#     tasks = []
#     RNBH_files = []   # 用于最终 all.RNBH 合并

#     # chro vs chro
#     for s1, s2 in combinations(sample_list, 2):
#         for chrom in chromosome_list:
#             blast_file1 = os.path.join(blastp_dir,f"{s1}_vs_{s2}_{chrom}{config['work_config']['blastp_suffix']}")
#             blast_file2 = os.path.join(blastp_dir,f"{s2}_vs_{s1}_{chrom}{config['work_config']['blastp_suffix']}")
#             out_path = os.path.join(output_dir,f"{s1}_vs_{s2}_{chrom}{config['work_config']['RBNH_suffix']}")
#             tasks.append((blast_file1, blast_file2, out_path))
#             RNBH_files.append(out_path)
#     # chro vs scaffold
#     for s1, s2 in combinations(sample_list, 2):
#         for chrom in chromosome_list:
#             blast_file1 = os.path.join(blastp_dir,f"{s1}_vs_{s2}_{chrom}_vs_scaffold{config['work_config']['blastp_suffix']}")
#             blast_file2 = os.path.join(blastp_dir,f"{s2}_vs_{s1}_{chrom}_vs_scaffold{config['work_config']['blastp_suffix']}")
#             out_path = os.path.join(output_dir,f"{s1}_vs_{s2}_{chrom}_vs_scaffold{config['work_config']['RBNH_suffix']}")
#             tasks.append((blast_file1, blast_file2, out_path))
#             RNBH_files.append(out_path)
#     # scaffold vs scaffold
#     for s1, s2 in combinations(sample_list, 2):
#         blast_file1 = os.path.join(blastp_dir,f"{s1}_vs_{s2}_scaffold{config['work_config']['blastp_suffix']}")
#         blast_file2 = os.path.join(blastp_dir,f"{s2}_vs_{s1}_scaffold{config['work_config']['blastp_suffix']}")
#         out_path = os.path.join(output_dir,f"{s1}_vs_{s2}_scaffold{config['work_config']['RBNH_suffix']}")
#         tasks.append((blast_file1, blast_file2, out_path))
#         RNBH_files.append(out_path)

#     # ===== 多线程执行 =====
#     with ThreadPoolExecutor(max_workers=thread) as executor:
#         futures = {
#             executor.submit(blastp2RBNH, f1, f2, op): (f1, f2, op)
#             for (f1, f2, op) in tasks
#         }
#         # ===== 输出执行状态 =====
#         for future in as_completed(futures):
#             _, _, out_path = futures[future]
#             tag = os.path.basename(out_path)
#             try:
#                 result = future.result()
#                 if result == "skip":
#                     console.print(f"[SKIP] MCL: {tag}", style="yellow")
#                 else:
#                     console.print(f"[OK] MCL: {tag}", style="green")
#             except Exception as e:
#                 console.print(f"[ERROR] MCL: {tag} -> {e}", style="red")
#     # ===== 汇总 all.RBNH =====
#     all_RNBH_files_path = os.path.join(output_dir,config["work_config"]["all_RBNH_file_name"])
#     if (not config["overwrite"]["RBNH_compute"]) and os.path.exists(all_RNBH_files_path):
#         console.print(f"[SKIP] MCL(all): {os.path.basename(all_RNBH_files_path)}", style="yellow")
#         return
#     with open(all_RNBH_files_path, "w") as outfile:
#         for f in RNBH_files:
#             if not os.path.exists(f):
#                 continue
#             if os.path.getsize(f) == 0:
#                 continue
#             with open(f, "r") as infile:
#                 for line in infile:
#                     outfile.write(line)
#     console.print(f"[OK] MCL(all): {os.path.basename(all_RNBH_files_path)}", style="green")

def compute_RBNH(sample_list: list[str], blastp_dir: str, output_dir: str,
                 chromosome_list: list[str], thread: int) -> None:
    """
    生成 RBNH（Reciprocal Best Non-Hit）输入文件 + 汇总 all.RBNH
    改进版特性：
        - 文件名构造与 DAGchainer 一致
        - 自动过滤不存在的 blast 文件
        - scaffold vs chr / chr vs scaffold 都支持
        - 更详细的日志输出
    """
    os.makedirs(output_dir, exist_ok=True)

    tasks = []
    RNBH_files = []
    total_pairs = 0

    def add_task(blast1, blast2, outpath):
        """安全添加任务：仅在文件存在时加入"""
        if not os.path.exists(blast1) or not os.path.exists(blast2):
            console.print(f"[MISS] Missing BLAST: "
                          f"{os.path.basename(blast1)} or {os.path.basename(blast2)}",
                          style="red")
            return
        tasks.append((blast1, blast2, outpath))
        RNBH_files.append(outpath)

    # ======================================================
    # 1. chr vs chr
    # ======================================================
    for s1, s2 in combinations(sample_list, 2):
        for chrom in chromosome_list:
            b1 = os.path.join(blastp_dir, f"{s1}_vs_{s2}_{chrom}{config['work_config']['blastp_suffix']}")
            b2 = os.path.join(blastp_dir, f"{s2}_vs_{s1}_{chrom}{config['work_config']['blastp_suffix']}")
            o = os.path.join(output_dir, f"{s1}_vs_{s2}_{chrom}{config['work_config']['RBNH_suffix']}")
            add_task(b1, b2, o)
            total_pairs += 1

    # ======================================================
    # 2. scaffold vs chr
    # ======================================================
    for s1, s2 in combinations(sample_list, 2):
        for chrom in chromosome_list:
            b1 = os.path.join(blastp_dir,
                              f"{s1}_vs_{s2}_scaffold_vs_{chrom}{config['work_config']['blastp_suffix']}")
            b2 = os.path.join(blastp_dir,
                              f"{s2}_vs_{s1}_{chrom}_vs_scaffold{config['work_config']['blastp_suffix']}")
            o = os.path.join(output_dir,
                             f"{s1}_vs_{s2}_scaffold_vs_{chrom}{config['work_config']['RBNH_suffix']}")
            add_task(b1, b2, o)
            total_pairs += 1

    # ======================================================
    # 3. chr vs scaffold
    # ======================================================
    for s1, s2 in combinations(sample_list, 2):
        for chrom in chromosome_list:
            b1 = os.path.join(blastp_dir,
                              f"{s1}_vs_{s2}_{chrom}_vs_scaffold{config['work_config']['blastp_suffix']}")
            b2 = os.path.join(blastp_dir,
                              f"{s2}_vs_{s1}_scaffold_vs_{chrom}{config['work_config']['blastp_suffix']}")
            o = os.path.join(output_dir,
                             f"{s1}_vs_{s2}_{chrom}_vs_scaffold{config['work_config']['RBNH_suffix']}")
            add_task(b1, b2, o)
            total_pairs += 1

    # ======================================================
    # 4. scaffold vs scaffold
    # ======================================================
    for s1, s2 in combinations(sample_list, 2):
        b1 = os.path.join(blastp_dir, f"{s1}_vs_{s2}_scaffold{config['work_config']['blastp_suffix']}")
        b2 = os.path.join(blastp_dir, f"{s2}_vs_{s1}_scaffold{config['work_config']['blastp_suffix']}")
        o = os.path.join(output_dir, f"{s1}_vs_{s2}_scaffold{config['work_config']['RBNH_suffix']}")
        add_task(b1, b2, o)
        total_pairs += 1

    console.print(f"[INFO] Total RBNH pairs prepared: {total_pairs}", style="cyan")
    console.print(f"[INFO] Valid RBNH tasks: {len(tasks)}", style="cyan")

    # ======================================================
    # 运行多线程
    # ======================================================
    success = 0
    skipped = 0

    with ThreadPoolExecutor(max_workers=thread) as executor:
        futures = {
            executor.submit(blastp2RBNH, f1, f2, op): (f1, f2, op)
            for (f1, f2, op) in tasks
        }

        for future in as_completed(futures):
            _, _, out_path = futures[future]
            tag = os.path.basename(out_path)
            try:
                result = future.result()
                if result == "skip":
                    skipped += 1
                    console.print(f"[SKIP] RBNH: {tag}", style="yellow")
                else:
                    success += 1
                    console.print(f"[OK] RBNH: {tag}", style="green")
            except Exception as e:
                console.print(f"[ERROR] RBNH: {tag} -> {e}", style="red")

    console.print(f"[INFO] RBNH Done: {success} OK, {skipped} skipped", style="cyan")

    # ======================================================
    # 合并 all.RBNH
    # ======================================================
    all_RNBH_files_path = os.path.join(output_dir, config["work_config"]["all_RBNH_file_name"])

    if (not config["overwrite"]["RBNH_compute"]) and os.path.exists(all_RNBH_files_path):
        console.print(f"[SKIP] RBNH(all): {os.path.basename(all_RNBH_files_path)}", style="yellow")
        return

    with open(all_RNBH_files_path, "w") as outfile:
        for f in RNBH_files:
            if not os.path.exists(f) or os.path.getsize(f) == 0:
                continue
            with open(f, "r") as infile:
                for line in infile:
                    outfile.write(line)

    console.print(f"[OK] RBNH(all): {os.path.basename(all_RNBH_files_path)}", style="green")



@dataclass
class Block:
    chromosome1: str = ""
    genelist1: list[str] = field(default_factory=list)
    positions1: list[int] = field(default_factory=list)

    chromosome2: str = ""
    genelist2: list[str] = field(default_factory=list)
    positions2: list[int] = field(default_factory=list)

    hits: int = 0
    size1: int = 0
    size2: int = 0
    reverse1: bool = False
    reverse2: bool = False

    def __post_init__(self):
        # 长度一致检查
        a = [len(self.genelist1), len(self.genelist2),
             len(self.positions1), len(self.positions2)]
        assert all(i == a[0] for i in a)

        self.hits = a[0]

        if self.hits == 0:
            return

        # 计算 size 与方向
        self.size1 = abs(self.positions1[-1] - self.positions1[0]) + 1
        self.size2 = abs(self.positions2[-1] - self.positions2[0]) + 1
        self.reverse1 = self.positions1[0] > self.positions1[-1]
        self.reverse2 = self.positions2[0] > self.positions2[-1]

    def __str__(self):
        return f"qry-chr:{self.chromosome1},qry-size:{self.size1},ref-chr:{self.chromosome2},ref-size:{self.size2},hits:{self.hits}"

    def __repr__(self):
        return str(self)

    def intersection(self,left_segments: list[tuple[int,int]],right_segments: list[tuple[int,int]])->"Block":
        # if self.is_normal()==False:
        #     raise ValueError("Block is not normal")
        
        left_segments_len=len(left_segments)
        right_segments_len=len(right_segments)
        if left_segments_len == 0 or right_segments_len==0:
            return Block(self.chromosome1,[],[],self.chromosome2,[],[])
        left_geneorder2index={position:index for index,position in enumerate(self.positions1)}
        right_geneorder2index={position:index  for index,position in enumerate(self.positions2)}
        
        left_index_list=[]
        right_index_list=[]
        for start,end in left_segments:
            start_index=left_geneorder2index[start]
            end_index=left_geneorder2index[end]
            left_index_list.extend(list(range(start_index,end_index+1)))
        for start,end in right_segments:
            start_index=right_geneorder2index[start]
            end_index=right_geneorder2index[end]
            right_index_list.extend(list(range(start_index,end_index+1)))
        left_index_set=set(left_index_list)
        right_index_set=set(right_index_list)
        final_index_list=sorted(left_index_set.intersection(right_index_set))
        new_genelist1=[self.genelist1[i] for i in final_index_list]
        new_positions1=[self.positions1[i] for i in final_index_list]
        new_genelist2=[self.genelist2[i] for i in final_index_list]
        new_positions2=[self.positions2[i] for i in final_index_list]
        return Block(self.chromosome1,new_genelist1,new_positions1,self.chromosome2,new_genelist2,new_positions2)
    
    def adjust_segments(self,left_segments: list[tuple[int,int]],right_segments: list[tuple[int,int]])->tuple[list[tuple[int,int]],list[tuple[int,int]]]:
        new_left_segments=[]
        for start,end in left_segments:
            a=[]
            for i in self.positions1:
                if i>=start and i<=end:
                    a.append(i)
            if len(a)==0:
                continue
            start=min(a)
            end=max(a)
            new_left_segments.append((start,end))
        new_right_segments=[]
        for start,end in right_segments:
            a=[]
            for i in self.positions2:
                if i>=start and i<=end:
                    a.append(i)
            if len(a)==0:
                continue
            start=min(a)
            end=max(a)
            new_right_segments.append((start,end))
        return new_left_segments,new_right_segments
    
    def is_normal(self)->bool:
        a=[len(self.genelist1),len(self.genelist2),len(self.positions1),len(self.positions2)]
        if all(i==a[0] for i in a):
            return True 
        else:
            return False

    def length(self)->int:
        is_normal=self.is_normal()
        if not is_normal:
            raise ValueError("Block is not normal")
        return self.genelist1.__len__()


class Gene:
    def __init__(self,chromosome: str,start: int,end: int,name: str,score:str=".",strand: str="+",index=None):
        """
        Initialize a Gene object
        
        Parameters:
        chromosome (str): The chromosome where the gene is located
        start (int): The start position of the gene
        end (int): The end position of the gene
        name (str): The name of the gene
        score (str): The score of the gene
        strand (str): The strand of the gene
        index (int): The index of the gene

        Returns:
        None
        """
        self.chromosome=chromosome
        self.start=start
        self.end=end
        self.name=name
        self.score=score
        self.strand=strand
        self.index=index

    def __str__(self):
        return f"{self.chromosome},{self.start},{self.end},{self.name},{self.score},{self.strand}"

    def __repr__(self):
        return str(self)

def read_bed(filename:str,)-> list[Gene]:
    """
    Read a BED file and return a list of Gene objects
    
    Parameters:
    filename (str): The name of the BED file to read
    
    Returns:
    list[Gene]: A list of Gene objects representing the contents of the BED file
    """
    result=[]
    if config["bed_format"]["bed6"]:
        with open(filename) as bed:
            for i in bed:
                i=i.strip()
                if not i:
                    continue
                line=i.split("\t")
                try:
                    chromosome,start,end,name,score,strand=line
                    start=int(start)
                    end=int(end)
                except :
                    raise ValueError(f"{line} format error in {filename}")
                result.append(Gene(chromosome,start,end,name,score,strand))
    if config["bed_format"]["bed4"]:
        with open(filename) as bed:
            for i in bed:
                i=i.strip()
                if not i:
                    continue
                line=i.split("\t")
                try:
                    chromosome,start,end,name=line
                    start=int(start)
                    end=int(end)
                except :
                    raise ValueError(f"{line} format error in {filename}")
                result.append(Gene(chromosome,start,end,name))
    result.sort(key=lambda x:(x.chromosome,x.start))

    thisChromosome=None
    thisIndex=1
    for gene in result:
        if thisChromosome!=gene.chromosome:
            thisChromosome=gene.chromosome
            thisIndex=1
        gene.index=thisIndex
        thisIndex+=1
    return result

def genelist2info(genelist: list[Gene])-> dict[str,Gene]:
    result={}
    for gene in genelist:
        result[gene.name]=gene
    return result


def run_diamond_makedb(diamond_path: str, ref_pep: str, db_name: str,args: list[str]) -> str:
    # 避免多个线程同时进入 makedb 创建相同 db
    with db_lock:
        if not config["overwrite"]["makedb"]:
            if os.path.exists(db_name):
                with print_lock:
                    print(f"[SKIP] DB exists: {db_name}")
                return db_name
        
        cmd_list = [diamond_path, "makedb", "--in", ref_pep, "-d", db_name, *args]
        with print_lock:
            print(f"[RUN] {' '.join(cmd_list)}")
    # 创建数据库
    result = subprocess.run(cmd_list, text=True)
    if result.returncode == 0:
        with print_lock:
            print(f"[OK] Created DB: {db_name}")
    else:
        with print_lock:
            print(f"[ERROR] Failed to create {db_name}")
    return db_name



def run_diamond_blastp(diamond_path: str, qry_pep: str, db_name: str,
                       output: str, args: list[str]) -> str:
    """Return: 'done' or 'skip' """
    # 如果输出文件已存在，直接跳过
    if not config["overwrite"]["blastp"]:
        if os.path.exists(output):
            # with print_lock:
            #     console.print(f"[SKIP] Output already exists: {output}", style="yellow")
            return "skip"
    cmd_list = [diamond_path, "blastp", "-q", qry_pep, "-d", db_name, *args, "-o", output]
    with print_lock:
        console.print(f"[RUN] {' '.join(cmd_list)}", style="blue")
    result = subprocess.run(cmd_list, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        # with print_lock:
            # console.print(f"[ERROR] blastp failed for {qry_pep} to {db_name}", style="red")
            # console.print(f"{result.stderr}", style="red")
        return "error"
    # with print_lock:
    #     console.print(f"[OK] blastp finished: {output}", style="green")
    return "done"

def run_DAG_chainer_pl(perl_path: str, dagchainer_path: str,
                       input_dagchainer: str, output: str,
                       args: list[str]) -> str:
    """
    Run DAGchainer perl script safely in a temporary directory.
    Returns: 'done' | 'skip' | 'error'
    """
    # 如果输出已经存在则跳过
    if not config["overwrite"]["dagchainer_output"] and os.path.exists(output):
        with print_lock:
            console.print(f"[SKIP] Output exists: {output}", style="yellow")
        return "skip"

    # 创建独立临时目录
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_input = os.path.join(tmpdir, os.path.basename(input_dagchainer))
        shutil.copy(input_dagchainer, tmp_input)

        # 构建命令
        cmd_list = [perl_path, dagchainer_path, "-i", tmp_input, *args]

        with print_lock:
            console.print(f"[RUN] {' '.join(cmd_list)}", style="blue")

        # 在临时目录执行
        result = subprocess.run(cmd_list, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tmpdir)

        if result.returncode != 0:
            with print_lock:
                console.print(f"[ERROR] run_DAG_chainer.pl failed: {input_dagchainer}", style="red")
                console.print(f"{result.stderr}", style="red")
            return "error"

        # 找到输出文件并移动到目标位置
        # 假设 DAGchainer 会生成 <input_basename>.dagchainer.aligncoords
        tmp_output_file = os.path.join(tmpdir, os.path.basename(output))
        if os.path.exists(tmp_output_file):
            shutil.move(tmp_output_file, output)
        else:
            # 如果找不到输出文件，也视为错误
            with print_lock:
                console.print(f"[ERROR] DAGchainer did not produce expected output: {output}", style="red")
            return "error"

    with print_lock:
        console.print(f"[OK] run_DAG_chainer.pl finished: {output}", style="green")
    return "done"



def read_synteny_dagchainer(filename: str, bed_dict_left: dict[str, "Gene"] = None, bed_dict_right: dict[str, "Gene"] = None)-> list[Block]:
    result = []
    for i in open(filename):
        i = i.strip()
        if not i:
            continue
        if i.startswith("#"):
            # 新的 synteny block
            result.append([])
        else:
            line = i.split("\t")
            chrom_left, gene_left = line[0], line[1]
            chrom_right, gene_right = line[4], line[5]
            # 添加到最后一个 block
            if not result:
                result.append([])  # 确保至少有一个 block
            result[-1].append([chrom_left, gene_left, chrom_right, gene_right])
    # 检查并根据基因来源决定是否交换
    for this_block in result:
        if not this_block:
            continue
        chrom_left, gene_left, chrom_right, gene_right = this_block[0]

        swap_flag = False
        if gene_left in bed_dict_left and gene_right in bed_dict_right:
            pass
        elif gene_left in bed_dict_right and gene_right in bed_dict_left:
            swap_flag = True
        else:
            raise KeyError(f"Gene {gene_left} or {gene_right} not found in bed file")
        if swap_flag:
            for this_line in this_block:
                this_line[0], this_line[2] = this_line[2], this_line[0]
                this_line[1], this_line[3] = this_line[3], this_line[1]
    
    final_result = []
    for this_block in result:
        chromosome1, chromosome2 = this_block[0][0], this_block[0][2]
        genelist1=[]
        positions1=[]
        genelist2=[]
        positions2=[]
        for line in this_block:
            genelist1.append(bed_dict_left[line[1]].name)
            positions1.append(bed_dict_left[line[1]].index)
            genelist2.append(bed_dict_right[line[3]].name)
            positions2.append(bed_dict_right[line[3]].index)
        final_result.append(Block(chromosome1,genelist1,positions1,chromosome2,genelist2,positions2))
    return final_result

def filter_synteny_blocks(synteny_blocks: list[Block],chromosome_list: list[str])-> list[tuple[str,str]]:#list[tuple[str,str,str,list[tuple[str,str]]]]:
    # sample1 、 sample2 、 chromosome 、 pairs
    chromosome_dict={chromosome:index for index,chromosome in enumerate(chromosome_list)}
    synteny_blocks.sort(key=lambda x: (
        not (x.chromosome1 in chromosome_dict and x.chromosome2 in chromosome_dict),  # 满足放前面
        chromosome_dict.get(x.chromosome1,float("inf")),
        chromosome_dict.get(x.chromosome2,float("inf")),
        -x.hits
        )
    )
    result_pair=[]
    block_tree_left={}
    block_tree_right={}

    for block_index,this_block in enumerate(synteny_blocks):
        if this_block.chromosome1 not in block_tree_left:
            block_tree_left[this_block.chromosome1]=Tree()
        if this_block.chromosome2 not in block_tree_right:
            block_tree_right[this_block.chromosome2]=Tree()

        this_block_left_start=min(this_block.positions1)
        this_block_left_end=max(this_block.positions1)
        this_block_right_start=min(this_block.positions2)
        this_block_right_end=max(this_block.positions2)
        
        if block_index==0:
            block_tree_left[this_block.chromosome1].add_closed(this_block_left_start, this_block_left_end)
            block_tree_right[this_block.chromosome2].add_closed(this_block_right_start, this_block_right_end)
            for gene1,gene2 in zip(this_block.genelist1,this_block.genelist2):
                result_pair.append((gene1,gene2))
            continue
        left_segments_uncovered = block_tree_left[this_block.chromosome1].get_uncovered(this_block_left_start, this_block_left_end)
        right_segments_uncovered = block_tree_right[this_block.chromosome2].get_uncovered(this_block_right_start, this_block_right_end)

        new_left_segments_uncovered,new_right_segments_uncovered=this_block.adjust_segments(left_segments_uncovered,right_segments_uncovered)

        new_block=this_block.intersection(new_left_segments_uncovered,new_right_segments_uncovered)
        block_tree_left[this_block.chromosome1].add_closed(this_block_left_start, this_block_left_end)
        block_tree_right[this_block.chromosome2].add_closed(this_block_right_start, this_block_right_end)
        if new_block.length()==0:
            continue
        for gene1,gene2 in zip(new_block.genelist1,new_block.genelist2):
            result_pair.append((gene1,gene2))
    return result_pair

def filter_synteny_blocks_normal(synteny_blocks: list[Block],chromosome_list: list[str])-> list[tuple[str,str]]:
    chromosome_dict={chromosome:index for index,chromosome in enumerate(chromosome_list)}
    synteny_blocks.sort(key=lambda x: (
        not (x.chromosome1 in chromosome_dict and x.chromosome2 in chromosome_dict),  # 满足放前面
        chromosome_dict.get(x.chromosome1,float("inf")),
        chromosome_dict.get(x.chromosome2,float("inf")),
        -x.hits
        )
    )
    left_gene_set=set()
    right_gene_set=set()

    result_pair=[]
    for this_block in synteny_blocks:
        if not this_block.is_normal():
            continue
        for gene1,gene2 in zip(this_block.genelist1,this_block.genelist2):
            if gene1 in left_gene_set or gene2 in right_gene_set:
                continue
            left_gene_set.add(gene1)
            right_gene_set.add(gene2)
            result_pair.append((gene1,gene2))
    return result_pair

def sample1_vs_sample2_pair(sample1:str,sample2:str,
                            bed_dir: str,
                            chromosome_list: list[str],
                            dagchainer_output_dir: str,
                            filter_method: str = config["gene_filter_method"]["default"],
                            )->tuple[str, str, list[tuple[str, str]]]:
    
    filtered_func = globals().get(filter_method)
    if filtered_func is None:
        raise ValueError(f"Unknown filter method: {filter_method}")
    
    ref_bed_list=[]
    qry_bed_list=[] 
    ref_blocks_list=[]
    qry_blocks_list=[]
    for chro in chromosome_list:
        bed_file1=os.path.join(bed_dir,f"{sample1}_{chro}.bed")
        bed_file2=os.path.join(bed_dir,f"{sample2}_{chro}.bed")
        ref_bed_list.extend(read_bed(bed_file1))
        qry_bed_list.extend(read_bed(bed_file2))
    bed_file1=os.path.join(bed_dir,f"{sample1}_scaffold.bed")
    bed_file2=os.path.join(bed_dir,f"{sample2}_scaffold.bed")
    ref_bed_list.extend(read_bed(bed_file1))
    qry_bed_list.extend(read_bed(bed_file2))

    ref_bed_dict=genelist2info(ref_bed_list)
    qry_bed_dict=genelist2info(qry_bed_list)
    
    for chro in chromosome_list:
        # chro vs chro 
        dagchainer_file1=os.path.join(dagchainer_output_dir,f"{sample1}_vs_{sample2}_{chro}{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
        dagchainer_file2=os.path.join(dagchainer_output_dir,f"{sample2}_vs_{sample1}_{chro}{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
        ref_blocks_list.extend(read_synteny_dagchainer(dagchainer_file1,bed_dict_left=ref_bed_dict,bed_dict_right=qry_bed_dict))
        qry_blocks_list.extend(read_synteny_dagchainer(dagchainer_file2,bed_dict_left=ref_bed_dict,bed_dict_right=qry_bed_dict))
        
        # chro vs scaffold
        dagchainer_file1=os.path.join(dagchainer_output_dir,f"{sample1}_vs_{sample2}_{chro}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
        dagchainer_file2=os.path.join(dagchainer_output_dir,f"{sample2}_vs_{sample1}_{chro}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
        ref_blocks_list.extend(read_synteny_dagchainer(dagchainer_file1,bed_dict_left=ref_bed_dict,bed_dict_right=qry_bed_dict))
        qry_blocks_list.extend(read_synteny_dagchainer(dagchainer_file2,bed_dict_left=ref_bed_dict,bed_dict_right=qry_bed_dict))

        # scaffold vs chro
        dagchainer_file1=os.path.join(dagchainer_output_dir,f"{sample1}_vs_{sample2}_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
        dagchainer_file2=os.path.join(dagchainer_output_dir,f"{sample2}_vs_{sample1}_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
        ref_blocks_list.extend(read_synteny_dagchainer(dagchainer_file1,bed_dict_left=ref_bed_dict,bed_dict_right=qry_bed_dict))
        qry_blocks_list.extend(read_synteny_dagchainer(dagchainer_file2,bed_dict_left=ref_bed_dict,bed_dict_right=qry_bed_dict))

    # scaffold vs scaffold
    dagchainer_file1=os.path.join(dagchainer_output_dir,f"{sample1}_vs_{sample2}_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
    dagchainer_file2=os.path.join(dagchainer_output_dir,f"{sample2}_vs_{sample1}_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
    ref_blocks_list.extend(read_synteny_dagchainer(dagchainer_file1,bed_dict_left=ref_bed_dict,bed_dict_right=qry_bed_dict))
    qry_blocks_list.extend(read_synteny_dagchainer(dagchainer_file2,bed_dict_left=ref_bed_dict,bed_dict_right=qry_bed_dict))

    filtered_blocks = filtered_func([*ref_blocks_list, *qry_blocks_list],chromosome_list)
    return (sample1, sample2, filtered_blocks)


def blastp2dagchainer(blast_file: str, 
                      sample1: str, sample2: str,
                      bed_dict1: dict[str,Gene], bed_dict2: dict[str,Gene],
                      out_path: str):
    """
    Convert BLASTP results to DAGchainer input format.
    Output columns:
    qry_chr qry_gene qry_start qry_end ref_chr ref_gene ref_start ref_end evalue
    """
    blast_dict = read_blast_file(blast_file)  # dict[frozenset, blast]

    # 将 dict 转为 list 以模仿 df.itertuples()
    rows = list(blast_dict.values())

    with open(out_path, "w") as out:
        for row in rows:
            qgene = row.q
            sgene = row.s

            if qgene not in bed_dict1 or sgene not in bed_dict2:
                continue

            qbed = bed_dict1[qgene]
            sbed = bed_dict2[sgene]

            new_qbed_chromosome = f"{sample1}-{qbed.chromosome}"
            new_sbed_chromosome = f"{sample2}-{sbed.chromosome}"

            out.write(
                f"{new_qbed_chromosome}\t{qgene}\t{qbed.index}\t{qbed.index}\t"
                f"{new_sbed_chromosome}\t{sgene}\t{sbed.index}\t{sbed.index}\t"
                f"{row.evalue}\n"
            )


def get_dagchainer_input(blast_file1: str, blast_file2: str,
                         sample1:str,sample2:str,
                         bed_file1: str, bed_file2: str,
                         out_path1: str, out_path2: str):
    if not config["overwrite"]["dagchainer_input"]:
        if os.path.exists(out_path1) and os.path.exists(out_path2):
            console.print(f"[SKIP] run_DAG_chainer.pl input exists: {out_path1}", style="yellow")
            console.print(f"[SKIP] run_DAG_chainer.pl input exists: {out_path2}", style="yellow")
            return "skip"

    bed_dict1 = {g.name: g for g in read_bed(bed_file1)}
    bed_dict2 = {g.name: g for g in read_bed(bed_file2)}

    blastp2dagchainer(blast_file1,sample1,sample2,bed_dict1, bed_dict2, out_path1)
    blastp2dagchainer(blast_file2,sample2,sample1,bed_dict2, bed_dict1, out_path2)
    return "done"


def read_pep(pep_file:str):
    pep_dict={}
    for record in SeqIO.parse(pep_file,"fasta"):
        pep_dict[record.id]=str(record.seq)
    return pep_dict


def split_pep_bed(bed_file: str, pep_file: str, pep_dir: str, bed_dir: str,
                  chromosome_list: list[str], sample_name: str) -> str:
    chromosome_set = set(chromosome_list)

    # 判断是否已存在全部输出文件 → 跳过
    if not config["overwrite"]["pep_bed_split"]: # 如果不允许覆盖
        pep_exist = all(os.path.exists(os.path.join(pep_dir, f"{sample_name}_{chromosome}.pep"))
                        for chromosome in chromosome_set)
        bed_exist = all(os.path.exists(os.path.join(bed_dir, f"{sample_name}_{chromosome}.bed"))
                        for chromosome in chromosome_set)
        if pep_exist and bed_exist:
            return "skip"
        
    # 解析数据
    pep_data = read_pep(pep_file)
    bed_data = read_bed(bed_file)
    bed2info = {gene.name: gene for gene in bed_data}
    pep_list = defaultdict(list)
    for genename, pep in pep_data.items():
        if genename not in bed2info:
            raise KeyError(f"Gene {genename} not found in bed file {bed_file}")
        chromosome = bed2info[genename].chromosome
        new_genename = f"{sample_name}_{genename}"
        pep_list[chromosome].append((new_genename, pep))

    # 写 scaffold pep
    pep_buffers = defaultdict(list)
    scaffold_pep = []

    for chromosome, genelist in pep_list.items():
        if chromosome not in chromosome_set:
            for genename, pep in genelist:
                scaffold_pep.append(f">{genename}\n{pep}\n")
        else:
            for genename, pep in genelist:
                pep_buffers[chromosome].append(f">{genename}\n{pep}\n")

    # 一次性写 pep
    with open(os.path.join(pep_dir, f"{sample_name}_scaffold.pep"), "w") as out:
        out.writelines(scaffold_pep)
    for chrom, lines in pep_buffers.items():
        with open(os.path.join(pep_dir, f"{sample_name}_{chrom}.pep"), "w") as out:
            out.writelines(lines)
    # 缓存 bed 输出
    bed_buffers = defaultdict(list)
    scaffold_bed = []
    if config["bed_format"]["bed6"]:
        for gene in bed_data:
            new_name = f"{sample_name}_{gene.name}"
            line = f"{gene.chromosome}\t{gene.start}\t{gene.end}\t{new_name}\t{gene.score}\t{gene.strand}\n"
            if gene.chromosome not in chromosome_set:
                scaffold_bed.append(line)
            else:
                bed_buffers[gene.chromosome].append(line)
    if config["bed_format"]["bed4"]:
        for gene in bed_data:
            new_name = f"{sample_name}_{gene.name}"
            line = f"{gene.chromosome}\t{gene.start}\t{gene.end}\t{new_name}\n"
            if gene.chromosome not in chromosome_set:
                scaffold_bed.append(line)
            else:
                bed_buffers[gene.chromosome].append(line)
    # 一次性写 bed
    with open(os.path.join(bed_dir, f"{sample_name}_scaffold.bed"), "w") as out:
        out.writelines(scaffold_bed)
    for chrom, lines in bed_buffers.items():
        with open(os.path.join(bed_dir, f"{sample_name}_{chrom}.bed"), "w") as out:
            out.writelines(lines)
    return "done"


def combine_SG(data: list[tuple[str, str, list[tuple[str, str]]]],chromosome_list: list[str],sample_list: list[str],bed_dir: str)-> list[tuple[str,str]]:
    # sample1 , sample2 , pairs
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    outfile_path=os.path.join(output_dir,f"{config['work_config']['geneindex_output_file_name']}")
    if not config["overwrite"]["geneindex_output"]:
        if os.path.exists(outfile_path):
            return []
    
    sample2index={sample:i for i,sample in enumerate(sample_list)}
    sample_size=len(sample_list)

    G=nx.Graph()
    for sample1,sample2,pairs in data:
        for gene1,gene2 in pairs:
            G.add_node(gene1,sample=sample1)
            G.add_node(gene2,sample=sample2)
            G.add_edge(gene1,gene2)
    visited=set()
    SG_list=[]
    result_edges=[]

    for node in G.nodes():
        # 深度优先搜索树
        if node in visited:
            continue
        tree = nx.dfs_tree(G, node)
        path_nodes = list(tree.nodes)
        path_edges=list(tree.edges)
        
        visited.update(path_nodes)

        SG_line=[[] for _ in range(sample_size)]
        for gene in path_nodes:
            sample = sample2index[G.nodes[gene]["sample"]]
            SG_line[sample].append(gene)
        
        flag_duplicate=False
        a=[]
        for i in SG_line:
            if not i:
                a.append("-")
            else:
                a.append(",".join(i))
                if len(i)>1:
                    flag_duplicate=True
        SG_list.append(a) # matched gene
        if flag_duplicate:
            result_edges.append(path_edges)

    # 对没有 blast 上的基因进行填充,也就是blastp能匹配上的都是有同源基因的，但是有些是没有同源基因的私有基因
    for sample,this_sample_gene_list in zip(sample_list,list(zip(*SG_list))):
        this_sample_gene_set=set() # 记录 matched gene
        for genes in this_sample_gene_list:
            if genes=="-":
                continue
            elif "," in genes:
                for gene in genes.split(","):
                    this_sample_gene_set.add(gene)
            else:
                this_sample_gene_set.add(genes)
        
        bed_data=[]
        for chromosome in chromosome_list:
            bed_path=os.path.join(bed_dir,f"{sample}_{chromosome}.bed")
            this_bed_data=read_bed(bed_path)
            bed_data.extend(this_bed_data)
        for gene in bed_data:
            if gene.name not in this_sample_gene_set:
                this_line=["-" for _ in range(sample_size)]
                this_line[sample2index[sample]]=gene.name
                SG_list.append(this_line) # unmatched gene


    # 写文件
    with open(outfile_path,"w") as out:
        for line in SG_list:
            out.write("\t".join(line)+"\n")
    return result_edges


def run_mcl(mcl_path: str, mcl_input_path: str, mcl_output_path: str, agrs: list[str]) -> str:
    """
    运行 MCL 聚类
    mcl_path         : mcl可执行路径或命令名 (如 "mcl")
    mcl_input_path   : 输入文件 (.mci / ABC格式)
    mcl_output_path  : 输出文件路径
    agrs             : 参数列表，如 ["--abc", "-I", "1.5"]
                       ⚠ 不修改参数设计，保持原样
    """
    # 跳过已有输出
    if (not config["overwrite"]["mcl_output"]) and os.path.exists(mcl_output_path):
        return "skip"
    mcl_args = " ".join(str(x) for x in agrs)
    cmd = f"{mcl_path} {mcl_input_path} {mcl_args} -o {mcl_output_path}"
    try:
        with print_lock:
            console.print(f"[RUN] {cmd}", style="blue")

        result = subprocess.run(
            cmd, shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        # 输出命令行内容
        with print_lock:
            if result.returncode == 0:
                console.print(result.stdout, style="green")
            else:
                console.print(result.stdout, style="red")
        return "done" if result.returncode == 0 else "error"
    except Exception as e:
        raise RuntimeError(f"MCL execution failed: {e}")

def read_mcl_cluster(mcl_output_path:str)->dict[str,int]:
    gene2cluster_dict={}
    with open(mcl_output_path) as f:
        for cluster_index,line in enumerate(f):
            line=line.strip()
            if not line:
                continue
            genes=line.split("\t")
            for gene in genes:
                gene2cluster_dict[gene]=cluster_index+1
    return gene2cluster_dict


def _process_cluster_chunk(all_RBNH_file_path: str,
                           clusters_chunk: list[list[tuple[str, str]]],
                           out_file_dir: str,
                           chunk_index: int) -> tuple[list[str], list[str]]:

    G = nx.Graph()
    with open(all_RBNH_file_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            node1, node2, weight = line.split(" ")
            G.add_edge(node1, node2, weight=float(weight))

    results = []
    logs = []

    base_idx = chunk_index * len(clusters_chunk)

    for i, cluster in enumerate(clusters_chunk):
        cluster_id = f"cluster-{base_idx + i + 1}"
        out_file = os.path.join(out_file_dir, f"{cluster_id}{config['work_config']['mcl_input_suffix']}")

        with open(out_file, "w") as f:
            for node1, node2 in cluster:
                w = G[node1][node2]["weight"]
                f.write(f"{node1} {node2} {w}\n")

        results.append(out_file)
        logs.append(f"[Chunk {chunk_index}] processed cluster {i+1}/{len(clusters_chunk)}")

    return results, logs


def get_duplicate_RBNH_parallel(all_RBNH_file_path, clusters_list, out_file_dir, n_processes):
    if not os.path.exists(out_file_dir):
        os.makedirs(out_file_dir)

    chunk_size = (len(clusters_list) + n_processes - 1) // n_processes
    chunks = [clusters_list[i:i+chunk_size] for i in range(0, len(clusters_list), chunk_size)]

    results = []
    total_chunks = len(chunks)

    queue = Queue()

    with ProcessPoolExecutor(max_workers=n_processes) as ex:
        futures = {ex.submit(_process_cluster_chunk, all_RBNH_file_path, chunk, out_file_dir, idx): idx
                for idx, chunk in enumerate(chunks)}
        completed_chunks = 0
        for future in as_completed(futures):
            files, logs = future.result()   # ✅ 现在返回 tuple[list[str], list[str]]
            results.extend(files)
            completed_chunks += 1
            for log in logs:
                console.print(log)
            console.print(f"[Progress] {completed_chunks}/{total_chunks} chunks completed "
                        f"({completed_chunks/total_chunks*100:.1f}%)")


    return results


def run_mcl_split(gene2cluster_dict: dict[str, int], geneindex_in_path: str, geneindex_out_path: str) -> None:
    with open(geneindex_in_path) as f, open(geneindex_out_path, "w") as fw:
        for line in f:
            line = line.strip()
            if not line:
                continue

            columns = line.split("\t")
            # 判断是否有重复基因
            if not any("," in genes for genes in columns):
                fw.write(line + "\n")
                continue

            cluster2gene = defaultdict(lambda: defaultdict(list))
            for idx, genes in enumerate(columns):
                if genes == "-":
                    continue
                elif "," in genes:
                    for gene in genes.split(","):
                        cluster_index = gene2cluster_dict.get(gene,gene)
                        cluster2gene[cluster_index][idx].append(gene)
                else:
                    cluster_index = gene2cluster_dict.get(genes,genes)
                    cluster2gene[cluster_index][idx].append(genes)
            all_line_write=[]
            for cluster_index, sg_dict in cluster2gene.items():
                all_line = ["-" for _ in range(len(columns))]
                for idx, gene_list in sg_dict.items():
                    all_line[idx] = ",".join(gene_list)
                # fw.write("\t".join(all_line) + "\n")
                all_line_write.append("\t".join(all_line))
            fw.write("\n".join(all_line_write) + "\n")


# def run_all_mcl_split(mcl_output_path: str,geneindex_output_dir: str,output_dir: str,thread: int) -> None:
#     os.makedirs(output_dir, exist_ok=True)
#     gene2cluster_dict = read_mcl_cluster(mcl_output_path)
    
#     geneindex_in_path = os.path.join(geneindex_output_dir,config['work_config']['geneindex_output_file_name'])
#     geneindex_out_path = os.path.join(output_dir,config['work_config']['geneindex_output_file_name'])
#     run_mcl_split(gene2cluster_dict, geneindex_in_path, geneindex_out_path)
#     console.print(f"[OK] mcl split finished: {geneindex_out_path}", style="green")


if __name__ == "__main__":

    helptext=r"""
    python synpan.py -i input_dir -o output_dir
    You need to provide four types of documents.

    Input file requirements: 
    1.amino acid sequence file(.pep): file1.pep file2.pep file3.pep ...
    >gene1
    MKTAYIAKQRQISFVKSHFSRQDILDLI
    >gene2
    GKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSRHPGNFGAF


    2.bed file: file1.bed file2.bed file3.bed ...
    bed6 :
    {chromosome}\t{start}\t{end}\t{genename}\t{socre}\t{strand}
    Chr1	6856	8726	01.col_col_AT1G01010	.	+
    Chr1	10012	11763	01.col_col_AT1G01020	.	-
    Chr1	14961	16037	01.col_col_AT1G01030	.	-
    bed4 :
    Chr1	6856	8726	01.col_col_AT1G01010
    Chr1	10012	11763	01.col_col_AT1G01020
    Chr1	14961	16037	01.col_col_AT1G01030

    3.The sample.list file (the order of the samples determines the order of the geneindex column)
    Bna_Darmor_v10.0
    Bna_Darmor_v5.0
    Bna_Express617_v1.0

    4.The chro.list file
    Chr1
    Chr2
    Chr3


    """
    
    import argparse
    parser=argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i","--input",required=True,help="输入目录")
    parser.add_argument("-o","--output",required=True,help="输出目录")
    
    parser.add_argument("-s","--samplelist",type=str,default="",help="样本列表文件")
    parser.add_argument("-c","--chromlist",type=str,default="",help="染色体列表文件")
    parser.add_argument("-t","--thread",type=int,default=10,help="线程数")

    parser.add_argument("-b","--bed",type=str,choices=["bed4","bed6"],default="bed6",help="bed文件格式选择")
    args=parser.parse_args()

    input_dir=args.input
    output_dir=args.output

    sample_list_path=args.samplelist
    chromosome_list_path=args.chromlist
    thread=args.thread

    if args.bed=="bed4":
        config["bed_format"]["bed4"]=True
        config["bed_format"]["bed6"]=False
    elif args.bed=="bed6":
        config["bed_format"]["bed4"]=False
        config["bed_format"]["bed6"]=True

    
    ######## 检查 input_dir 是否存在 ########
    if not os.path.exists(input_dir):
        console.print(f"input directory {input_dir} does not exist",style="red")
        exit(1)
    
    ######## 检查 sample.list 是否存在 ########
    if not sample_list_path:
        sample_list_file=config["input_file_type"]["sample_file_name"]
        sample_list_path=os.path.join(input_dir,sample_list_file)
    if not os.path.exists(sample_list_path):
        console.print(f"sample list file {sample_list_path} does not exist",style="red")
        exit(1)
    else:
        console.print(f"sample list file {sample_list_path} exists",style="green")
    
    ######## 检查 chromosome.list 是否存在 ########
    if not chromosome_list_path:
        chromosome_list_file=config["input_file_type"]["chromosome_file_name"]
        chromosome_list_path=os.path.join(input_dir,chromosome_list_file)
    if not os.path.exists(chromosome_list_path):
        console.print(f"chromosome list file {chromosome_list_path} does not exist",style="red")
        exit(1)
    else:
        console.print(f"chromosome list file {chromosome_list_path} exists",style="green")
    
    ######## 读取 chromosome.list ########
    chromosome_list=[]
    with open(chromosome_list_path) as f:
        for line in f:
            line=line.strip()
            if not line:
                continue
            chromosome_list.append(line)

    ######## 读取 sample.list ########
    sample_list=[]
    with open(sample_list_path) as f:
        for line in f:
            line=line.strip()
            if not line:
                continue
            sample_list.append(line)
    
    ######## 检查 bed 文件是否存在 ########
    all_files_exist=[]
    for sample in sample_list:
        bed_file=os.path.join(input_dir,sample+config["input_file_type"]["bed_suffix"])
        if not os.path.exists(bed_file):
            all_files_exist.append(False)
        else:
            all_files_exist.append(True)
    for i,sample in enumerate(sample_list):
        bed_file=os.path.join(input_dir,sample+config["input_file_type"]["bed_suffix"])
        if not all_files_exist[i]:
            console.print(f"bed file {bed_file} does not exist",style="red")
        else:
            console.print(f"bed file {bed_file} exists",style="green")
    if not all(all_files_exist):
        console.print("Some bed files are missing, please check!",style="red")
        exit(1)

    ######## 检查 pep 文件是否存在 ########
    all_files_exist=[]
    for sample in sample_list:
        pep_file=os.path.join(input_dir,sample+config["input_file_type"]["pep_suffix"])
        if not os.path.exists(pep_file):
            all_files_exist.append(False)
        else:
            all_files_exist.append(True)
    for i,sample in enumerate(sample_list):
        pep_file=os.path.join(input_dir,sample+config["input_file_type"]["pep_suffix"])
        if not all_files_exist[i]:
            console.print(f"pep file {pep_file} does not exist",style="red")
        else:
            console.print(f"pep file {pep_file} exists",style="green")
    if not all(all_files_exist):
        console.print("Some pep files are missing, please check!",style="red")
        exit(1)
        
    ######## 检查 pep 和 bed 的输出文件夹是否存在 ########
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    working_dir=os.path.join(output_dir,config["work_config"]["working_dir_name"])
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    working_bed_dir=os.path.join(working_dir,config["work_config"]["chromosome_bed_dir_name"])
    working_pep_dir=os.path.join(working_dir,config["work_config"]["chromosome_pep_dir_name"])
    if not os.path.exists(working_bed_dir):
        os.makedirs(working_bed_dir)
    if not os.path.exists(working_pep_dir):
        os.makedirs(working_pep_dir)
    
    ######## 多线程运行分离 pep 和 bed ########
    console.print(f"Divide the BED and PEP files according to chromosomes...")
    with ThreadPoolExecutor(max_workers=thread) as executor:
        futures = {}
        for sample in sample_list:
            bed_file = os.path.join(input_dir, sample + config["input_file_type"]["bed_suffix"])
            pep_file = os.path.join(input_dir, sample + config["input_file_type"]["pep_suffix"])
            # 直接提交 split_pep_bed，并把 sample 记在映射中，用于日志
            fut = executor.submit(split_pep_bed, bed_file, pep_file,
                                working_pep_dir, working_bed_dir,
                                chromosome_list, sample)
            futures[fut] = sample

        for fut in as_completed(futures):
            sample = futures[fut]
            try:
                status = fut.result()   # split_pep_bed 应返回 "skip" 或 "done"
                if status == "skip":
                    console.print(f"[SKIP] {sample} already processed, skipped.", style="yellow")
                elif status == "done":
                    console.print(f"[OK]   {sample} split completed.", style="green")
                else:
                    console.print(f"[WARN] {sample} returned unknown status: {status}", style="yellow")
            except Exception as e:
                console.print(f"[ERROR] {sample} split failed: {e}", style="red")
                exit(1)
    
    ######## 检查 diamond 程序存在与否 ########
    diamond_path=config["diamond_path"]
    if not diamond_path:
        diamond_path=shutil.which("diamond")
    if not diamond_path:
        console.print(f"diamond not found in PATH",style="red")
        console.print(f"please install diamond and add it to PATH",style="red")
        exit(1)

    ######## 按照染色体进行批量运行 diamond makedb 和 blastp ########
    console.print("Run diamond makedb...")
    working_diamond_dir=os.path.join(working_dir,config["work_config"]["blastp_dir_name"])
    if not os.path.exists(working_diamond_dir):
        os.makedirs(working_diamond_dir)
    tasks_diamond_makedb=[]
    args_makedb = [str(item) for kv in config["diamond_makedb"].items() for item in kv]
    for chromosome in [*chromosome_list,"scaffold"]:
        for sample in sample_list:
            pep_path=os.path.join(working_pep_dir,sample+"_"+chromosome+".pep")
            diamond_db_path=os.path.join(working_diamond_dir,sample+"_"+chromosome+config["work_config"]["blastp_db_suffix"])
            tasks_diamond_makedb.append((diamond_path,pep_path,diamond_db_path,args_makedb))
    with ThreadPoolExecutor(max_workers=thread) as executor:
        futures = {executor.submit(run_diamond_makedb, *t): t for t in tasks_diamond_makedb}
        for fut in as_completed(futures):
            t = futures[fut]
            try:
                fut.result()
            except Exception as e:
                print(f"[ERROR] diamond makedb failed {t}: {e}")
                exit(1)
    ######## all vs all 运行 blastp 程序 ########
    console.print("Run diamond blastp...")
    tasks_diamond_blastp=[]
    for sample1 in sample_list: 
        for sample2 in sample_list:
            if sample1 == sample2:
                continue
            # chro vs chro 比对 -> 1 vs 1 、 2 vs 2 、 3 vs 3
            for chromosome in chromosome_list:
                qry = os.path.join(working_pep_dir,f"{sample1}_{chromosome}.pep")
                db = os.path.join(working_diamond_dir,f"{sample2}_{chromosome}{config['work_config']['blastp_db_suffix']}")
                out = os.path.join(working_diamond_dir,f"{sample1}_vs_{sample2}_{chromosome}{config['work_config']['blastp_suffix']}")
                tasks_diamond_blastp.append((qry, db, out))
            # scaffold vs chro 比对 
            # -> scaffold vs 1 、 scaffold vs 2 、 scaffold vs 3
            # -> 1 vs scaffold 、 2 vs scaffold 、 3 vs scaffold
            for chromosome in chromosome_list:
                qry=os.path.join(working_pep_dir,f"{sample1}_{chromosome}.pep") # A A01 vs B scaffold 、 B A01 vs A scaffold
                db=os.path.join(working_diamond_dir,f"{sample2}_scaffold{config['work_config']['blastp_db_suffix']}")
                out = os.path.join(working_diamond_dir,f"{sample1}_vs_{sample2}_{chromosome}_vs_scaffold{config['work_config']['blastp_suffix']}")
                tasks_diamond_blastp.append((qry, db, out))
                
                qry=os.path.join(working_pep_dir,f"{sample1}_scaffold.pep")# A scaffold vs B A01 、 B scaffold vs A A01
                db=os.path.join(working_diamond_dir,f"{sample2}_{chromosome}{config['work_config']['blastp_db_suffix']}")
                out = os.path.join(working_diamond_dir,f"{sample1}_vs_{sample2}_scaffold_vs_{chromosome}{config['work_config']['blastp_suffix']}")
                tasks_diamond_blastp.append((qry, db, out))

            # scaffold vs scaffold 比对 -> scaffold vs scaffold
            qry=os.path.join(working_pep_dir,f"{sample1}_scaffold.pep")# A scaffold vs B scaffold
            db=os.path.join(working_diamond_dir,f"{sample2}_scaffold{config['work_config']['blastp_db_suffix']}")
            out = os.path.join(working_diamond_dir,f"{sample1}_vs_{sample2}_scaffold{config['work_config']['blastp_suffix']}")
            tasks_diamond_blastp.append((qry, db, out))
            
    args_blastp = []
    for k, v in config["diamond_blastp"].items():
        args_blastp.append(k)
        if isinstance(v, list):
            args_blastp.extend(v)
        else:
            args_blastp.append(str(v))

    console.print(f"[INFO] Total blastp tasks: {len(tasks_diamond_blastp)}", style="cyan")
    console.print(f"[INFO] Using threads: {thread}", style="cyan")
    with ThreadPoolExecutor(max_workers=thread) as executor:
        futures = {
            executor.submit(run_diamond_blastp, diamond_path, qry, db, out, args_blastp): (qry, out)
            for qry, db, out in tasks_diamond_blastp
        }
        for future in as_completed(futures):
            qry, out = futures[future]
            try:
                status = future.result()
                if status == "skip":
                    console.print(f"[SKIP] {out} skipped, output exists.", style="yellow")
                elif status == "error":
                    console.print(f"[FAIL] {out} failed.", style="red")
                else:
                    console.print(f"[DONE] {out} completed.", style="green")
            except Exception as e:
                console.print(f"[CRASH] {out} exception: {e}", style="red")
                exit(1)

    ######## 获取 run_DAG_chainer.pl 输入文件 ########
    working_dagchainer_dir=os.path.join(working_dir,config["work_config"]["dagchainer_input_dir_name"])
    if not os.path.exists(working_dagchainer_dir):
        os.makedirs(working_dagchainer_dir)
    tasks_dagchainer_input=[]
    for sample1 in sample_list: 
        for sample2 in sample_list:
            if sample1==sample2:
                continue
            # chro vs chro
            for chromosome in chromosome_list:
                blastp_out1=os.path.join(working_diamond_dir,f"{sample1}_vs_{sample2}_{chromosome}{config['work_config']['blastp_suffix']}")
                blastp_out2=os.path.join(working_diamond_dir,f"{sample2}_vs_{sample1}_{chromosome}{config['work_config']['blastp_suffix']}")
                bed_file1=os.path.join(working_bed_dir,f"{sample1}_{chromosome}.bed")
                bed_file2=os.path.join(working_bed_dir,f"{sample2}_{chromosome}.bed")
                dagchainer_input1=os.path.join(working_dagchainer_dir,f"{sample1}_vs_{sample2}_{chromosome}{config['work_config']['dagchainer_input_suffix']}")
                dagchainer_input2=os.path.join(working_dagchainer_dir,f"{sample2}_vs_{sample1}_{chromosome}{config['work_config']['dagchainer_input_suffix']}")
                tasks_dagchainer_input.append((blastp_out1,blastp_out2,sample1, sample2,bed_file1,bed_file2,dagchainer_input1,dagchainer_input2))
            # scaffold vs chro
            for chromosome in chromosome_list:
                blastp_out1=os.path.join(working_diamond_dir,f"{sample1}_vs_{sample2}_scaffold_vs_{chromosome}{config['work_config']['blastp_suffix']}")
                blastp_out2=os.path.join(working_diamond_dir,f"{sample2}_vs_{sample1}_{chromosome}_vs_scaffold{config['work_config']['blastp_suffix']}")
                bed_file1=os.path.join(working_bed_dir,f"{sample1}_scaffold.bed")
                bed_file2=os.path.join(working_bed_dir,f"{sample2}_{chromosome}.bed")
                dagchainer_input1=os.path.join(working_dagchainer_dir,f"{sample1}_vs_{sample2}_scaffold_vs_{chromosome}{config['work_config']['dagchainer_input_suffix']}")
                dagchainer_input2=os.path.join(working_dagchainer_dir,f"{sample2}_vs_{sample1}_{chromosome}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}")
                tasks_dagchainer_input.append((blastp_out1,blastp_out2,sample1, sample2,bed_file1,bed_file2,dagchainer_input1,dagchainer_input2))
            # chro vs scaffold
            for chromosome in chromosome_list:
                blastp_out1=os.path.join(working_diamond_dir,f"{sample1}_vs_{sample2}_{chromosome}_vs_scaffold{config['work_config']['blastp_suffix']}")
                blastp_out2=os.path.join(working_diamond_dir,f"{sample2}_vs_{sample1}_scaffold_vs_{chromosome}{config['work_config']['blastp_suffix']}")
                bed_file1=os.path.join(working_bed_dir,f"{sample1}_{chromosome}.bed")
                bed_file2=os.path.join(working_bed_dir,f"{sample2}_scaffold.bed")
                dagchainer_input1=os.path.join(working_dagchainer_dir,f"{sample1}_vs_{sample2}_{chromosome}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}")
                dagchainer_input2=os.path.join(working_dagchainer_dir,f"{sample2}_vs_{sample1}_scaffold_vs_{chromosome}{config['work_config']['dagchainer_input_suffix']}")
                tasks_dagchainer_input.append((blastp_out1,blastp_out2,sample1, sample2,bed_file1,bed_file2,dagchainer_input1,dagchainer_input2))
            # scaffold vs scaffold
            blastp_out1=os.path.join(working_diamond_dir,f"{sample1}_vs_{sample2}_scaffold{config['work_config']['blastp_suffix']}")
            blastp_out2=os.path.join(working_diamond_dir,f"{sample2}_vs_{sample1}_scaffold{config['work_config']['blastp_suffix']}")
            bed_file1=os.path.join(working_bed_dir,f"{sample1}_scaffold.bed")
            bed_file2=os.path.join(working_bed_dir,f"{sample2}_scaffold.bed")
            dagchainer_input1=os.path.join(working_dagchainer_dir,f"{sample1}_vs_{sample2}_scaffold{config['work_config']['dagchainer_input_suffix']}")
            dagchainer_input2=os.path.join(working_dagchainer_dir,f"{sample2}_vs_{sample1}_scaffold{config['work_config']['dagchainer_input_suffix']}")
            tasks_dagchainer_input.append((blastp_out1,blastp_out2,sample1, sample2,bed_file1,bed_file2,dagchainer_input1,dagchainer_input2))
    with ThreadPoolExecutor(max_workers=thread) as executor:
        futures = {
            executor.submit(get_dagchainer_input, b1, b2, sample1, sample2, bed1, bed2, o1, o2): 
            (b1, b2, sample1, sample2, bed1, bed2, o1, o2)
            for (b1, b2, sample1, sample2, bed1, bed2, o1, o2) in tasks_dagchainer_input
        }
        for future in as_completed(futures):
            b1, b2, sample1, sample2, bed1, bed2, o1, o2 = futures[future]
            try:
                status = future.result()
                if status == "skip":
                    console.print(f"[SKIP] run_DAG_chainer.pl input file generation skipped: {o1}", style="yellow")
                    console.print(f"[SKIP] run_DAG_chainer.pl input file generation skipped: {o2}", style="yellow")
                else:
                    console.print(f"[OK] run_DAG_chainer.pl input file generation completed: {o1}", style="green")
                    console.print(f"[OK] run_DAG_chainer.pl input file generation completed: {o2}", style="green")
            except Exception as e:
                console.print(f"[ERROR] run_DAG_chainer.pl input file generation failed: {o2} : {e}", style="red")
                console.print(f"[ERROR] run_DAG_chainer.pl input file generation failed: {o2} : {e}", style="red")
                exit(1)
    
    ######## 运行 run_DAG_chainer.pl 和 perl ########
    perl_path = config["perl_path"]
    if not perl_path:
        perl_path = shutil.which("perl")
    if not perl_path:
        console.print(f"perl not found in PATH", style="red")
        console.print(f"please install perl and add it to PATH", style="red")
        exit(1)

    rundagchainer_pl_path = config["run_DAG_chainer_path"]
    if not rundagchainer_pl_path:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        rundagchainer_pl_path = os.path.join(script_dir, "run_DAG_chainer.pl")

    if not rundagchainer_pl_path:
        console.print(f"run_DAG_chainer.pl not found", style="red")
        exit(1)

    working_dagchainer_dir = os.path.join(working_dir, config["work_config"]["dagchainer_input_dir_name"])
    working_dagchainer_output_dir = os.path.join(working_dir, config["work_config"]["dagchainer_output_dir_name"])
    if not os.path.exists(working_dagchainer_output_dir):
        os.makedirs(working_dagchainer_output_dir)

    rundagchainer_pl_args = [str(item) for kv in config["DAGchainer_config"].items() for item in kv]

    tasks_dagchainer_run = []

    for sample1 in sample_list:
        for sample2 in sample_list:
            if sample1 == sample2:
                continue

            # --------------------------
            # 1. chr vs chr
            # --------------------------
            for chromosome in chromosome_list:
                inp = os.path.join(
                    working_dagchainer_dir,
                    f"{sample1}_vs_{sample2}_{chromosome}{config['work_config']['dagchainer_input_suffix']}"
                )
                out = f"{inp}{config['work_config']['dagchainer_output_suffix']}"
                tasks_dagchainer_run.append((inp, out))

            # --------------------------
            # 2. scaffold vs chr
            # --------------------------
            for chromosome in chromosome_list:
                inp = os.path.join(
                    working_dagchainer_dir,
                    f"{sample1}_vs_{sample2}_scaffold_vs_{chromosome}{config['work_config']['dagchainer_input_suffix']}"
                )
                out = f"{inp}{config['work_config']['dagchainer_output_suffix']}"
                tasks_dagchainer_run.append((inp, out))

            # --------------------------
            # 3. chr vs scaffold
            # --------------------------
            for chromosome in chromosome_list:
                inp = os.path.join(
                    working_dagchainer_dir,
                    f"{sample1}_vs_{sample2}_{chromosome}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}"
                )
                out = f"{inp}{config['work_config']['dagchainer_output_suffix']}"
                tasks_dagchainer_run.append((inp, out))

            # --------------------------
            # 4. scaffold vs scaffold
            # --------------------------
            inp = os.path.join(
                working_dagchainer_dir,
                f"{sample1}_vs_{sample2}_scaffold{config['work_config']['dagchainer_input_suffix']}"
            )
            out = f"{inp}{config['work_config']['dagchainer_output_suffix']}"
            tasks_dagchainer_run.append((inp, out))


    # --------------------------
    # 并行运行 DAGchainer
    # --------------------------
    with ThreadPoolExecutor(max_workers=thread) as executor:
        futures = {
            executor.submit(
                run_DAG_chainer_pl,
                perl_path, rundagchainer_pl_path, inp, out, rundagchainer_pl_args
            ): (inp, out)
            for inp, out in tasks_dagchainer_run
        }

        for future in as_completed(futures):
            inp, out = futures[future]
            try:
                status = future.result()
                if status == "skip":
                    console.print(f"[SKIP] {out}", style="yellow")
                elif status == "error":
                    console.print(f"[FAIL] {out}", style="red")
                else:
                    console.print(f"[DONE] {out}", style="green")
            except Exception as e:
                console.print(f"[CRASH] {out} exception: {e}", style="red")
                exit(1)

    ######## 过滤 重复基因 ########
    result_SG=[]
    for sample1,sample2 in combinations(sample_list,2):
        result_SG.append(sample1_vs_sample2_pair(sample1,sample2,
                                bed_dir=working_bed_dir,
                                chromosome_list=chromosome_list,
                                dagchainer_output_dir=working_dagchainer_output_dir,
                                )
                        )

    ######## 合并 SG 形成 geneindex ########
    '''results_SG: [ (sample1,sample2,chromosome,[(gene1,gene2),...]), ... ] '''
    console.print(f"combine SG into geneindex", style="green")
    edges_list=combine_SG(result_SG,chromosome_list,sample_list,bed_dir=working_bed_dir)
    result_SG.clear()

    ######## 均一化 blastp 的得分为 RBNH : 获取 mcl 输入文件 ########
    RBNH_output_dir_name=os.path.join(working_dir,config["work_config"]["RBNH_output_dir_name"])
    if not os.path.exists(RBNH_output_dir_name):
        os.makedirs(RBNH_output_dir_name)
    compute_RBNH(sample_list,working_diamond_dir,RBNH_output_dir_name,chromosome_list,thread)

    ######## 使用 mcl 的聚类结果对 geneindex 行进行拆分 ########
    mcl_path=config["mcl_path"]
    if not mcl_path:
        mcl_path=shutil.which("mcl")
    if not mcl_path:
        console.print(f"mcl not found in PATH",style="red")
        console.print(f"please install mcl and add it to PATH",style="red")
        exit(1)

    all_RBNH_file_path=os.path.join(RBNH_output_dir_name,config["work_config"]["all_RBNH_file_name"])
    mcl_input_dir=os.path.join(working_dir,config["work_config"]["mcl_input_dir_name"])
    console.print(f"get duplicate RBNH parallel", style="green")
    duplicate_gene_pair_file_list=get_duplicate_RBNH_parallel(all_RBNH_file_path,edges_list,mcl_input_dir,thread)
    
    run_mcl_task_list=[]
    mcl_output_dir=os.path.join(working_dir,config["work_config"]["mcl_output_dir_name"])
    if not os.path.exists(mcl_output_dir):
        os.makedirs(mcl_output_dir)
    mcl_args_list=[str(item) for kv in config["mcl_config"].items() for item in kv ]
    for mcl_input_file_path in duplicate_gene_pair_file_list:
        mcl_input_file_basename=os.path.basename(mcl_input_file_path)
        mcl_input_file_basename_no_suffix=mcl_input_file_basename.split(config["work_config"]["mcl_input_suffix"])[0]
        out_file_path=os.path.join(mcl_output_dir,f"{mcl_input_file_basename_no_suffix}{config['work_config']['mcl_output_suffix']}")
        run_mcl_task_list.append((mcl_path,mcl_input_file_path,out_file_path,mcl_args_list))
    with ThreadPoolExecutor(max_workers=thread) as executor:
        future_map = {executor.submit(run_mcl, *task): task for task in run_mcl_task_list}
        for future in as_completed(future_map):
            mcl_path, inp, out, _ = future_map[future]
            try:
                res = future.result()
                console.print(f"[MCL] {os.path.basename(inp)} -> {res}", style="yellow")
            except Exception as e:
                console.print(f"[ERROR] {inp}: {e}", style="red")
    
    ######## 归并 mcl 的多个结果为一个文件 ########
    mcl_all_result_path=os.path.join(mcl_output_dir,config["work_config"]["mcl_output_merge_file_name"])
    if config["overwrite"]["merge_mcl_output"] or not os.path.exists(mcl_all_result_path):
        with open(mcl_all_result_path,"w") as f:
            for this_task in run_mcl_task_list:
                this_mcl_output_path=this_task[2]
                with open(this_mcl_output_path) as f2:
                    for line in f2:
                        line=line.strip()
                        if line:
                            f.write(f"{line}\n")

    working_split_mcl_output_dir=os.path.join(output_dir,config["work_config"]["mcl_split_geneindex_dir_name"])
    if not os.path.exists(working_split_mcl_output_dir):
        os.makedirs(working_split_mcl_output_dir)
    
    gene2cluster_dict = read_mcl_cluster(mcl_all_result_path)
    geneindex_in_path = os.path.join(output_dir,config['work_config']['geneindex_output_file_name'])
    geneindex_out_path = os.path.join(output_dir,config['work_config']['geneindex_output_file_final_name'])
    run_mcl_split(gene2cluster_dict, geneindex_in_path, geneindex_out_path)
    console.print(f"[OK] all finished, results are in {output_dir}", style="green")

