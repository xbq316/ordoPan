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
from rich.logging import RichHandler
from rich.progress import Progress
import logging
from intervaltree import IntervalTree
import pandas as pd
import numpy as np
import math
from typing import Dict, List, Tuple, Any
from Bio import SeqIO
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from threading import Lock
import shutil
from itertools import combinations
import sys
import tempfile
import networkx as nx
from multiprocessing import cpu_count
import argparse
from dataclasses import dataclass, field
from datetime import datetime
import traceback

# 配置日志系统
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)]
)
logger = logging.getLogger("synpan")

# 全局配置
config = {
    "diamond_blastp": {
        "-p": 1,
        "-e": 1e-20,
        "-f": ["6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"],
        "--id": 40,
        "-k": 3
    },
    "diamond_makedb": {"--threads": 1},
    "diamond_path": "",
    "run_DAG_chainer_path": "",
    "perl_path": "",
    "mcl_path": "",
    "work_config": {
        "working_dir_name": "working",
        "blastp_dir_name": "blastp",
        "blastp_db_suffix": ".dmnd",
        "blastp_suffix": ".blastp",
        "chromosome_pep_dir_name": "pep",
        "chromosome_bed_dir_name": "bed",
        "dagchainer_input_dir_name": "dagchainer",
        "dagchainer_output_dir_name": "dagchainer",
        "dagchainer_input_suffix": ".dagchainer",
        "dagchainer_output_suffix": ".aligncoords",
        "geneindex_output_file_name": "geneindex.txt",
        "geneindex_without_sample_name": "geneindex.final.txt",
        "geneindex_output_file_split_name": "geneindex.split.txt",
        "geneindex_output_file_split_name_without_sample": "geneindex.split.final.txt",
        "RBNH_output_dir_name": "RBNH",
        "RBNH_suffix": ".RBNH",
        "all_RBNH_file_name": "all_sample_chromosome.RBNH",
        "mcl_output_dir_name": "mcl",
        "mcl_output_merge_file_name": "all.cluster",
        "mcl_input_dir_name": "mcl",
        "mcl_input_suffix": ".mcl",
        "mcl_output_suffix": ".cluster",
    },
    "gene_filter_method": {
        "methods": ["filter_synteny_blocks", "filter_synteny_blocks_normal"],
        "default": "filter_synteny_blocks",
    },
    "bed_format": {"bed4": False, "bed6": True},
    "overwrite": {
        "pep_bed_split": False,
        "makedb": False,
        "blastp": False,
        "mcl": False,
        "dagchainer_input": False,
        "dagchainer_output": False,
        "geneindex_output": True,
        "RBNH_compute": False,
        "mcl_output": False,
        "merge_mcl_output": False
    },
    "mcl_config": {"-I": 1.5, "--abc": ""},
    "input_file_type": {
        "pep_suffix": ".pep",
        "bed_suffix": ".bed",
        "chromosome_file_name": "chro.list",
        "sample_file_name": "sample.list"
    },
    "DAGchainer_config": {"-D": 10, "-g": 1, "-A": 5},
    "synteny_config": {"fotmat": "dagchainer"},
}

# 配置验证
if sum([config["bed_format"]["bed4"], config["bed_format"]["bed6"]]) != 1:
    raise ValueError("Config error: choose exactly one type: bed4 or bed6.")

# 全局锁和控制台
print_lock = Lock()
db_lock = Lock()
console = Console(soft_wrap=False)


class PipelineLogger:
    """统一的日志管理类"""
    
    @staticmethod
    def start_pipeline():
        """开始运行提示"""
        logger.info("🚀 Starting SynPan Pipeline")
        logger.info(f"📅 Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("=" * 60)
    
    @staticmethod
    def end_pipeline():
        """结束运行提示"""
        logger.info("=" * 60)
        logger.info("✅ SynPan Pipeline Completed Successfully")
        logger.info(f"📅 End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    @staticmethod
    def start_stage(stage_name: str):
        """开始一个处理阶段"""
        logger.info(f"🔄 Starting stage: {stage_name}")
    
    @staticmethod
    def end_stage(stage_name: str):
        """结束一个处理阶段"""
        logger.info(f"✅ Completed stage: {stage_name}")
    
    @staticmethod
    def info(message: str):
        """信息日志"""
        logger.info(f"ℹ️  {message}")
    
    @staticmethod
    def warning(message: str):
        """警告日志"""
        logger.warning(f"⚠️  {message}")
    
    @staticmethod
    def error(message: str):
        """错误日志"""
        logger.error(f"❌ {message}")
    
    @staticmethod
    def success(message: str):
        """成功日志"""
        logger.info(f"✅ {message}")
    
    @staticmethod
    def skip(message: str):
        """跳过日志"""
        logger.info(f"⏭️  {message}")


class FileProcessor:
    """文件处理模块"""
    
    @staticmethod
    def validate_input_files(input_dir: str, sample_list: List[str], 
                           file_suffix: str, file_type: str) -> bool:
        """验证输入文件是否存在"""
        missing_files = []
        for sample in sample_list:
            file_path = os.path.join(input_dir, f"{sample}{file_suffix}")
            if not os.path.exists(file_path):
                missing_files.append(file_path)
        
        if missing_files:
            for file_path in missing_files:
                PipelineLogger.error(f"{file_type} file missing: {file_path}")
            return False
        
        for sample in sample_list:
            file_path = os.path.join(input_dir, f"{sample}{file_suffix}")
            PipelineLogger.success(f"{file_type} file exists: {file_path}")
        
        return True
    
    @staticmethod
    def ensure_directory(directory: str) -> str:
        """确保目录存在，返回目录路径"""
        os.makedirs(directory, exist_ok=True)
        return directory
    
    @staticmethod
    def read_list_file(file_path: str, file_type: str) -> List[str]:
        """读取列表文件"""
        if not os.path.exists(file_path):
            PipelineLogger.error(f"{file_type} file does not exist: {file_path}")
            raise FileNotFoundError(f"{file_type} file not found: {file_path}")
        
        items = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    items.append(line)
        
        PipelineLogger.success(f"{file_type} file loaded: {file_path} ({len(items)} items)")
        return items


class ProcessManager:
    """进程/线程管理类"""
    
    @staticmethod
    def get_optimal_workers(task_count: int, max_workers: int = None) -> int:
        """计算最优的工作线程/进程数"""
        if max_workers is None:
            max_workers = min(cpu_count(), task_count)
        return min(max_workers, task_count)
    
    @staticmethod
    def run_parallel_tasks(tasks: List[Tuple], task_function: callable, 
                          max_workers: int = None, use_processes: bool = False,
                          task_description: str = "Processing", 
                          continue_on_error: bool = True) -> Dict[str, Any]:
        """
        并行执行任务，返回统计信息
        """
        total_tasks = len(tasks)
        if total_tasks == 0:
            return {"success": 0, "skipped": 0, "errors": 0, "results": []}
            
        workers = ProcessManager.get_optimal_workers(total_tasks, max_workers)
        executor_class = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
        
        PipelineLogger.info(f"{task_description}: {total_tasks} tasks, using {workers} {'processes' if use_processes else 'threads'}")
        
        results = []
        success_count = 0
        error_count = 0
        skip_count = 0
        
        with Progress() as progress:
            task_id = progress.add_task(f"{task_description}...", total=total_tasks)
            
            with executor_class(max_workers=workers) as executor:
                # 提交所有任务
                future_to_task = {
                    executor.submit(task_function, *task_args): task_args 
                    for task_args in tasks
                }
                
                # 处理完成的任务
                for future in as_completed(future_to_task):
                    task_args = future_to_task[future]
                    try:
                        result = future.result()
                        results.append(result)
                        
                        if result == "skip":
                            skip_count += 1
                        elif result == "error":
                            error_count += 1
                        else:
                            success_count += 1
                            
                    except Exception as e:
                        PipelineLogger.error(f"Task failed with exception: {task_args} -> {e}")
                        error_count += 1
                        if continue_on_error:
                            results.append(None)
                        else:
                            # 如果不允许继续，重新抛出异常
                            raise
                    
                    progress.update(task_id, advance=1)
        
        # 统计信息
        PipelineLogger.info(f"Task completion: {success_count} success, {skip_count} skipped, {error_count} errors")
        
        return {
            "success": success_count,
            "skipped": skip_count, 
            "errors": error_count,
            "results": results
        }

    @staticmethod
    def is_stage_successful(stats: Dict[str, Any]) -> bool:
        """判断阶段是否成功：只要有成功任务或只有跳过，就认为是成功的"""
        return stats["success"] > 0 or (stats["errors"] == 0 and stats["skipped"] > 0)


class DiamondManager:
    """Diamond 比对管理模块"""
    
    @staticmethod
    def get_diamond_path() -> str:
        """获取diamond路径"""
        diamond_path = config["diamond_path"]
        if not diamond_path:
            diamond_path = shutil.which("diamond")
        if not diamond_path:
            PipelineLogger.warning("diamond not found in PATH")
            return None
        return diamond_path
    
    @staticmethod
    def run_diamond_makedb(diamond_path: str, ref_pep: str, db_name: str, args: List[str]) -> str:
        """运行 diamond makedb"""
        with db_lock:
            if not config["overwrite"]["makedb"] and os.path.exists(db_name) and os.path.getsize(db_name) > 0:
                PipelineLogger.skip(f"DB exists: {db_name}")
                return "skip"
            
            # 检查输入文件
            if not os.path.exists(ref_pep) or os.path.getsize(ref_pep) == 0:
                PipelineLogger.warning(f"PEP file not found or empty: {ref_pep}")
                return "skip"

            cmd_list = [diamond_path, "makedb", "--in", ref_pep, "-d", db_name, *args]
            PipelineLogger.info(f"Creating database: {' '.join(cmd_list)}")
        
        result = subprocess.run(cmd_list, capture_output=True, text=True)
        if result.returncode == 0:
            PipelineLogger.success(f"Created DB: {db_name}")
            return "done"
        else:
            PipelineLogger.error(f"Failed to create {db_name}: {result.stderr}")
            return "error"
    
    @staticmethod
    def run_diamond_blastp(diamond_path: str, qry_pep: str, db_name: str,
                          output: str, args: List[str]) -> str:
        """运行 diamond blastp"""
        if not config["overwrite"]["blastp"] and os.path.exists(output) and os.path.getsize(output) > 0:
            PipelineLogger.skip(f"BLASTP output exists: {output}")
            return "skip"
        
        # 检查输入文件是否存在
        if not os.path.exists(qry_pep) or os.path.getsize(qry_pep) == 0:
            PipelineLogger.warning(f"Query PEP file not found or empty: {qry_pep}")
            return "skip"
            
        if not os.path.exists(db_name) or os.path.getsize(db_name) == 0:
            PipelineLogger.warning(f"Database file not found or empty: {db_name}")
            return "skip"
        
        cmd_list = [diamond_path, "blastp", "-q", qry_pep, "-d", db_name, *args, "-o", output]
        PipelineLogger.info(f"Running BLASTP: {os.path.basename(qry_pep)} vs {os.path.basename(db_name)}")
        
        result = subprocess.run(cmd_list, capture_output=True, text=True)
        if result.returncode == 0:
            PipelineLogger.success(f"BLASTP completed: {output}")
            return "done"
        else:
            PipelineLogger.error(f"BLASTP failed: {output} -> {result.stderr}")
            return "error"


class DAGchainerManager:
    """DAGchainer 管理模块"""
    
    @staticmethod
    def get_perl_path() -> str:
        """获取perl路径"""
        perl_path = config["perl_path"]
        if not perl_path:
            perl_path = shutil.which("perl")
        if not perl_path:
            PipelineLogger.warning("perl not found in PATH")
            return None
        return perl_path
    
    @staticmethod
    def get_dagchainer_path() -> str:
        """获取DAGchainer路径"""
        dagchainer_path = config["run_DAG_chainer_path"]
        if not dagchainer_path:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            dagchainer_path = os.path.join(script_dir, "run_DAG_chainer.pl")
        
        if not dagchainer_path or not os.path.exists(dagchainer_path):
            PipelineLogger.warning(f"DAGchainer script not found: {dagchainer_path}")
            return None
        
        return dagchainer_path
    
    @staticmethod
    def run_DAG_chainer_pl(perl_path: str, dagchainer_path: str,
                          input_dagchainer: str, output: str,
                          args: List[str]) -> str:
        """运行DAGchainer"""
        if not config["overwrite"]["dagchainer_output"] and os.path.exists(output) and os.path.getsize(output) > 0:
            PipelineLogger.skip(f"DAGchainer output exists: {output}")
            return "skip"

        # 检查输入文件
        if not os.path.exists(input_dagchainer) or os.path.getsize(input_dagchainer) == 0:
            PipelineLogger.warning(f"DAGchainer input file not found or empty: {input_dagchainer}")
            return "skip"

        PipelineLogger.info(f"Running DAGchainer: {os.path.basename(input_dagchainer)}")
        
        # 创建临时目录执行
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_input = os.path.join(tmpdir, os.path.basename(input_dagchainer))
            shutil.copy(input_dagchainer, tmp_input)

            cmd_list = [perl_path, dagchainer_path, "-i", tmp_input, *args]
            
            result = subprocess.run(cmd_list, capture_output=True, text=True, cwd=tmpdir)

            if result.returncode != 0:
                PipelineLogger.error(f"DAGchainer failed: {input_dagchainer}")
                PipelineLogger.error(f"Error: {result.stderr}")
                return "error"

            # 移动输出文件
            tmp_output_file = os.path.join(tmpdir, os.path.basename(output))
            if os.path.exists(tmp_output_file):
                shutil.move(tmp_output_file, output)
            else:
                PipelineLogger.error(f"DAGchainer did not produce expected output: {output}")
                return "error"

        PipelineLogger.success(f"DAGchainer completed: {output}")
        return "done"


class MCLManager:
    """MCL 聚类管理模块"""
    
    @staticmethod
    def get_mcl_path() -> str:
        """获取mcl路径"""
        mcl_path = config["mcl_path"]
        if not mcl_path:
            mcl_path = shutil.which("mcl")
        if not mcl_path:
            PipelineLogger.warning("mcl not found in PATH")
            return None
        return mcl_path
    
    @staticmethod
    def run_mcl(mcl_path: str, mcl_input_path: str, mcl_output_path: str, args: List[str]) -> str:
        """运行MCL聚类"""
        if (not config["overwrite"]["mcl_output"]) and os.path.exists(mcl_output_path) and os.path.getsize(mcl_output_path) > 0:
            PipelineLogger.skip(f"MCL output exists: {mcl_output_path}")
            return "skip"
        
        # 检查输入文件
        if not os.path.exists(mcl_input_path) or os.path.getsize(mcl_input_path) == 0:
            PipelineLogger.warning(f"MCL input file not found or empty: {mcl_input_path}")
            return "skip"
        
        mcl_args = " ".join(str(x) for x in args)
        cmd = f"{mcl_path} {mcl_input_path} {mcl_args} -o {mcl_output_path}"
        
        PipelineLogger.info(f"Running MCL: {os.path.basename(mcl_input_path)}")
        
        try:
            result = subprocess.run(
                cmd, shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )
            
            if result.returncode == 0:
                PipelineLogger.success(f"MCL completed: {mcl_output_path}")
                return "done"
            else:
                PipelineLogger.error(f"MCL failed: {mcl_output_path} -> {result.stdout}")
                return "error"
        except Exception as e:
            PipelineLogger.error(f"MCL execution failed: {e}")
            return "error"


# 原有的数据类定义保持不变
class Tree(IntervalTree):
    """闭区间的区间树"""
    def __init__(self, intervals=None):
        super().__init__(intervals)

    def add_closed(self, start: int, end: int):
        if isinstance(start, int) and isinstance(end, int):
            pass
        else:
            raise ValueError("start and end must be int")            
        if start > end:
            raise ValueError("start must <= end")

        self.addi(start, end + 1)
        self.merge_overlaps()

    def remove_closed(self, start: int, end: int):
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

    def get_uncovered(self, start: int, end: int) -> List[Tuple[int, int]]:
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
        val = (self.qlen or 0) * (self.slen or 0)
        if val is None or not np.isfinite(val) or val <= 0:
            val = 1e-8
        self.L_qh = val

    def compute_Bprime(self, a: float, b: float):
        L = self.L_qh if (self.L_qh is not None and np.isfinite(self.L_qh) and self.L_qh > 0) else 1e-8
        bits = self.bits if (self.bits is not None and np.isfinite(self.bits) and self.bits > 0) else 1e-8
        self.Bprime = bits / (10 ** b * (L ** a))

    def __str__(self):
        return f"{self.q}\t{self.s}\t{self.Bprime}"


@dataclass
class Block:
    chromosome1: str = ""
    genelist1: List[str] = field(default_factory=list)
    positions1: List[int] = field(default_factory=list)

    chromosome2: str = ""
    genelist2: List[str] = field(default_factory=list)
    positions2: List[int] = field(default_factory=list)

    hits: int = 0
    size1: int = 0
    size2: int = 0
    reverse1: bool = False
    reverse2: bool = False

    def __post_init__(self):
        a = [len(self.genelist1), len(self.genelist2),
             len(self.positions1), len(self.positions2)]
        assert all(i == a[0] for i in a)

        self.hits = a[0]

        if self.hits == 0:
            return

        self.size1 = abs(self.positions1[-1] - self.positions1[0]) + 1
        self.size2 = abs(self.positions2[-1] - self.positions2[0]) + 1
        self.reverse1 = self.positions1[0] > self.positions1[-1]
        self.reverse2 = self.positions2[0] > self.positions2[-1]

    def __str__(self):
        return f"qry-chr:{self.chromosome1},qry-size:{self.size1},ref-chr:{self.chromosome2},ref-size:{self.size2},hits:{self.hits}"

    def __repr__(self):
        return str(self)

    def intersection(self, left_segments: List[Tuple[int,int]], right_segments: List[Tuple[int,int]]) -> "Block":
        left_segments_len = len(left_segments)
        right_segments_len = len(right_segments)
        if left_segments_len == 0 or right_segments_len == 0:
            return Block(self.chromosome1,[],[],self.chromosome2,[],[])
        
        left_geneorder2index = {position:index for index,position in enumerate(self.positions1)}
        right_geneorder2index = {position:index for index,position in enumerate(self.positions2)}
        
        left_index_list = []
        right_index_list = []
        for start,end in left_segments:
            start_index = left_geneorder2index[start]
            end_index = left_geneorder2index[end]
            left_index_list.extend(list(range(start_index, end_index+1)))
        for start,end in right_segments:
            start_index = right_geneorder2index[start]
            end_index = right_geneorder2index[end]
            right_index_list.extend(list(range(start_index, end_index+1)))
        
        left_index_set = set(left_index_list)
        right_index_set = set(right_index_list)
        final_index_list = sorted(left_index_set.intersection(right_index_set))
        
        new_genelist1 = [self.genelist1[i] for i in final_index_list]
        new_positions1 = [self.positions1[i] for i in final_index_list]
        new_genelist2 = [self.genelist2[i] for i in final_index_list]
        new_positions2 = [self.positions2[i] for i in final_index_list]
        
        return Block(self.chromosome1, new_genelist1, new_positions1, self.chromosome2, new_genelist2, new_positions2)
    
    def adjust_segments(self, left_segments: List[Tuple[int,int]], right_segments: List[Tuple[int,int]]) -> Tuple[List[Tuple[int,int]], List[Tuple[int,int]]]:
        new_left_segments = []
        for start,end in left_segments:
            a = []
            for i in self.positions1:
                if i >= start and i <= end:
                    a.append(i)
            if len(a) == 0:
                continue
            start = min(a)
            end = max(a)
            new_left_segments.append((start, end))
        
        new_right_segments = []
        for start,end in right_segments:
            a = []
            for i in self.positions2:
                if i >= start and i <= end:
                    a.append(i)
            if len(a) == 0:
                continue
            start = min(a)
            end = max(a)
            new_right_segments.append((start, end))
        
        return new_left_segments, new_right_segments
    
    def is_normal(self) -> bool:
        a = [len(self.genelist1), len(self.genelist2), len(self.positions1), len(self.positions2)]
        return all(i == a[0] for i in a)

    def length(self) -> int:
        if not self.is_normal():
            raise ValueError("Block is not normal")
        return len(self.genelist1)


class Gene:
    def __init__(self, chromosome: str, start: int, end: int, name: str, score: str = ".", strand: str = "+", index=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.index = index

    def __str__(self):
        return f"{self.chromosome},{self.start},{self.end},{self.name},{self.score},{self.strand}"

    def __repr__(self):
        return str(self)


# 核心功能函数
def read_blast_file(path: str) -> Dict[frozenset[str], blast]:
    """读取BLAST文件"""
    result: Dict[frozenset[str], blast] = {}
    with open(path, 'r') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12:
                continue
            try:
                q = parts[0]
                s = parts[1]
                bits = float(parts[11]) if parts[11] != '' else float('nan')
                evalue = float(parts[10]) if parts[10] != '' else float('nan')
                if len(parts) >= 14:
                    qlen = int(float(parts[12])) if parts[12] != '' else 0
                    slen = int(float(parts[13])) if parts[13] != '' else 0
                else:
                    qlen = int(float(parts[3])) if parts[3] != '' else 0
                    slen = int(float(parts[3])) if parts[3] != '' else 0
            except Exception:
                continue

            key = frozenset([q, s])
            obj = blast(q, s, bits, qlen, slen, evalue)
            existing = result.get(key)
            if existing is None or (np.isfinite(obj.bits) and np.isfinite(existing.bits) and obj.bits > existing.bits):
                result[key] = obj
            elif existing is not None and (not np.isfinite(existing.bits)) and np.isfinite(obj.bits):
                result[key] = obj
    return result


def merge_max_bits(blast_dict1: Dict[frozenset[str], blast], blast_dict2: Dict[frozenset[str], blast]) -> List[blast]:
    """合并两个方向的比对结果"""
    merge_result: List[blast] = []
    keys1 = set(blast_dict1.keys())
    keys2 = set(blast_dict2.keys())
    inter = keys1 & keys2
    for k in inter:
        b1 = blast_dict1[k]
        b2 = blast_dict2[k]
        b1_bits = b1.bits if np.isfinite(b1.bits) else -np.inf
        b2_bits = b2.bits if np.isfinite(b2.bits) else -np.inf
        merge_result.append(b1 if b1_bits >= b2_bits else b2)
    for k in keys1 - keys2:
        merge_result.append(blast_dict1[k])
    for k in keys2 - keys1:
        merge_result.append(blast_dict2[k])
    return merge_result


def compute_L_qh(blastp_list: List[blast]) -> None:
    for obj in blastp_list:
        obj.compute_L_qh()


def bin_top_hits(blastp_list: List[blast], bin_size=1000, top_pct=0.05):
    """分箱选择top hits"""
    rows = []
    for obj in blastp_list:
        L = obj.L_qh if (obj.L_qh is not None and np.isfinite(obj.L_qh) and obj.L_qh > 0) else 1e-8
        bits = obj.bits if (obj.bits is not None and np.isfinite(obj.bits) and obj.bits > 0) else 1e-8
        rows.append({'q': obj.q, 's': obj.s, 'bits': bits, 'L_qh': L})
    df = pd.DataFrame(rows)
    if df.empty:
        return df

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


def fit_log_model(sel_df: pd.DataFrame) -> Tuple[float, float]:
    if sel_df is None or len(sel_df) == 0:
        return 0.0, 0.0
    x = np.log10(sel_df['L_qh'])
    y = np.log10(sel_df['bits'])
    if len(x) < 2:
        return 0.0, 0.0
    a, b = np.polyfit(x, y, 1)
    return a, b


def blastp2RBNH(blast_file1: str, blast_file2: str, out_path: str, bin_size=1000, top_pct=0.05, min_Bprime=0.0) -> str:
    """将BLASTP结果转换为RBNH格式"""
    try:
        if (not config["overwrite"]["RBNH_compute"]) and os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            PipelineLogger.skip(f"RBNH file exists: {out_path}")
            return "skip"
        if not (os.path.exists(blast_file1) and os.path.exists(blast_file2)):
            PipelineLogger.warning(f"BLAST file missing: {blast_file1} or {blast_file2}")
            return "skip"
        
        PipelineLogger.info(f"Processing RBNH: {os.path.basename(out_path)}")
        
        blastp_dict1 = read_blast_file(blast_file1)
        blastp_dict2 = read_blast_file(blast_file2)
        blastp_merge_list = merge_max_bits(blastp_dict1, blastp_dict2)

        compute_L_qh(blastp_merge_list)

        df_binned = bin_top_hits(blastp_merge_list, bin_size, top_pct)
        a_coef, b_coef = fit_log_model(df_binned)

        for obj in blastp_merge_list:
            obj.compute_Bprime(a_coef, b_coef)

        seen = set()
        with open(out_path, 'w') as out_f:
            for obj in blastp_merge_list:
                if obj.q == obj.s:
                    continue
                key = frozenset([obj.q, obj.s])
                if key in seen:
                    continue
                bprime_val = obj.Bprime if (obj.Bprime is not None and np.isfinite(obj.Bprime)) else 0.0
                if bprime_val <= min_Bprime:
                    continue
                u, v = sorted([obj.q, obj.s])
                out_f.write(f"{u} {v} {bprime_val}\n")
                seen.add(key)

        PipelineLogger.success(f"RBNH completed: {os.path.basename(out_path)}")
        return "done"
    except Exception as e:
        PipelineLogger.error(f"RBNH failed ({os.path.basename(out_path)}): {e}")
        return "error"


def compute_RBNH(sample_list: List[str], blastp_dir: str, output_dir: str,
                 chromosome_list: List[str], thread: int) -> bool:
    """计算RBNH，返回是否成功"""
    PipelineLogger.start_stage("RBNH Computation")
    
    output_dir = FileProcessor.ensure_directory(output_dir)
    tasks = []
    RNBH_files = []

    def add_task(blast1, blast2, outpath):
        if not os.path.exists(blast1) or not os.path.exists(blast2):
            PipelineLogger.warning(f"Missing BLAST: {os.path.basename(blast1)} or {os.path.basename(blast2)}")
            return
        tasks.append((blast1, blast2, outpath))
        RNBH_files.append(outpath)

    # chr vs chr
    for s1, s2 in combinations(sample_list, 2):
        for chrom in chromosome_list:
            b1 = os.path.join(blastp_dir, f"{s1}_vs_{s2}_{chrom}{config['work_config']['blastp_suffix']}")
            b2 = os.path.join(blastp_dir, f"{s2}_vs_{s1}_{chrom}{config['work_config']['blastp_suffix']}")
            o = os.path.join(output_dir, f"{s1}_vs_{s2}_{chrom}{config['work_config']['RBNH_suffix']}")
            add_task(b1, b2, o)

    # scaffold vs chr
    for s1, s2 in combinations(sample_list, 2):
        for chrom in chromosome_list:
            b1 = os.path.join(blastp_dir, f"{s1}_vs_{s2}_scaffold_vs_{chrom}{config['work_config']['blastp_suffix']}")
            b2 = os.path.join(blastp_dir, f"{s2}_vs_{s1}_{chrom}_vs_scaffold{config['work_config']['blastp_suffix']}")
            o = os.path.join(output_dir, f"{s1}_vs_{s2}_scaffold_vs_{chrom}{config['work_config']['RBNH_suffix']}")
            add_task(b1, b2, o)

    # chr vs scaffold
    for s1, s2 in combinations(sample_list, 2):
        for chrom in chromosome_list:
            b1 = os.path.join(blastp_dir, f"{s1}_vs_{s2}_{chrom}_vs_scaffold{config['work_config']['blastp_suffix']}")
            b2 = os.path.join(blastp_dir, f"{s2}_vs_{s1}_scaffold_vs_{chrom}{config['work_config']['blastp_suffix']}")
            o = os.path.join(output_dir, f"{s1}_vs_{s2}_{chrom}_vs_scaffold{config['work_config']['RBNH_suffix']}")
            add_task(b1, b2, o)

    # scaffold vs scaffold
    for s1, s2 in combinations(sample_list, 2):
        b1 = os.path.join(blastp_dir, f"{s1}_vs_{s2}_scaffold{config['work_config']['blastp_suffix']}")
        b2 = os.path.join(blastp_dir, f"{s2}_vs_{s1}_scaffold{config['work_config']['blastp_suffix']}")
        o = os.path.join(output_dir, f"{s1}_vs_{s2}_scaffold{config['work_config']['RBNH_suffix']}")
        add_task(b1, b2, o)

    if not tasks:
        PipelineLogger.warning("No valid RBNH tasks to run")
        PipelineLogger.end_stage("RBNH Computation")
        return False

    PipelineLogger.info(f"Total RBNH pairs: {len(tasks)}")
    
    results = ProcessManager.run_parallel_tasks(
        tasks, 
        blastp2RBNH, 
        max_workers=thread,
        use_processes=False,
        task_description="RBNH Computation"
    )
    
    successful_tasks = ProcessManager.is_stage_successful(results)
    
    # 合并结果文件
    all_RNBH_files_path = os.path.join(output_dir, config["work_config"]["all_RBNH_file_name"])
    if not config["overwrite"]["RBNH_compute"] and os.path.exists(all_RNBH_files_path) and os.path.getsize(all_RNBH_files_path) > 0:
        PipelineLogger.skip(f"All RBNH file exists: {all_RNBH_files_path}")
    else:
        with open(all_RNBH_files_path, "w") as outfile:
            for f in RNBH_files:
                if not os.path.exists(f) or os.path.getsize(f) == 0:
                    continue
                with open(f, "r") as infile:
                    for line in infile:
                        outfile.write(line)
        PipelineLogger.success(f"All RBNH file created: {all_RNBH_files_path}")
    
    PipelineLogger.end_stage("RBNH Computation")
    return successful_tasks


def read_bed(filename: str) -> List[Gene]:
    """读取BED文件，支持更灵活的格式处理"""
    result = []
    try:
        with open(filename) as bed:
            for line_num, line in enumerate(bed, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split("\t")
                
                # 检查列数
                if len(parts) < 4:
                    PipelineLogger.warning(f"Skipping line {line_num} in {filename}: insufficient columns ({len(parts)})")
                    continue
                
                try:
                    chromosome = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3]
                    
                    # 处理可选的列
                    score = "."
                    strand = "+"
                    
                    if len(parts) >= 6:
                        score = parts[4] if parts[4] != "." else "."
                        strand = parts[5] if parts[5] in ["+", "-"] else "+"
                    elif len(parts) >= 5:
                        score = parts[4] if parts[4] != "." else "."
                    
                    result.append(Gene(chromosome, start, end, name, score, strand))
                    
                except (ValueError, IndexError) as e:
                    PipelineLogger.warning(f"Skipping line {line_num} in {filename}: {e}")
                    continue
                    
    except Exception as e:
        PipelineLogger.error(f"Error reading BED file {filename}: {e}")
        raise
    
    # 排序和索引
    result.sort(key=lambda x: (x.chromosome, x.start))

    thisChromosome = None
    thisIndex = 1
    for gene in result:
        if thisChromosome != gene.chromosome:
            thisChromosome = gene.chromosome
            thisIndex = 1
        gene.index = thisIndex
        thisIndex += 1
        
    PipelineLogger.info(f"Loaded {len(result)} genes from {filename}")
    return result


def genelist2info(genelist: List[Gene]) -> Dict[str, Gene]:
    result = {}
    for gene in genelist:
        result[gene.name] = gene
    return result


def split_pep_bed(bed_file: str, pep_file: str, pep_dir: str, bed_dir: str,
                  chromosome_list: List[str], sample_name: str) -> str:
    """分割PEP和BED文件"""
    chromosome_set = set(chromosome_list)

    if not config["overwrite"]["pep_bed_split"]:
        pep_exist = all(
            os.path.exists(os.path.join(pep_dir, f"{sample_name}_{chromosome}.pep")) and 
            os.path.getsize(os.path.join(pep_dir, f"{sample_name}_{chromosome}.pep")) > 0
            for chromosome in chromosome_set
        )
        bed_exist = all(
            os.path.exists(os.path.join(bed_dir, f"{sample_name}_{chromosome}.bed")) and 
            os.path.getsize(os.path.join(bed_dir, f"{sample_name}_{chromosome}.bed")) > 0
            for chromosome in chromosome_set
        )
        if pep_exist and bed_exist:
            PipelineLogger.skip(f"PEP/BED files already split for {sample_name}")
            return "skip"
    
    # 检查输入文件是否存在
    if not os.path.exists(pep_file) or not os.path.exists(bed_file):
        PipelineLogger.error(f"Input file missing: {pep_file} or {bed_file}")
        return "error"
    
    # 读取PEP数据
    pep_data = {}
    try:
        for record in SeqIO.parse(pep_file, "fasta"):
            pep_data[record.id] = str(record.seq)
    except Exception as e:
        PipelineLogger.error(f"Error reading PEP file {pep_file}: {e}")
        return "error"
    
    # 读取BED数据
    try:
        bed_data = read_bed(bed_file)
    except Exception as e:
        PipelineLogger.error(f"Error reading BED file {bed_file}: {e}")
        return "error"
    
    bed2info = {gene.name: gene for gene in bed_data}
    
    # 处理PEP文件
    pep_buffers = defaultdict(list)
    scaffold_pep = []

    for genename, pep in pep_data.items():
        if genename not in bed2info:
            PipelineLogger.warning(f"Gene {genename} not found in bed file {bed_file}")
            continue
        chromosome = bed2info[genename].chromosome
        new_genename = f"{sample_name}_{genename}"
        if chromosome not in chromosome_set:
            scaffold_pep.append(f">{new_genename}\n{pep}\n")
        else:
            pep_buffers[chromosome].append(f">{new_genename}\n{pep}\n")

    # 写PEP文件
    try:
        with open(os.path.join(pep_dir, f"{sample_name}_scaffold.pep"), "w") as out:
            out.writelines(scaffold_pep)
        for chrom, lines in pep_buffers.items():
            with open(os.path.join(pep_dir, f"{sample_name}_{chrom}.pep"), "w") as out:
                out.writelines(lines)
    except Exception as e:
        PipelineLogger.error(f"Error writing PEP files: {e}")
        return "error"
    
    # 处理BED文件
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
    elif config["bed_format"]["bed4"]:
        for gene in bed_data:
            new_name = f"{sample_name}_{gene.name}"
            line = f"{gene.chromosome}\t{gene.start}\t{gene.end}\t{new_name}\n"
            if gene.chromosome not in chromosome_set:
                scaffold_bed.append(line)
            else:
                bed_buffers[gene.chromosome].append(line)
    
    # 写BED文件
    try:
        with open(os.path.join(bed_dir, f"{sample_name}_scaffold.bed"), "w") as out:
            out.writelines(scaffold_bed)
        for chrom, lines in bed_buffers.items():
            with open(os.path.join(bed_dir, f"{sample_name}_{chrom}.bed"), "w") as out:
                out.writelines(lines)
    except Exception as e:
        PipelineLogger.error(f"Error writing BED files: {e}")
        return "error"
    
    return "done"


def blastp2dagchainer(blast_file: str, 
                      sample1: str, sample2: str,
                      bed_dict1: Dict[str, Gene], bed_dict2: Dict[str, Gene],
                      out_path: str):
    """Convert BLASTP results to DAGchainer input format"""
    blast_dict = read_blast_file(blast_file)

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
                         sample1: str, sample2: str,
                         bed_file1: str, bed_file2: str,
                         out_path1: str, out_path2: str) -> str:
    """生成DAGchainer输入文件"""
    if not config["overwrite"]["dagchainer_input"]:
        if os.path.exists(out_path1) and os.path.exists(out_path2) and os.path.getsize(out_path1) > 0 and os.path.getsize(out_path2) > 0:
            PipelineLogger.skip(f"DAGchainer input exists: {out_path1}")
            PipelineLogger.skip(f"DAGchainer input exists: {out_path2}")
            return "skip"

    PipelineLogger.info(f"Creating DAGchainer input: {os.path.basename(out_path1)}")
    
    # 检查输入文件是否存在
    if not os.path.exists(blast_file1) or not os.path.exists(blast_file2):
        PipelineLogger.warning(f"BLAST file missing: {blast_file1} or {blast_file2}")
        return "skip"
    if not os.path.exists(bed_file1) or not os.path.exists(bed_file2):
        PipelineLogger.warning(f"BED file missing: {bed_file1} or {bed_file2}")
        return "skip"
    
    try:
        bed_dict1 = {g.name: g for g in read_bed(bed_file1)}
        bed_dict2 = {g.name: g for g in read_bed(bed_file2)}
    except Exception as e:
        PipelineLogger.error(f"Error reading BED files: {e}")
        return "error"

    try:
        blastp2dagchainer(blast_file1, sample1, sample2, bed_dict1, bed_dict2, out_path1)
        blastp2dagchainer(blast_file2, sample2, sample1, bed_dict2, bed_dict1, out_path2)
    except Exception as e:
        PipelineLogger.error(f"Error creating DAGchainer input: {e}")
        return "error"
    
    PipelineLogger.success(f"DAGchainer input created: {os.path.basename(out_path1)}")
    return "done"


def read_synteny_dagchainer(filename: str, bed_dict_left: Dict[str, Gene] = None, bed_dict_right: Dict[str, Gene] = None) -> List[Block]:
    """读取DAGchainer输出的同源区块"""
    result = []
    try:
        for i in open(filename):
            i = i.strip()
            if not i:
                continue
            if i.startswith("#"):
                result.append([])
            else:
                line = i.split("\t")
                chrom_left, gene_left = line[0], line[1]
                chrom_right, gene_right = line[4], line[5]
                if not result:
                    result.append([])
                result[-1].append([chrom_left, gene_left, chrom_right, gene_right])
    except Exception as e:
        PipelineLogger.error(f"Error reading DAGchainer file {filename}: {e}")
        return []
    
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
            PipelineLogger.warning(f"Gene {gene_left} or {gene_right} not found in bed file")
            continue
        
        if swap_flag:
            for this_line in this_block:
                this_line[0], this_line[2] = this_line[2], this_line[0]
                this_line[1], this_line[3] = this_line[3], this_line[1]
    
    final_result = []
    for this_block in result:
        if not this_block:
            continue
        chromosome1, chromosome2 = this_block[0][0], this_block[0][2]
        genelist1 = []
        positions1 = []
        genelist2 = []
        positions2 = []
        for line in this_block:
            if line[1] not in bed_dict_left or line[3] not in bed_dict_right:
                continue
            genelist1.append(bed_dict_left[line[1]].name)
            positions1.append(bed_dict_left[line[1]].index)
            genelist2.append(bed_dict_right[line[3]].name)
            positions2.append(bed_dict_right[line[3]].index)
        if genelist1 and genelist2:
            final_result.append(Block(chromosome1, genelist1, positions1, chromosome2, genelist2, positions2))
    return final_result


def filter_synteny_blocks(synteny_blocks: List[Block], chromosome_list: List[str]) -> List[Tuple[str, str]]:
    """过滤同源区块"""
    chromosome_dict = {chromosome: index for index, chromosome in enumerate(chromosome_list)}
    synteny_blocks.sort(key=lambda x: (
        not (x.chromosome1 in chromosome_dict and x.chromosome2 in chromosome_dict),
        chromosome_dict.get(x.chromosome1, float("inf")),
        chromosome_dict.get(x.chromosome2, float("inf")),
        -x.hits
    ))
    
    result_pair = []
    block_tree_left = {}
    block_tree_right = {}

    for block_index, this_block in enumerate(synteny_blocks):
        if this_block.chromosome1 not in block_tree_left:
            block_tree_left[this_block.chromosome1] = Tree()
        if this_block.chromosome2 not in block_tree_right:
            block_tree_right[this_block.chromosome2] = Tree()

        this_block_left_start = min(this_block.positions1)
        this_block_left_end = max(this_block.positions1)
        this_block_right_start = min(this_block.positions2)
        this_block_right_end = max(this_block.positions2)
        
        if block_index == 0:
            block_tree_left[this_block.chromosome1].add_closed(this_block_left_start, this_block_left_end)
            block_tree_right[this_block.chromosome2].add_closed(this_block_right_start, this_block_right_end)
            for gene1, gene2 in zip(this_block.genelist1, this_block.genelist2):
                result_pair.append((gene1, gene2))
            continue
        
        left_segments_uncovered = block_tree_left[this_block.chromosome1].get_uncovered(this_block_left_start, this_block_left_end)
        right_segments_uncovered = block_tree_right[this_block.chromosome2].get_uncovered(this_block_right_start, this_block_right_end)

        new_left_segments_uncovered, new_right_segments_uncovered = this_block.adjust_segments(left_segments_uncovered, right_segments_uncovered)

        new_block = this_block.intersection(new_left_segments_uncovered, new_right_segments_uncovered)
        block_tree_left[this_block.chromosome1].add_closed(this_block_left_start, this_block_left_end)
        block_tree_right[this_block.chromosome2].add_closed(this_block_right_start, this_block_right_end)
        
        if new_block.length() == 0:
            continue
        for gene1, gene2 in zip(new_block.genelist1, new_block.genelist2):
            result_pair.append((gene1, gene2))
    
    return result_pair

def geneindex_removeSampleName(geneindex:list[list[str]],sample_list: list[str])-> list[list[str]]:
    # 检查数据
    if len(geneindex[0])!=len(sample_list):
        raise Exception(f"geneindex columns {len(geneindex[0])} not equal to sample num {len(sample_list)}")
    sample_name_len=[len(sample_name) for sample_name in sample_list]
    result=[]
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
    return result


def filter_synteny_blocks_normal(synteny_blocks: List[Block], chromosome_list: List[str]) -> List[Tuple[str, str]]:
    """正常过滤同源区块"""
    chromosome_dict = {chromosome: index for index, chromosome in enumerate(chromosome_list)}
    synteny_blocks.sort(key=lambda x: (
        not (x.chromosome1 in chromosome_dict and x.chromosome2 in chromosome_dict),
        chromosome_dict.get(x.chromosome1, float("inf")),
        chromosome_dict.get(x.chromosome2, float("inf")),
        -x.hits
    ))
    
    left_gene_set = set()
    right_gene_set = set()
    result_pair = []
    
    for this_block in synteny_blocks:
        if not this_block.is_normal():
            continue
        for gene1, gene2 in zip(this_block.genelist1, this_block.genelist2):
            if gene1 in left_gene_set or gene2 in right_gene_set:
                continue
            left_gene_set.add(gene1)
            right_gene_set.add(gene2)
            result_pair.append((gene1, gene2))
    
    return result_pair


def sample1_vs_sample2_pair(sample1: str, sample2: str,
                           bed_dir: str,
                           chromosome_list: List[str],
                           dagchainer_output_dir: str,
                           filter_method: str = config["gene_filter_method"]["default"]) -> Tuple[str, str, List[Tuple[str, str]]]:
    """样本对之间的同源基因对分析，增加错误处理"""
    try:
        filtered_func = globals().get(filter_method)
        if filtered_func is None:
            raise ValueError(f"Unknown filter method: {filter_method}")
        
        ref_bed_list = []
        qry_bed_list = [] 
        ref_blocks_list = []
        qry_blocks_list = []
        
        PipelineLogger.info(f"Processing pair: {sample1} vs {sample2}")
        
        # 读取BED文件，增加文件存在性检查
        for chro in chromosome_list:
            bed_file1 = os.path.join(bed_dir, f"{sample1}_{chro}.bed")
            bed_file2 = os.path.join(bed_dir, f"{sample2}_{chro}.bed")
            
            if os.path.exists(bed_file1):
                try:
                    ref_bed_list.extend(read_bed(bed_file1))
                except Exception as e:
                    PipelineLogger.warning(f"Error reading {bed_file1}: {e}")
            else:
                PipelineLogger.warning(f"BED file not found: {bed_file1}")
                
            if os.path.exists(bed_file2):
                try:
                    qry_bed_list.extend(read_bed(bed_file2))
                except Exception as e:
                    PipelineLogger.warning(f"Error reading {bed_file2}: {e}")
            else:
                PipelineLogger.warning(f"BED file not found: {bed_file2}")
        
        # 读取scaffold文件
        scaffold_bed1 = os.path.join(bed_dir, f"{sample1}_scaffold.bed")
        scaffold_bed2 = os.path.join(bed_dir, f"{sample2}_scaffold.bed")
        
        if os.path.exists(scaffold_bed1):
            try:
                ref_bed_list.extend(read_bed(scaffold_bed1))
            except Exception as e:
                PipelineLogger.warning(f"Error reading {scaffold_bed1}: {e}")
        if os.path.exists(scaffold_bed2):
            try:
                qry_bed_list.extend(read_bed(scaffold_bed2))
            except Exception as e:
                PipelineLogger.warning(f"Error reading {scaffold_bed2}: {e}")

        if not ref_bed_list or not qry_bed_list:
            PipelineLogger.warning(f"No BED data found for {sample1} or {sample2}")
            return (sample1, sample2, [])

        ref_bed_dict = genelist2info(ref_bed_list)
        qry_bed_dict = genelist2info(qry_bed_list)
        
        PipelineLogger.info(f"Loaded {len(ref_bed_dict)} genes from {sample1}, {len(qry_bed_dict)} genes from {sample2}")
        
        # 读取DAGchainer结果，增加文件存在性检查
        dagchainer_files_processed = 0
        
        for chro in chromosome_list:
            # chro vs chro
            dagchainer_file1 = os.path.join(dagchainer_output_dir, f"{sample1}_vs_{sample2}_{chro}{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
            dagchainer_file2 = os.path.join(dagchainer_output_dir, f"{sample2}_vs_{sample1}_{chro}{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
            
            if os.path.exists(dagchainer_file1):
                try:
                    ref_blocks_list.extend(read_synteny_dagchainer(dagchainer_file1, bed_dict_left=ref_bed_dict, bed_dict_right=qry_bed_dict))
                    dagchainer_files_processed += 1
                except Exception as e:
                    PipelineLogger.warning(f"Error reading {dagchainer_file1}: {e}")
                    
            if os.path.exists(dagchainer_file2):
                try:
                    qry_blocks_list.extend(read_synteny_dagchainer(dagchainer_file2, bed_dict_left=ref_bed_dict, bed_dict_right=qry_bed_dict))
                    dagchainer_files_processed += 1
                except Exception as e:
                    PipelineLogger.warning(f"Error reading {dagchainer_file2}: {e}")
            
            # chro vs scaffold
            dagchainer_file1 = os.path.join(dagchainer_output_dir, f"{sample1}_vs_{sample2}_{chro}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
            dagchainer_file2 = os.path.join(dagchainer_output_dir, f"{sample2}_vs_{sample1}_{chro}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
            
            if os.path.exists(dagchainer_file1):
                try:
                    ref_blocks_list.extend(read_synteny_dagchainer(dagchainer_file1, bed_dict_left=ref_bed_dict, bed_dict_right=qry_bed_dict))
                    dagchainer_files_processed += 1
                except Exception as e:
                    PipelineLogger.warning(f"Error reading {dagchainer_file1}: {e}")
                    
            if os.path.exists(dagchainer_file2):
                try:
                    qry_blocks_list.extend(read_synteny_dagchainer(dagchainer_file2, bed_dict_left=ref_bed_dict, bed_dict_right=qry_bed_dict))
                    dagchainer_files_processed += 1
                except Exception as e:
                    PipelineLogger.warning(f"Error reading {dagchainer_file2}: {e}")

            # scaffold vs chro
            dagchainer_file1 = os.path.join(dagchainer_output_dir, f"{sample1}_vs_{sample2}_scaffold_vs_{chro}{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
            dagchainer_file2 = os.path.join(dagchainer_output_dir, f"{sample2}_vs_{sample1}_scaffold_vs_{chro}{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
            
            if os.path.exists(dagchainer_file1):
                try:
                    ref_blocks_list.extend(read_synteny_dagchainer(dagchainer_file1, bed_dict_left=ref_bed_dict, bed_dict_right=qry_bed_dict))
                    dagchainer_files_processed += 1
                except Exception as e:
                    PipelineLogger.warning(f"Error reading {dagchainer_file1}: {e}")
                    
            if os.path.exists(dagchainer_file2):
                try:
                    qry_blocks_list.extend(read_synteny_dagchainer(dagchainer_file2, bed_dict_left=ref_bed_dict, bed_dict_right=qry_bed_dict))
                    dagchainer_files_processed += 1
                except Exception as e:
                    PipelineLogger.warning(f"Error reading {dagchainer_file2}: {e}")

        # scaffold vs scaffold
        dagchainer_file1 = os.path.join(dagchainer_output_dir, f"{sample1}_vs_{sample2}_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
        dagchainer_file2 = os.path.join(dagchainer_output_dir, f"{sample2}_vs_{sample1}_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
        
        if os.path.exists(dagchainer_file1):
            try:
                ref_blocks_list.extend(read_synteny_dagchainer(dagchainer_file1, bed_dict_left=ref_bed_dict, bed_dict_right=qry_bed_dict))
                dagchainer_files_processed += 1
            except Exception as e:
                PipelineLogger.warning(f"Error reading {dagchainer_file1}: {e}")
                
        if os.path.exists(dagchainer_file2):
            try:
                qry_blocks_list.extend(read_synteny_dagchainer(dagchainer_file2, bed_dict_left=ref_bed_dict, bed_dict_right=qry_bed_dict))
                dagchainer_files_processed += 1
            except Exception as e:
                PipelineLogger.warning(f"Error reading {dagchainer_file2}: {e}")

        PipelineLogger.info(f"Processed {dagchainer_files_processed} DAGchainer files for {sample1} vs {sample2}")

        filtered_blocks = filtered_func([*ref_blocks_list, *qry_blocks_list], chromosome_list)
        PipelineLogger.info(f"Found {len(filtered_blocks)} synteny pairs for {sample1} vs {sample2}")
        
        return (sample1, sample2, filtered_blocks)
        
    except Exception as e:
        PipelineLogger.error(f"Error in sample1_vs_sample2_pair for {sample1} vs {sample2}: {e}")
        # 返回空结果而不是抛出异常，让流程继续
        return (sample1, sample2, [])


def combine_SG(data: List[Tuple[str, str, List[Tuple[str, str]]]], 
               chromosome_list: List[str], sample_list: List[str], 
               bed_dir: str, output_dir: str) -> List[Tuple[str, str]]:
    """合并同源基因对生成geneindex，增加空数据检查"""
    PipelineLogger.start_stage("GeneIndex Generation")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    outfile_path = os.path.join(output_dir, config['work_config']['geneindex_output_file_name'])
    outfile_without_sample_path = os.path.join(output_dir, config['work_config']['geneindex_without_sample_name'])
    if not config["overwrite"]["geneindex_output"] and os.path.exists(outfile_path) and os.path.getsize(outfile_path) > 0 and \
       os.path.exists(outfile_without_sample_path) and os.path.getsize(outfile_without_sample_path) > 0:
        PipelineLogger.skip(f"GeneIndex file exists: {outfile_path} and {outfile_without_sample_path}")
        return []

    # 检查是否有有效数据
    if not data:
        PipelineLogger.warning("No synteny data available for GeneIndex generation")
        # 创建空的geneindex文件
        with open(outfile_path, "w") as f:
            f.write("# Empty GeneIndex - no synteny pairs found\n")
        PipelineLogger.end_stage("GeneIndex Generation")
        return []

    sample2index = {sample: i for i, sample in enumerate(sample_list)}
    sample_size = len(sample_list)

    # 构建图
    G = nx.Graph()
    total_pairs = 0
    for sample1, sample2, pairs in data:
        total_pairs += len(pairs)
        for gene1, gene2 in pairs:
            G.add_node(gene1, sample=sample1)
            G.add_node(gene2, sample=sample2)
            G.add_edge(gene1, gene2)

    PipelineLogger.info(f"Built graph with {len(G.nodes())} nodes and {len(G.edges())} edges from {total_pairs} pairs")

    visited = set()
    SG_list = []
    result_edges = []

    # 寻找连通分量
    for node in G.nodes():
        if node in visited:
            continue
        try:
            tree = nx.dfs_tree(G, node)
            path_nodes = list(tree.nodes)
            path_edges = list(tree.edges)
            
            visited.update(path_nodes)

            SG_line = [[] for _ in range(sample_size)]
            for gene in path_nodes:
                sample = G.nodes[gene]["sample"]
                SG_line[sample2index[sample]].append(gene)
            
            flag_duplicate = False
            a = []
            for i in SG_line:
                if not i:
                    a.append("-")
                else:
                    a.append(",".join(i))
                    if len(i) > 1:
                        flag_duplicate = True
            SG_list.append(a)
            if flag_duplicate:
                result_edges.append(path_edges)
                
        except Exception as e:
            PipelineLogger.warning(f"Error processing node {node}: {e}")
            continue

    # 添加未匹配的基因
    PipelineLogger.info("Adding unmatched genes...")
    for sample_index, sample in enumerate(sample_list):
        this_sample_gene_set = set()
        for genes in [line[sample_index] for line in SG_list]:
            if genes == "-":
                continue
            elif "," in genes:
                for gene in genes.split(","):
                    this_sample_gene_set.add(gene)
            else:
                this_sample_gene_set.add(genes)
        
        bed_data = []
        for chromosome in chromosome_list:
            bed_path = os.path.join(bed_dir, f"{sample}_{chromosome}.bed")
            if os.path.exists(bed_path):
                try:
                    bed_data.extend(read_bed(bed_path))
                except Exception as e:
                    PipelineLogger.warning(f"Error reading {bed_path}: {e}")
                    
        scaffold_bed_path = os.path.join(bed_dir, f"{sample}_scaffold.bed")
        if os.path.exists(scaffold_bed_path):
            try:
                bed_data.extend(read_bed(scaffold_bed_path))
            except Exception as e:
                PipelineLogger.warning(f"Error reading {scaffold_bed_path}: {e}")
                
        unmatched_count = 0
        for gene in bed_data:
            if gene.name not in this_sample_gene_set:
                this_line = ["-" for _ in range(sample_size)]
                this_line[sample_index] = gene.name
                SG_list.append(this_line)
                unmatched_count += 1
        
        PipelineLogger.info(f"Added {unmatched_count} unmatched genes for {sample}")

    # 写文件
    with open(outfile_path, "w") as out:
        for line in SG_list:
            out.write("\t".join(line) + "\n")
    with open(outfile_without_sample_path, "w") as out:
        SG_list_no_sample = geneindex_removeSampleName(SG_list, sample_list)
        for line in SG_list_no_sample:
            out.write("\t".join(line) + "\n")

    PipelineLogger.success(f"GeneIndex created: {outfile_path} with {len(SG_list)} rows")
    PipelineLogger.end_stage("GeneIndex Generation")
    

    return result_edges

def read_mcl_cluster(mcl_output_path: str) -> Dict[str, int]:
    """读取MCL聚类结果"""
    gene2cluster_dict = {}
    with open(mcl_output_path) as f:
        for cluster_index, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            genes = line.split("\t")
            for gene in genes:
                gene2cluster_dict[gene] = cluster_index + 1
    return gene2cluster_dict


def _process_cluster_chunk(args):
    """处理聚类分块的辅助函数"""
    all_RBNH_file_path, clusters_chunk, out_file_dir, chunk_index = args
    
    # 读取RBNH图
    G = nx.Graph()
    try:
        with open(all_RBNH_file_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                node1, node2, weight = line.split(" ")
                G.add_edge(node1, node2, weight=float(weight))
    except Exception as e:
        PipelineLogger.error(f"Error reading RBNH file: {e}")
        return []

    results = []
    base_idx = chunk_index * len(clusters_chunk)

    for i, cluster in enumerate(clusters_chunk):
        cluster_id = f"cluster-{base_idx + i + 1}"
        out_file = os.path.join(out_file_dir, f"{cluster_id}{config['work_config']['mcl_input_suffix']}")

        try:
            with open(out_file, "w") as f:
                for node1, node2 in cluster:
                    if G.has_edge(node1, node2):
                        w = G[node1][node2]["weight"]
                        f.write(f"{node1} {node2} {w}\n")
            results.append(out_file)
        except Exception as e:
            PipelineLogger.error(f"Error writing MCL input file {out_file}: {e}")

    return results


def get_duplicate_RBNH_parallel(all_RBNH_file_path: str, clusters_list: List[List[Tuple[str, str]]], 
                               out_file_dir: str, n_processes: int) -> List[str]:
    """并行处理重复的RBNH"""
    PipelineLogger.start_stage("Duplicate RBNH Processing")
    
    if not os.path.exists(out_file_dir):
        os.makedirs(out_file_dir)

    if not clusters_list:
        PipelineLogger.warning("No clusters to process")
        PipelineLogger.end_stage("Duplicate RBNH Processing")
        return []

    chunk_size = (len(clusters_list) + n_processes - 1) // n_processes
    chunks = [clusters_list[i:i+chunk_size] for i in range(0, len(clusters_list), chunk_size)]

    tasks = [(all_RBNH_file_path, chunk, out_file_dir, idx) for idx, chunk in enumerate(chunks)]
    
    PipelineLogger.info(f"Processing {len(clusters_list)} clusters with {n_processes} processes")
    
    results = []
    try:
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            future_to_chunk = {executor.submit(_process_cluster_chunk, task): i 
                            for i, task in enumerate(tasks)}
            
            for future in as_completed(future_to_chunk):
                chunk_results = future.result()
                results.extend(chunk_results)
                PipelineLogger.info(f"Processed chunk {future_to_chunk[future] + 1}/{len(tasks)}")
    except Exception as e:
        PipelineLogger.error(f"Error in parallel processing: {e}")
        return []

    PipelineLogger.success(f"Generated {len(results)} MCL input files")
    PipelineLogger.end_stage("Duplicate RBNH Processing")
    
    return results


def run_mcl_split(gene2cluster_dict: Dict[str, int], geneindex_in_path: str, geneindex_out_path: str,geneindex_out_path2: str,sample_list:list[str]) -> None:
    """使用MCL聚类结果分割 geneindex"""
    PipelineLogger.info(f"Splitting GeneIndex using MCL clusters")

    result_geneindex_split=[]
    try:
        with open(geneindex_in_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                columns = line.split("\t")
                result_geneindex_split.append(columns)
                if not any("," in genes for genes in columns):
                    continue

                cluster2gene = defaultdict(lambda: defaultdict(list))
                for idx, genes in enumerate(columns):
                    if genes == "-":
                        continue
                    elif "," in genes:
                        for gene in genes.split(","):
                            cluster_index = gene2cluster_dict.get(gene, gene)
                            cluster2gene[cluster_index][idx].append(gene)
                    else:
                        cluster_index = gene2cluster_dict.get(genes, genes)
                        cluster2gene[cluster_index][idx].append(genes)
                
                all_line_write = []
                for cluster_index, sg_dict in cluster2gene.items():
                    all_line = ["-" for _ in range(len(columns))]
                    for idx, gene_list in sg_dict.items():
                        all_line[idx] = ",".join(gene_list)
                    all_line_write.append(all_line)                
                result_geneindex_split.extend(all_line_write)
        # 写入分割后的geneindex文件
        with open(geneindex_out_path, "w") as out:
            for line in result_geneindex_split:
                out.write("\t".join(line) + "\n")

    except Exception as e:
        PipelineLogger.error(f"Error splitting GeneIndex: {e}")
        return
    result_geneindex_split_without_sample=geneindex_removeSampleName(result_geneindex_split, sample_list)
    with open(geneindex_out_path2,"w")as out2:
        for line in result_geneindex_split_without_sample:
            out2.write("\t".join(line)+"\n")
    

    PipelineLogger.success(f"Split GeneIndex created: {geneindex_out_path} and {geneindex_out_path2}")


def run_all_mcl_split(mcl_output_path: str, geneindex_output_dir: str,sample_list:list[str]) -> bool:
    """运行完整的MCL分割流程，返回是否成功"""
    PipelineLogger.start_stage("MCL-based GeneIndex Splitting")
    
    if not os.path.exists(mcl_output_path) or os.path.getsize(mcl_output_path) == 0:
        PipelineLogger.warning(f"MCL output file not found or empty: {mcl_output_path}")
        PipelineLogger.end_stage("MCL-based GeneIndex Splitting")
        return False
    
    try:
        gene2cluster_dict = read_mcl_cluster(mcl_output_path)
    except Exception as e:
        PipelineLogger.error(f"Error reading MCL clusters: {e}")
        PipelineLogger.end_stage("MCL-based GeneIndex Splitting")
        return False
    
    geneindex_in_path = os.path.join(geneindex_output_dir, config['work_config']['geneindex_output_file_name'])
    geneindex_out_path = os.path.join(geneindex_output_dir, config['work_config']['geneindex_output_file_split_name'])
    geneindex_out_path2= os.path.join(geneindex_output_dir, config['work_config']['geneindex_output_file_split_name_without_sample'])

    if not os.path.exists(geneindex_in_path):
        PipelineLogger.error(f"GeneIndex input file not found: {geneindex_in_path}")
        PipelineLogger.end_stage("MCL-based GeneIndex Splitting")
        return False
    
    run_mcl_split(gene2cluster_dict, geneindex_in_path, geneindex_out_path,geneindex_out_path2, sample_list)
    
    PipelineLogger.end_stage("MCL-based GeneIndex Splitting")
    return True


# 主要流程函数
def run_pep_bed_split(input_dir: str, working_pep_dir: str, working_bed_dir: str,
                     chromosome_list: List[str], sample_list: List[str], thread: int) -> bool:
    """运行PEP和BED文件分割，返回是否成功"""
    PipelineLogger.start_stage("PEP/BED File Splitting")
    
    tasks = []
    for sample in sample_list:
        bed_file = os.path.join(input_dir, sample + config["input_file_type"]["bed_suffix"])
        pep_file = os.path.join(input_dir, sample + config["input_file_type"]["pep_suffix"])
        tasks.append((bed_file, pep_file, working_pep_dir, working_bed_dir, chromosome_list, sample))
    
    results = ProcessManager.run_parallel_tasks(
        tasks, split_pep_bed, max_workers=thread, 
        task_description="Splitting PEP/BED files"
    )
    
    successful_tasks = ProcessManager.is_stage_successful(results)
    PipelineLogger.info(f"PEP/BED splitting: {results['success']} success, {results['skipped']} skipped, {results['errors']} errors")
    
    PipelineLogger.end_stage("PEP/BED File Splitting")
    return successful_tasks


def run_diamond_analysis(working_pep_dir: str, working_diamond_dir: str,
                        chromosome_list: List[str], sample_list: List[str], 
                        diamond_path: str, thread: int) -> bool:
    """运行Diamond分析，返回是否成功（跳过不算失败）"""
    PipelineLogger.start_stage("Diamond Analysis")
    
    # 检查diamond是否可用
    if not diamond_path:
        PipelineLogger.warning("Diamond not available, skipping Diamond analysis")
        PipelineLogger.end_stage("Diamond Analysis")
        return True  # 工具不可用视为跳过，不是错误
    
    # 创建数据库任务
    PipelineLogger.info("Creating Diamond databases...")
    tasks_makedb = []
    args_makedb = [str(item) for kv in config["diamond_makedb"].items() for item in kv]
    
    for chromosome in [*chromosome_list, "scaffold"]:
        for sample in sample_list:
            pep_path = os.path.join(working_pep_dir, f"{sample}_{chromosome}.pep")
            # 检查PEP文件是否存在
            if not os.path.exists(pep_path) or os.path.getsize(pep_path) == 0:
                PipelineLogger.warning(f"PEP file not found or empty: {pep_path}, skipping")
                continue
            diamond_db_path = os.path.join(working_diamond_dir, f"{sample}_{chromosome}{config['work_config']['blastp_db_suffix']}")
            tasks_makedb.append((diamond_path, pep_path, diamond_db_path, args_makedb))
    
    if not tasks_makedb:
        PipelineLogger.warning("No valid PEP files found for Diamond database creation")
        PipelineLogger.end_stage("Diamond Analysis")
        return True  # 没有任务视为跳过，不是错误
    
    results_makedb = ProcessManager.run_parallel_tasks(
        tasks_makedb, DiamondManager.run_diamond_makedb, max_workers=thread,
        task_description="Creating Diamond databases"
    )
    
    # 检查数据库创建结果 - 只要有成功或跳过，没有错误，就认为是成功的
    diamond_db_success = ProcessManager.is_stage_successful(results_makedb)
    
    if not diamond_db_success:
        PipelineLogger.error("Diamond database creation failed with errors")
        PipelineLogger.end_stage("Diamond Analysis")
        return False
    
    # BLASTP任务
    PipelineLogger.info("Running Diamond BLASTP...")
    tasks_blastp = []
    args_blastp = []
    for k, v in config["diamond_blastp"].items():
        args_blastp.append(k)
        if isinstance(v, list):
            args_blastp.extend(v)
        else:
            args_blastp.append(str(v))
    
    # 准备BLASTP任务
    for sample1 in sample_list: 
        for sample2 in sample_list:
            if sample1 == sample2:
                continue
            
            # chro vs chro
            for chromosome in chromosome_list:
                qry = os.path.join(working_pep_dir, f"{sample1}_{chromosome}.pep")
                db = os.path.join(working_diamond_dir, f"{sample2}_{chromosome}{config['work_config']['blastp_db_suffix']}")
                if os.path.exists(qry) and os.path.exists(db):
                    out = os.path.join(working_diamond_dir, f"{sample1}_vs_{sample2}_{chromosome}{config['work_config']['blastp_suffix']}")
                    tasks_blastp.append((diamond_path, qry, db, out, args_blastp))
            
            # scaffold vs chro
            for chromosome in chromosome_list:
                qry = os.path.join(working_pep_dir, f"{sample1}_{chromosome}.pep")
                db = os.path.join(working_diamond_dir, f"{sample2}_scaffold{config['work_config']['blastp_db_suffix']}")
                if os.path.exists(qry) and os.path.exists(db):
                    out = os.path.join(working_diamond_dir, f"{sample1}_vs_{sample2}_{chromosome}_vs_scaffold{config['work_config']['blastp_suffix']}")
                    tasks_blastp.append((diamond_path, qry, db, out, args_blastp))
                    
                qry = os.path.join(working_pep_dir, f"{sample1}_scaffold.pep")
                db = os.path.join(working_diamond_dir, f"{sample2}_{chromosome}{config['work_config']['blastp_db_suffix']}")
                if os.path.exists(qry) and os.path.exists(db):
                    out = os.path.join(working_diamond_dir, f"{sample1}_vs_{sample2}_scaffold_vs_{chromosome}{config['work_config']['blastp_suffix']}")
                    tasks_blastp.append((diamond_path, qry, db, out, args_blastp))

            # scaffold vs scaffold
            qry = os.path.join(working_pep_dir, f"{sample1}_scaffold.pep")
            db = os.path.join(working_diamond_dir, f"{sample2}_scaffold{config['work_config']['blastp_db_suffix']}")
            if os.path.exists(qry) and os.path.exists(db):
                out = os.path.join(working_diamond_dir, f"{sample1}_vs_{sample2}_scaffold{config['work_config']['blastp_suffix']}")
                tasks_blastp.append((diamond_path, qry, db, out, args_blastp))
    
    if not tasks_blastp:
        PipelineLogger.warning("No valid BLASTP tasks to run")
        PipelineLogger.end_stage("Diamond Analysis")
        return True  # 没有BLASTP任务视为跳过，不是错误
    
    results_blastp = ProcessManager.run_parallel_tasks(
        tasks_blastp, 
        lambda diamond, qry, db, out, args: DiamondManager.run_diamond_blastp(diamond, qry, db, out, args),
        max_workers=thread,
        task_description="Running Diamond BLASTP"
    )
    
    # 检查BLASTP结果 - 只要有成功或跳过，没有错误，就认为是成功的
    blastp_success = ProcessManager.is_stage_successful(results_blastp)
    
    PipelineLogger.info(f"Diamond analysis completed: {results_blastp['success']} success, {results_blastp['skipped']} skipped, {results_blastp['errors']} errors")
    PipelineLogger.end_stage("Diamond Analysis")
    
    return blastp_success


def run_dagchainer_analysis(working_diamond_dir: str, working_bed_dir: str, working_dagchainer_dir: str,
                           chromosome_list: List[str], sample_list: List[str], thread: int) -> str:
    """运行DAGchainer分析，返回输出目录路径"""
    PipelineLogger.start_stage("DAGchainer Analysis")
    
    perl_path = DAGchainerManager.get_perl_path()
    dagchainer_path = DAGchainerManager.get_dagchainer_path()
    
    # 检查必要工具是否可用
    if not perl_path or not dagchainer_path:
        PipelineLogger.error("Perl or DAGchainer not available, skipping DAGchainer analysis")
        PipelineLogger.end_stage("DAGchainer Analysis")
        return FileProcessor.ensure_directory(os.path.join(working_dagchainer_dir, "output"))
    
    # 生成DAGchainer输入文件
    PipelineLogger.info("Generating DAGchainer input files...")
    tasks_dagchainer_input = []
    
    for sample1 in sample_list: 
        for sample2 in sample_list:
            if sample1 == sample2:
                continue
            # chro vs chro
            for chromosome in chromosome_list:
                blastp_out1 = os.path.join(working_diamond_dir, f"{sample1}_vs_{sample2}_{chromosome}{config['work_config']['blastp_suffix']}")
                blastp_out2 = os.path.join(working_diamond_dir, f"{sample2}_vs_{sample1}_{chromosome}{config['work_config']['blastp_suffix']}")
                bed_file1 = os.path.join(working_bed_dir, f"{sample1}_{chromosome}.bed")
                bed_file2 = os.path.join(working_bed_dir, f"{sample2}_{chromosome}.bed")
                dagchainer_input1 = os.path.join(working_dagchainer_dir, f"{sample1}_vs_{sample2}_{chromosome}{config['work_config']['dagchainer_input_suffix']}")
                dagchainer_input2 = os.path.join(working_dagchainer_dir, f"{sample2}_vs_{sample1}_{chromosome}{config['work_config']['dagchainer_input_suffix']}")
                tasks_dagchainer_input.append((blastp_out1, blastp_out2, sample1, sample2, bed_file1, bed_file2, dagchainer_input1, dagchainer_input2))
            
            # scaffold vs chro
            for chromosome in chromosome_list:
                blastp_out1 = os.path.join(working_diamond_dir, f"{sample1}_vs_{sample2}_scaffold_vs_{chromosome}{config['work_config']['blastp_suffix']}")
                blastp_out2 = os.path.join(working_diamond_dir, f"{sample2}_vs_{sample1}_{chromosome}_vs_scaffold{config['work_config']['blastp_suffix']}")
                bed_file1 = os.path.join(working_bed_dir, f"{sample1}_scaffold.bed")
                bed_file2 = os.path.join(working_bed_dir, f"{sample2}_{chromosome}.bed")
                dagchainer_input1 = os.path.join(working_dagchainer_dir, f"{sample1}_vs_{sample2}_scaffold_vs_{chromosome}{config['work_config']['dagchainer_input_suffix']}")
                dagchainer_input2 = os.path.join(working_dagchainer_dir, f"{sample2}_vs_{sample1}_{chromosome}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}")
                tasks_dagchainer_input.append((blastp_out1, blastp_out2, sample1, sample2, bed_file1, bed_file2, dagchainer_input1, dagchainer_input2))
            
            # chro vs scaffold
            for chromosome in chromosome_list:
                blastp_out1 = os.path.join(working_diamond_dir, f"{sample1}_vs_{sample2}_{chromosome}_vs_scaffold{config['work_config']['blastp_suffix']}")
                blastp_out2 = os.path.join(working_diamond_dir, f"{sample2}_vs_{sample1}_scaffold_vs_{chromosome}{config['work_config']['blastp_suffix']}")
                bed_file1 = os.path.join(working_bed_dir, f"{sample1}_{chromosome}.bed")
                bed_file2 = os.path.join(working_bed_dir, f"{sample2}_scaffold.bed")
                dagchainer_input1 = os.path.join(working_dagchainer_dir, f"{sample1}_vs_{sample2}_{chromosome}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}")
                dagchainer_input2 = os.path.join(working_dagchainer_dir, f"{sample2}_vs_{sample1}_scaffold_vs_{chromosome}{config['work_config']['dagchainer_input_suffix']}")
                tasks_dagchainer_input.append((blastp_out1, blastp_out2, sample1, sample2, bed_file1, bed_file2, dagchainer_input1, dagchainer_input2))

            # scaffold vs scaffold
            blastp_out1 = os.path.join(working_diamond_dir, f"{sample1}_vs_{sample2}_scaffold{config['work_config']['blastp_suffix']}")
            blastp_out2 = os.path.join(working_diamond_dir, f"{sample2}_vs_{sample1}_scaffold{config['work_config']['blastp_suffix']}")
            bed_file1 = os.path.join(working_bed_dir, f"{sample1}_scaffold.bed")
            bed_file2 = os.path.join(working_bed_dir, f"{sample2}_scaffold.bed")
            dagchainer_input1 = os.path.join(working_dagchainer_dir, f"{sample1}_vs_{sample2}_scaffold{config['work_config']['dagchainer_input_suffix']}")
            dagchainer_input2 = os.path.join(working_dagchainer_dir, f"{sample2}_vs_{sample1}_scaffold{config['work_config']['dagchainer_input_suffix']}")
            tasks_dagchainer_input.append((blastp_out1, blastp_out2, sample1, sample2, bed_file1, bed_file2, dagchainer_input1, dagchainer_input2))
    
    if not tasks_dagchainer_input:
        PipelineLogger.warning("No valid DAGchainer input tasks to run")
        PipelineLogger.end_stage("DAGchainer Analysis")
        return FileProcessor.ensure_directory(os.path.join(working_dagchainer_dir, "output"))
    
    results_input = ProcessManager.run_parallel_tasks(
        tasks_dagchainer_input,
        get_dagchainer_input,
        max_workers=thread,
        task_description="Generating DAGchainer input files"
    )
    
    successful_input = ProcessManager.is_stage_successful(results_input)
    PipelineLogger.info(f"DAGchainer input generation: {results_input['success']} success, {results_input['skipped']} skipped, {results_input['errors']} errors")
    
    # 运行DAGchainer
    PipelineLogger.info("Running DAGchainer...")
    working_dagchainer_output_dir = FileProcessor.ensure_directory(os.path.join(working_dagchainer_dir, "output"))
    
    tasks_dagchainer_run = []
    dagchainer_args = [str(item) for kv in config["DAGchainer_config"].items() for item in kv]
    
    for sample1 in sample_list:
        for sample2 in sample_list:
            if sample1 == sample2:
                continue

            # chro vs chro
            for chromosome in chromosome_list:
                inp = os.path.join(working_dagchainer_dir, f"{sample1}_vs_{sample2}_{chromosome}{config['work_config']['dagchainer_input_suffix']}")
                out = os.path.join(working_dagchainer_output_dir, f"{sample1}_vs_{sample2}_{chromosome}{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
                if os.path.exists(inp) and os.path.getsize(inp) > 0:
                    tasks_dagchainer_run.append((perl_path, dagchainer_path, inp, out, dagchainer_args))

            # scaffold vs chro
            for chromosome in chromosome_list:
                inp = os.path.join(working_dagchainer_dir, f"{sample1}_vs_{sample2}_scaffold_vs_{chromosome}{config['work_config']['dagchainer_input_suffix']}")
                out = os.path.join(working_dagchainer_output_dir, f"{sample1}_vs_{sample2}_scaffold_vs_{chromosome}{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
                if os.path.exists(inp) and os.path.getsize(inp) > 0:
                    tasks_dagchainer_run.append((perl_path, dagchainer_path, inp, out, dagchainer_args))

            # chro vs scaffold
            for chromosome in chromosome_list:
                inp = os.path.join(working_dagchainer_dir, f"{sample1}_vs_{sample2}_{chromosome}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}")
                out = os.path.join(working_dagchainer_output_dir, f"{sample1}_vs_{sample2}_{chromosome}_vs_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
                if os.path.exists(inp) and os.path.getsize(inp) > 0:
                    tasks_dagchainer_run.append((perl_path, dagchainer_path, inp, out, dagchainer_args))

            # scaffold vs scaffold
            inp = os.path.join(working_dagchainer_dir, f"{sample1}_vs_{sample2}_scaffold{config['work_config']['dagchainer_input_suffix']}")
            out = os.path.join(working_dagchainer_output_dir, f"{sample1}_vs_{sample2}_scaffold{config['work_config']['dagchainer_input_suffix']}{config['work_config']['dagchainer_output_suffix']}")
            if os.path.exists(inp) and os.path.getsize(inp) > 0:
                tasks_dagchainer_run.append((perl_path, dagchainer_path, inp, out, dagchainer_args))
    
    if not tasks_dagchainer_run:
        PipelineLogger.warning("No valid DAGchainer run tasks")
        PipelineLogger.end_stage("DAGchainer Analysis")
        return working_dagchainer_output_dir
    
    results_run = ProcessManager.run_parallel_tasks(
        tasks_dagchainer_run,
        lambda perl, dag, inp, out, args: DAGchainerManager.run_DAG_chainer_pl(perl, dag, inp, out, args),
        max_workers=thread,
        task_description="Running DAGchainer"
    )
    
    successful_run = ProcessManager.is_stage_successful(results_run)
    PipelineLogger.info(f"DAGchainer execution: {results_run['success']} success, {results_run['skipped']} skipped, {results_run['errors']} errors")
    
    PipelineLogger.end_stage("DAGchainer Analysis")
    
    return working_dagchainer_output_dir


def run_synteny_analysis(working_bed_dir: str, dagchainer_output_dir: str,
                        chromosome_list: List[str], sample_list: List[str], 
                        output_dir: str, thread: int) -> List[Tuple[str, str, List[Tuple[str, str]]]]:
    """运行同源分析，增加错误处理"""
    PipelineLogger.start_stage("Synteny Analysis")
    
    tasks = []
    for sample1, sample2 in combinations(sample_list, 2):
        tasks.append((sample1, sample2, working_bed_dir, chromosome_list, dagchainer_output_dir))
    
    results = ProcessManager.run_parallel_tasks(
        tasks,
        lambda s1, s2, bed, chro, dag: sample1_vs_sample2_pair(s1, s2, bed, chro, dag),
        max_workers=thread,
        task_description="Synteny Analysis"
    )
    
    # 过滤掉None结果
    valid_results = [r for r in results["results"] if r is not None]
    
    PipelineLogger.info(f"Synteny analysis completed: {len(valid_results)} valid pairs out of {len(tasks)} total")
    PipelineLogger.end_stage("Synteny Analysis")
    
    return valid_results


def run_mcl_analysis(all_RBNH_file_path: str, edges_list: List[List[Tuple[str, str]]], 
                    working_dir: str, thread: int) -> str:
    """运行MCL聚类分析，返回结果文件路径"""
    PipelineLogger.start_stage("MCL Clustering Analysis")
    
    mcl_path = MCLManager.get_mcl_path()
    if not mcl_path:
        PipelineLogger.error("MCL not available, skipping MCL analysis")
        PipelineLogger.end_stage("MCL Clustering Analysis")
        return ""
    
    mcl_input_dir = FileProcessor.ensure_directory(os.path.join(working_dir, config["work_config"]["mcl_input_dir_name"]))
    mcl_output_dir = FileProcessor.ensure_directory(os.path.join(working_dir, config["work_config"]["mcl_output_dir_name"]))
    
    # 检查RBNH文件是否存在
    if not os.path.exists(all_RBNH_file_path) or os.path.getsize(all_RBNH_file_path) == 0:
        PipelineLogger.warning(f"RBNH file not found or empty: {all_RBNH_file_path}")
        PipelineLogger.end_stage("MCL Clustering Analysis")
        return ""
    
    # 生成重复RBNH文件
    duplicate_gene_pair_file_list = get_duplicate_RBNH_parallel(
        all_RBNH_file_path, edges_list, mcl_input_dir, thread
    )
    
    if not duplicate_gene_pair_file_list:
        PipelineLogger.warning("No duplicate RBNH files generated")
        PipelineLogger.end_stage("MCL Clustering Analysis")
        return ""
    
    # 运行MCL聚类
    PipelineLogger.info("Running MCL clustering...")
    run_mcl_task_list = []
    mcl_args_list = [str(item) for kv in config["mcl_config"].items() for item in kv]
    
    for mcl_input_file_path in duplicate_gene_pair_file_list:
        mcl_input_file_basename = os.path.basename(mcl_input_file_path)
        mcl_input_file_basename_no_suffix = mcl_input_file_basename.split(config["work_config"]["mcl_input_suffix"])[0]
        out_file_path = os.path.join(mcl_output_dir, f"{mcl_input_file_basename_no_suffix}{config['work_config']['mcl_output_suffix']}")
        run_mcl_task_list.append((mcl_path, mcl_input_file_path, out_file_path, mcl_args_list))
    
    results_mcl = ProcessManager.run_parallel_tasks(
        run_mcl_task_list,
        lambda mcl, inp, out, args: MCLManager.run_mcl(mcl, inp, out, args),
        max_workers=thread,
        task_description="MCL Clustering"
    )
    
    successful_mcl = ProcessManager.is_stage_successful(results_mcl)
    PipelineLogger.info(f"MCL clustering: {results_mcl['success']} success, {results_mcl['skipped']} skipped, {results_mcl['errors']} errors")
    
    # 合并MCL结果
    mcl_all_result_path = os.path.join(mcl_output_dir, config["work_config"]["mcl_output_merge_file_name"])
    if config["overwrite"]["merge_mcl_output"] or not os.path.exists(mcl_all_result_path) or (os.path.exists(mcl_all_result_path) and os.path.getsize(mcl_all_result_path) == 0):
        try:
            with open(mcl_all_result_path, "w") as f:
                for this_task in run_mcl_task_list:
                    this_mcl_output_path = this_task[2]
                    if os.path.exists(this_mcl_output_path) and os.path.getsize(this_mcl_output_path) > 0:
                        with open(this_mcl_output_path) as f2:
                            for line in f2:
                                line = line.strip()
                                if line:
                                    f.write(f"{line}\n")
            PipelineLogger.success(f"Merged MCL results: {mcl_all_result_path}")
        except Exception as e:
            PipelineLogger.error(f"Error merging MCL results: {e}")
            mcl_all_result_path = ""
    else:
        PipelineLogger.skip(f"Merged MCL results already exist: {mcl_all_result_path}")
    
    PipelineLogger.end_stage("MCL Clustering Analysis")
    
    return mcl_all_result_path


def main():
    """主函数"""
    try:
        # 参数解析
        helptext = r"""
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

        3.The sample.list file (the order of the samples determines the order of the geneindex column)
        Bna_Darmor_v10.0
        Bna_Darmor_v5.0
        Bna_Express617_v1.0

        4.The chro.list file
        Chr1
        Chr2
        Chr3
        """
        
        parser = argparse.ArgumentParser(description=helptext, formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument("-i", "--input", required=True, help="输入目录")
        parser.add_argument("-o", "--output", required=True, help="输出目录")
        parser.add_argument("-s", "--samplelist", type=str, default="", help="样本列表文件")
        parser.add_argument("-c", "--chromlist", type=str, default="", help="染色体列表文件")
        parser.add_argument("-t", "--thread", type=int, default=10, help="线程数")
        parser.add_argument("-b", "--bed", type=str, choices=["bed4", "bed6"], default="bed6", help="bed文件格式选择")
        parser.add_argument("--split", action="store_true", help="是否进行最终的 GeneIndex 分割")
        args = parser.parse_args()
        
        PipelineLogger.start_pipeline()
        
        # 输入验证
        if not os.path.exists(args.input):
            PipelineLogger.error(f"Input directory does not exist: {args.input}")
            sys.exit(1)
        
        # 读取样本和染色体列表
        sample_list = FileProcessor.read_list_file(
            args.samplelist or os.path.join(args.input, config["input_file_type"]["sample_file_name"]),
            "Sample list"
        )
        chromosome_list = FileProcessor.read_list_file(
            args.chromlist or os.path.join(args.input, config["input_file_type"]["chromosome_file_name"]),
            "Chromosome list"
        )
        
        # 验证输入文件
        if not FileProcessor.validate_input_files(args.input, sample_list, 
                                                config["input_file_type"]["bed_suffix"], "BED"):
            PipelineLogger.warning("Some BED files are missing, but continuing...")
        
        if not FileProcessor.validate_input_files(args.input, sample_list,
                                                config["input_file_type"]["pep_suffix"], "PEP"):
            PipelineLogger.warning("Some PEP files are missing, but continuing...")
        
        # 创建输出目录
        output_dir = FileProcessor.ensure_directory(args.output)
        working_dir = FileProcessor.ensure_directory(os.path.join(output_dir, config["work_config"]["working_dir_name"]))
        working_bed_dir = FileProcessor.ensure_directory(os.path.join(working_dir, config["work_config"]["chromosome_bed_dir_name"]))
        working_pep_dir = FileProcessor.ensure_directory(os.path.join(working_dir, config["work_config"]["chromosome_pep_dir_name"]))
        
        # 执行各个处理阶段
        pep_bed_success = run_pep_bed_split(args.input, working_pep_dir, working_bed_dir, chromosome_list, sample_list, args.thread)
        
        if not pep_bed_success:
            PipelineLogger.error("PEP/BED splitting failed with errors, stopping pipeline")
            PipelineLogger.end_pipeline()
            return
        
        # Diamond 分析
        working_diamond_dir = FileProcessor.ensure_directory(os.path.join(working_dir, config["work_config"]["blastp_dir_name"]))
        diamond_path = DiamondManager.get_diamond_path()
        
        diamond_success = run_diamond_analysis(working_pep_dir, working_diamond_dir, chromosome_list, sample_list, diamond_path, args.thread)
        
        if not diamond_success:
            PipelineLogger.error("Diamond analysis failed with errors, stopping pipeline")
            PipelineLogger.end_pipeline()
            return
        
        # DAGchainer 分析
        working_dagchainer_dir = FileProcessor.ensure_directory(os.path.join(working_dir, config["work_config"]["dagchainer_input_dir_name"]))
        dagchainer_output_dir = run_dagchainer_analysis(working_diamond_dir, working_bed_dir, working_dagchainer_dir, 
                                                       chromosome_list, sample_list, args.thread)
        
        # 同源分析
        synteny_results = run_synteny_analysis(working_bed_dir, dagchainer_output_dir, 
                                              chromosome_list, sample_list, output_dir, args.thread)
        
        # 生成 GeneIndex
        edges_list = combine_SG(synteny_results, chromosome_list, sample_list, working_bed_dir, output_dir)
        
        # 拆分 geneindex
        if args.split:
            # RBNH计算
            RBNH_output_dir = FileProcessor.ensure_directory(os.path.join(working_dir, config["work_config"]["RBNH_output_dir_name"]))
            rbnh_success = compute_RBNH(sample_list, working_diamond_dir, RBNH_output_dir, chromosome_list, args.thread)
            if not rbnh_success:
                PipelineLogger.warning("RBNH computation had issues, but continuing...")
            # MCL聚类分析
            mcl_result_path = ""
            if rbnh_success:
                all_RBNH_file_path = os.path.join(RBNH_output_dir, config["work_config"]["all_RBNH_file_name"])
                if os.path.exists(all_RBNH_file_path) and os.path.getsize(all_RBNH_file_path) > 0:
                    mcl_result_path = run_mcl_analysis(all_RBNH_file_path, edges_list, working_dir, args.thread)
                else:
                    PipelineLogger.warning("RBNH file not available, skipping MCL analysis")
            else:
                PipelineLogger.warning("Skipping MCL analysis due to RBNH issues")
            # 最终GeneIndex分割
            if mcl_result_path and os.path.exists(mcl_result_path) and os.path.getsize(mcl_result_path) > 0:
                run_all_mcl_split(mcl_result_path, output_dir,sample_list)
            else:
                PipelineLogger.warning("MCL results not available, skipping GeneIndex splitting")
        PipelineLogger.end_pipeline()
    except Exception as e:
        PipelineLogger.error(f"Pipeline failed: {e}")
        PipelineLogger.error(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()
