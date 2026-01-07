#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: ordopan.py
# Created time: 2025/12/26 17:21:56
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python ordopan.py -h


import argparse
import subprocess
import sys
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent

SCRIPT_MAP = {
    # main tools
    "synpan": BASE_DIR / "src" / "synpan.py",
    "ordo": BASE_DIR / "src" / "ordo.py",

    # sub tools
    "cluster": BASE_DIR / "src" / "cluster.py",
    "sort": BASE_DIR / "src" / "Destruction-Fusion-Cycle.py",
    "resort": BASE_DIR / "src" / "reorder_geneindex_by_sorted_num.py",

    # utils
    "bed2wgdi": BASE_DIR / "src" / "utils" / "bed2wgdi.py",
    "kendall": BASE_DIR / "src" / "utils" / "kendallTau.py",
    "gfa": BASE_DIR / "src" / "utils" / "geneindex2gfa.py",
    "score": BASE_DIR / "src" / "utils" / "num2adj-loss.py",
    "old": BASE_DIR / "src" / "utils" / "geneindex2old.py",
    "pep2bed": BASE_DIR / "src" / "utils" / "pep2bed.py",
    "imitate": BASE_DIR / "src" / "utils" / "imitate.py",
}


def build_command(script: Path, args: list[str]) -> list[str]:
    """Python 脚本用当前解释器，其它直接执行"""
    if script.suffix == ".py":
        return [sys.executable, str(script)] + args
    else:
        return [str(script)] + args


def main():
    argv = sys.argv[1:]

    # ========== 手动处理顶层 -h / -v ==========
    if not argv or argv == ["-h"] or argv == ["--help"]:
        print("usage: ordopan <tool> [args...]")
        print("available tools:\n"+",".join(SCRIPT_MAP.keys()))
        sys.exit(0)

    if argv == ["-v"] or argv == ["--version"]:
        print("ordopan version: 1.0")
        sys.exit(0)

    # ========== argparse 只做“结构解析”，无 help ==========
    parser = argparse.ArgumentParser(prog="ordopan",add_help=False)
    parser.add_argument("tool", choices=SCRIPT_MAP.keys())
    parser.add_argument("args", nargs=argparse.REMAINDER)

    args = parser.parse_args(argv)

    script = SCRIPT_MAP[args.tool]
    cmd = build_command(script, args.args)
    subprocess.run(cmd)

if __name__ == "__main__":
    main()
