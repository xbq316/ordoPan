#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: sort_geneindex.py
# Created time: 2025/07/01 19:02:34
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python sort_geneindex.py -h



import multiprocessing
import os
import subprocess
from rich.console import Console
import argparse

console = Console()


def check_scipt(script_path):
    if not os.path.exists(script_path):
        console.print(f"{script} does not exist",style="red")
        exit(1)
        
def check_success(result,desc:str=""):
    sep="*"
    console.print(sep*50)
    console.print(" ".join(result.args),style="green")
    if result.returncode != 0:
        console.print(f"Error: {result.stderr}",style="red")
        exit(1)
    else:
        # 输出命令行
        console.print(desc,style="green")
    console.print(sep*50)
    console.print()


def process_chro(chro):
    # 聚类
    numfile=os.path.join(output_dir,f"{chro}.num")
    cluster_file=os.path.join(output_dir,f"{chro}.cluster")
    run_cluster_cmd =["python3",cluster_script,"-i",numfile,"-o",cluster_file]
    cluster_result = subprocess.run(run_cluster_cmd,text=True)
    check_success(cluster_result,f"{cluster_script} 运行成功. {chro}.cluster 文件已生成")

    # 对cluster排序
    if args.cycle:
        # 运行 Destruction-Fusion-Cycle.py 脚本进行 破坏-融合
        sort_num_file_prefix = os.path.join(output_dir, f"{chro}")
        fusion_cmd = ["python3", Destruction_Fusion_Cycle_script, "-i", cluster_file, "-o", sort_num_file_prefix, "-r"] + list(map(str,args.rate)) + ["-n", str(args.max_iter), "-m", str(args.minimum_genome_length)]
        fusion_result = subprocess.run(fusion_cmd, text=True)
        check_success(fusion_result,f"{Destruction_Fusion_Cycle_script} 运行成功. {chro}.new.cluster {chro}.sorted.num 文件已生成")
    else:
        sort_num_file=os.path.join(output_dir,f"{chro}.sorted.num")
        run_sort_cmd =[sort_exec,"-i",cluster_file,"-o",sort_num_file,"-n",str(iteration)]
        sort_result = subprocess.run(run_sort_cmd,text=True)
        check_success(sort_result,f"{sort_exec} 运行成功. {chro}.sorted.num 文件已生成")

    # 还原排序路径
    path_file=os.path.join(output_dir,f"{chro}.path")
    run_path_cmd =["python3",path_script,"-n",numfile,"-s",sort_num_file,"-o",path_file]
    path_result = subprocess.run(run_path_cmd,text=True)
    check_success(path_result,f"{path_script} 运行成功. {chro}.path 文件已生成")

    # 将 annotation 按照 排序后的 num 进行排序
    annotation_file=os.path.join(output_dir,f"{chro}.annotation.txt")
    annotation_sort_file=os.path.join(output_dir,f"{chro}.sort.annotation.txt")
    resort_cmd =["python3",resort_script,"-i",annotation_file,"-o",annotation_sort_file,"-p",path_file]
    resort_result = subprocess.run(resort_cmd,text=True)
    check_success(resort_result,f"{resort_script} 运行成功. {annotation_sort_file} 文件已生成")

    # 绘制 html
    html_file=os.path.join(output_dir,f"{chro}.html")
    html_cmd =["python3",html_script,"-i",sort_num_file,"-a",annotation_sort_file,"-o",html_file]
    html_result = subprocess.run(html_cmd,text=True)
    check_success(html_result,f"{html_script} 运行成功. {html_file} 文件已生成")

if __name__ =="__main__":

    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)

    parser=argparse.ArgumentParser()
    parser.add_argument("-i","--input",required=True,type=str,help="输入文件目录")
    parser.add_argument("-t","--thread",default=1,type=int,help="线程数")
    parser.add_argument("-n","--iteration",default=10,type=int,help="迭代次数(单位:万次)")
    parser.add_argument("-c","--chro",nargs='+',type=str,help="染色体列表")

    parser.add_argument("--cycle",action='store_true',help="是否调用 Destruction-Fusion-Cycle.py 脚本进行 破坏-融合 ,默认调用")
    parser.add_argument("-r","--rate",nargs='+',type=float,default=[0.5,0.6,0.7,0.8],help="融合阈值(默认: [0.5,0.6,0.7,0.8]),越大越严格，则融合越少")
    parser.add_argument("-m", "--minimum-genome-length", default=100,help="最小基因组长度(用于打分使用，避免有的样本基因个数太少，导致打分得到NA或是不准确)")
    parser.add_argument("--max-iter", type=int, default=5, help="每个阈值最大迭代次数 (默认: 5)")
    args=parser.parse_args()

    # 脚本所在目录
    output_dir = args.input
    THREAD=args.thread
    chros=args.chro
    iteration=args.iteration


    script="cluster"
    cluster_script = os.path.join(script_dir, script)
    check_scipt(cluster_script)

    script="sort2cluster.py"
    sort2cluster_script = os.path.join(script_dir, script)
    check_scipt(sort2cluster_script)

    script="Destruction-Fusion-Cycle.py"
    Destruction_Fusion_Cycle_script = os.path.join(script_dir, script)
    check_scipt(Destruction_Fusion_Cycle_script)

    script="clusterSort"
    sort_exec = os.path.join(script_dir, "clusterSort")
    check_scipt(sort_exec)

    script="num_sortnum2path.py"
    path_script = os.path.join(script_dir, script)
    check_scipt(path_script)

    script="resort.py"
    resort_script = os.path.join(script_dir, script)
    check_scipt(resort_script)

    script="draw_html.py"
    html_script = os.path.join(script_dir, script)
    check_scipt(html_script)

    if len(chros)<=THREAD:
        THREAD=len(chros)
    with multiprocessing.Pool(THREAD) as pool:
            pool.map(process_chro, chros)
