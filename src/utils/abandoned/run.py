#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: run.py
# Created time: 2025/06/25 20:36:58
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python run.py -h


import shutil
import subprocess
import os
from rich.console import Console
import argparse
import multiprocessing

helptext=r"""
输入文件要求: 
1.氨基酸序列文件(.prot): file1.prot file2.prot file3.prot ...
>gene1
MKTAYIAKQRQISFVKSHFSRQDILDLI
>gene2
GKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSRHPGNFGAF

2.bed 文件： file1.bed file2.bed file3.bed ...
{chromosome}\t{start}\t{end}\t{genename}
A01	831	1437	BnaA01g00010D
A01	1487	2436	BnaA01g00020D
A01	2665	5455	BnaA01g00030D
A01	8421	9623	BnaA01g00040D

3.sample.list 文件(样本的顺序决定了 geneindex 列的顺序)
Bna_Darmor_v10.0
Bna_Darmor_v5.0
Bna_Express617_v1.0

4.chro.list 文件必须存在(常染色体文件,非常染色体将只会被统一编为 Scaffold)
A01
A02
A03
A04
A05
A06
...

"""


parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i","--input",required=True,type=str,help="输入文件目录")
parser.add_argument("-o","--output",required=True,type=str,help="输出文件目录")
parser.add_argument("-n","--iteration",default=10,type=int,help="迭代次数(单位:万次)")
parser.add_argument("-w","--overwrite",action="store_true",help="是否覆盖已存在的blast结果文件(默认不覆盖)")
parser.add_argument("-t","--thread",type=int,default=10,help="线程数")

parser.add_argument("--cycle",action='store_true',help="是否调用 Destruction-Fusion-Cycle.py 脚本进行 破坏-融合 ,默认调用")
parser.add_argument("-r","--rate",nargs='+',type=float,default=[0.5,0.6,0.7,0.8],help="融合阈值(默认: [0.5,0.6,0.7,0.8]),越大越严格，则融合越少")
parser.add_argument("-m", "--minimum-genome-length", default=100,help="最小基因组长度(用于打分使用，避免有的样本基因个数太少，导致打分得到NA或是不准确)")
parser.add_argument("--max-iter", type=int, default=5, help="每个阈值最大迭代次数 (默认: 5)")

args=parser.parse_args()
console = Console()


input_dir=args.input
output_dir=args.output
overWrite=args.overwrite
THREAD=args.thread
iteration=args.iteration

# 脚本所在目录
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

# 检查 input_dir 是否存在
if not os.path.exists(input_dir):
    console.print(f"input directory {input_dir} does not exist",style="yellow")
    exit(1)
else:
    console.print(f"input directory {input_dir} exists",style="green")

# 检查 output_dir 是否存在
if not os.path.exists(output_dir):
    console.print(f"output directory {output_dir} does not exist",style="yellow")
    console.print(f"creating output directory {output_dir}",style="yellow")
    os.mkdir(output_dir)
else:
    console.print(f"output directory {output_dir} exists",style="green")
    # 检查是否有写权限
    if not os.access(output_dir, os.W_OK):
        console.print(f"[快速检查] 无写权限: {dir_path}",style="red")

# 开始检测输入目录中的文件情况： file.bed file.prot chro.list sample.list
# 1.检测 sample.list
if not os.path.exists(os.path.join(input_dir,"sample.list")):
    console.print(f"sample.list does not exist",style="yellow")
    # creating sample.list
    console.print(f"creating sample.list",style="yellow")
    
    bedset=set()
    protset=set()
    for file in os.listdir(input_dir):
        if file.endswith(".bed"):
            bedset.add(file.split(".bed")[0])
        elif file.endswith(".prot"):
            protset.add(file.split(".prot")[0])
    if bedset==protset:
        sample_list=open(os.path.join(input_dir,"sample.list"),"w")
        for file in bedset:
            sample_list.write(f"{file}\n")
        sample_list.close()
    else:
        console.print(f"sample.list 创建失败",style="red")
        a=bedset-protset
        for file in a:
            console.print(f"存在 {file}.bed 但不存在 {file}.prot",style="red")
        b=protset-bedset
        for file in b:
            console.print(f"存在 {file}.prot 但不存在 {file}.bed",style="red")
        exit(1)
else:
    console.print(f"sample.list exists",style="green")

# 2.检测 chro.list
if not os.path.exists(os.path.join(input_dir,"chro.list")):
    console.print(f"chro.list(常染色体文件) does not exist",style="red")
    exit(1)
else:
    console.print(f"chro.list exists",style="green")


# 是否覆盖已存在的 blast 文件
if overWrite:
    console.print(f"若是blast结果存在，则覆盖",style="yellow")
else:
    console.print(f"若是blast结果存在，则跳过",style="yellow")


#################### 检测 脚本是否存在 ####################

# run_DAG_chainer.pl  dagchainer
if not os.path.exists(os.path.join(script_dir,"run_DAG_chainer.pl")):
    console.print(f"run_DAG_chainer.pl does not exist",style="red")
    exit(1)
if not os.path.exists(os.path.join(script_dir,"dagchainer")):
    console.print(f"dagchainer does not exist",style="red")
    exit(1)

# 检查 DIAMOND 可执行文件是否存在于 PATH 中
# 如果不存在，提示用户手动安装
DIAMOND = "diamond"
diamond_path = shutil.which(DIAMOND)
if not diamond_path:
    console.print(f"{DIAMOND} not found in PATH",style="red")
    console.print(f"please install {DIAMOND} and add it to PATH",style="red")
    exit(1)

# 检测 perl 是否存在
PERL = "perl"
perl_path = shutil.which(PERL)
if not perl_path:
    console.print(f"{PERL} not found in PATH",style="red")
    console.print(f"please install {PERL} and add it to PATH",style="red")
    exit(1)

def check_scipt(script_path):
    if not os.path.exists(script_path):
        console.print(f"{script} does not exist",style="red")
        exit(1)

script="synpan.pl"
script_path = os.path.join(script_dir, script)
check_scipt(script_path)

script="synpan2geneindex.py"
script_path = os.path.join(script_dir, script)
check_scipt(script_path)

script="annotation.py"
script_path = os.path.join(script_dir, script)
check_scipt(script_path)

script="vote.py"
script_path = os.path.join(script_dir, script)
check_scipt(script_path)

script="start.py"
script_path = os.path.join(script_dir, script)
check_scipt(script_path)

script="split.py"
script_path = os.path.join(script_dir, script)
check_scipt(script_path)

script="reindex.py"
script_path = os.path.join(script_dir, script)
check_scipt(script_path)

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

#################### 检测 脚本是否存在 ####################

# 检测运行状态
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



# 运行 synpan
script="synpan.pl"
script_path = os.path.join(script_dir, script)
# 检测最终文件是否存在，考虑是否覆盖
if os.path.exists(os.path.join(output_dir,"sample.list.SG.pan")):
    console.print(f"sample.list.SG.pan 文件已存在,跳过 synpan.pl",style="green")
else:
    synpan_result=None
    if overWrite:
        run_synpan_cmd =["perl",script_path,"sample.list",output_dir,str(THREAD),"1"]
        synpan_result = subprocess.run(run_synpan_cmd,cwd=input_dir,text=True)
    else:
        run_synpan_cmd =["perl",script_path,"sample.list",output_dir,str(THREAD),"0"]
        synpan_result = subprocess.run(run_synpan_cmd,cwd=input_dir,text=True)
    if synpan_result!=None:
        check_success(synpan_result,"synpan 运行成功")

#获得 geneindex
script="synpan2geneindex.py"
script_path = os.path.join(script_dir, script)
synpan_file = os.path.join(output_dir, "sample.list.SG.pan")
geneindex_file = os.path.join(output_dir, "geneindex.txt")
sample_list_file=os.path.join(input_dir, "sample.list")
run_geneindex_cmd =["python3",script_path,"-l",sample_list_file,"-s",synpan_file,"-o",geneindex_file]
geneindex_result = subprocess.run(run_geneindex_cmd,text=True)
check_success(geneindex_result,f"{script} 运行成功. {geneindex_file} 文件已生成")

# 获取染色体注释文件/获取注释文件
script="annotation.py"
script_path = os.path.join(script_dir, script)
run_annotation_cmd =["python3",script_path,"-b",input_dir,"-i",geneindex_file,"-o",output_dir]
annotation_result = subprocess.run(run_annotation_cmd,text=True)
check_success(annotation_result,f"{script} 运行成功. geneindex.chro.txt 和 geneindex.annotation.txt 文件已生成")
geneindex_chro=os.path.join(output_dir,"geneindex.chro.txt")
geneindex_annotation=os.path.join(output_dir,"geneindex.annotation.txt")

# 获取染色体投票文件
script="vote.py"
script_path = os.path.join(script_dir, script)
chro_vote_file=os.path.join(output_dir,"chro.vote.txt")
vote_prefix = os.path.join(output_dir, "sample")
run_vote_cmd =["python3",script_path,"-i",geneindex_chro,"-o",chro_vote_file]
vote_result = subprocess.run(run_vote_cmd,text=True)
check_success(vote_result,f"{script} 运行成功. {chro_vote_file} 文件已生成")

# 获取每个基因的位置(start位置)
script="start.py"
script_path = os.path.join(script_dir, script)
start_file=os.path.join(output_dir,"start.txt")
run_start_cmd = ["python3",script_path,"-b",input_dir,"-i",geneindex_file,"-o",start_file,"-v",chro_vote_file]
start_result = subprocess.run(run_start_cmd,text=True)
check_success(start_result,f"{script} 运行成功. {start_file} 文件已生成")

# 按照染色体划分为多个 start 文件
script="split.py"
script_path = os.path.join(script_dir, script)
chro_file=os.path.join(input_dir,"chro.list")#常染色体文件
run_split_cmd =["python3",script_path,
    "-i",start_file,
    "-o",output_dir,
    "-v",chro_vote_file,
    "-s",".start.txt",
    "-c",chro_file]
split_result = subprocess.run(run_split_cmd,text=True)
check_success(split_result,f"{script} 运行成功. start 文件 已划分为多个文件")

# 按照染色体划分为多个 annotation 文件
run_split_cmd =["python3",script_path,
    "-i",geneindex_annotation,
    "-o",output_dir,
    "-s",".annotation.txt",
    "-v",chro_vote_file,
    "-c",chro_file]
split_result = subprocess.run(run_split_cmd,text=True)
check_success(split_result,f"{script} 运行成功. 注释文件 已划分为多个文件")

# 常染色体文件
chros=[]
for chro in open(chro_file):
    chro=chro.strip()
    chros.append(chro)
chros.sort()

# 重新打上 索引标签
script="reindex.py"
script_path = os.path.join(script_dir, script)
for chro in chros:
    start_file = os.path.join(output_dir,f"{chro}.start.txt")
    num_file=os.path.join(output_dir,f"{chro}.num")
    run_reindex_cmd =["python3",script_path,"-i",start_file,"-o",num_file]
    reindex_result = subprocess.run(run_reindex_cmd,text=True)
    check_success(reindex_result,f"reindex.py 运行成功. {chro}.num 文件已生成")


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
        sort_num_file=os.path.join(output_dir,f"{chro}.sorted.num")
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

with multiprocessing.Pool(THREAD) as pool:
        pool.map(process_chro, chros)
