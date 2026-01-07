#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: runPan.py
# Created time: 2025/06/25 20:36:58
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python runPan.py -h


import shutil
import subprocess
import os
from rich.console import Console
import argparse

helptext=r"""
Input file requirements: 
1.amino acid sequence file(.prot): file1.prot file2.prot file3.prot ...
>gene1
MKTAYIAKQRQISFVKSHFSRQDILDLI
>gene2
GKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSRHPGNFGAF


2.bed file: file1.bed file2.bed file3.bed ...
{chromosome}\t{start}\t{end}\t{genename}
A01	831	1437	BnaA01g00010D
A01	1487	2436	BnaA01g00020D
A01	2665	5455	BnaA01g00030D
A01	8421	9623	BnaA01g00040D

3.The sample.list file (the order of the samples determines the order of the geneindex column)
Bna_Darmor_v10.0
Bna_Darmor_v5.0
Bna_Express617_v1.0
...
"""


parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i","--input",required=True,type=str,help="Input file directory")
parser.add_argument("-o","--output",required=True,type=str,help="Output file directory")

parser.add_argument("-w","--overwrite",action="store_true",help="Whether to overwrite existing BLAST result files (default is no overwrite, meaning if a BLAST result file exists, the BLAST step will be skipped).")
parser.add_argument("-t","--thread",type=int,default=20,help="diamond blastp thread count")
args=parser.parse_args()

console = Console()


input_dir=os.path.abspath(args.input)
output_dir=os.path.abspath(args.output)
overWrite=args.overwrite
BLASTTHREAD=args.thread

# Script directory
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

#################### Check if the script exists ####################

# run_DAG_chainer.pl  dagchainer
if not os.path.exists(os.path.join(script_dir,"run_DAG_chainer.pl")):
    console.print(f"run_DAG_chainer.pl does not exist",style="red")
    exit(1)
if not os.path.exists(os.path.join(script_dir,"dagchainer")):
    console.print(f"dagchainer does not exist",style="red")
    exit(1)

# Check if the DIAMOND executable file exists in the PATH
# If it does not exist, prompt the user to install it manually.DIAMOND = "diamond"
diamond_path = shutil.which(DIAMOND)
if not diamond_path:
    console.print(f"{DIAMOND} not found in PATH",style="red")
    console.print(f"please install {DIAMOND} and add it to PATH",style="red")
    exit(1)

# Check if perl exists
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
#################### Check if the script exists ####################

# Check if input_dir exists.
if not os.path.exists(input_dir):
    console.print(f"input directory {input_dir} does not exist",style="yellow")
    exit(1)
else:
    console.print(f"input directory {input_dir} exists",style="green")

# Check if output_dir exists.
if not os.path.exists(output_dir):
    console.print(f"output directory {output_dir} does not exist",style="yellow")
    console.print(f"creating output directory {output_dir}",style="yellow")
    os.mkdir(output_dir)
else:
    console.print(f"output directory {output_dir} exists",style="green")
    # Check for write permissions.
    if not os.access(output_dir, os.W_OK):
        console.print(f"[快速检查] 无写权限: {output_dir}",style="red")


# Start checking the files in the input directory: file.bed file.prot chro.list sample.list
# 1.Detection sample.list
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
    
# 是否Overwrite existing blast files
if overWrite:
    console.print(f"If the blast results exist, then overwrite.",style="yellow")
else:
    console.print(f"If the blast result exists, skip it.",style="yellow")



# Detect running status
def check_success(result,desc:str=""):
    sep="*"
    console.print(sep*50)
    console.print(" ".join(result.args),style="green")
    if result.returncode != 0:
        console.print(f"Error: {result.stderr}",style="red")
        exit(1)
    else:
        # Output command line
        console.print(desc,style="green")
    console.print(sep*50)
    console.print()


# run synpan
script="synpan.pl"
script_path = os.path.join(script_dir, script)
# Check if the final file exists, and consider whether to overwrite it.
if os.path.exists(os.path.join(output_dir,"sample.list.SG.pan")):
    console.print(f"The file sample.list.SG.pan already exists; skip synpan.pl.",style="green")
else:
    synpan_result=None
    if overWrite:
        run_synpan_cmd =["perl",script_path,"sample.list",output_dir,str(BLASTTHREAD),"1"]
        synpan_result = subprocess.run(run_synpan_cmd,cwd=input_dir,text=True)
    else:
        run_synpan_cmd =["perl",script_path,"sample.list",output_dir,str(BLASTTHREAD),"0"]
        synpan_result = subprocess.run(run_synpan_cmd,cwd=input_dir,text=True)
    if synpan_result!=None:
        check_success(synpan_result,"synpan Run successfully")

#获得 geneindex
script="synpan2geneindex.py"
script_path = os.path.join(script_dir, script)
synpan_file = os.path.join(output_dir, "sample.list.SG.pan")
geneindex_file = os.path.join(output_dir, "geneindex.txt")
sample_list_file=os.path.join(input_dir, "sample.list")
run_geneindex_cmd =["python3",script_path,"-l",sample_list_file,"-s",synpan_file,"-o",geneindex_file]
geneindex_result = subprocess.run(run_geneindex_cmd,text=True)
check_success(geneindex_result,f"{script} run successfully. {geneindex_file} The file has been generated.")
print("Done.")