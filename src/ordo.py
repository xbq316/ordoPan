#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: OrdoPan.py
# Created time: 2025/09/23 10:53:59
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python OrdoPan.py -h



from rich.console import Console
import os
import subprocess
from multiprocessing.dummy import Pool
import time

console = Console()
NONE_GENE_SET=set([".","-",""])
SEP="\t"

# Adjust the column order of the bed file
CHRO_INDEX=0
START_INDEX=1
END_INDEX=2
GENE_INDEX=3

helptext=r"""
OrdoPan: A Linear Coordinate System for Pangenome Navigation

bed : Directory containing BED files for each sample(sample1.bed,sample2.bed,sample3.bed)
A01	2982	7224	Bnir.A01G000100	.	+
A01	8588	9545	Bnir.A01G000200	.	+
A01	9174	9584	Bnir.A01G000300	.	-

sample : File containing the list of sample names(The order of the samples must match the order of samples in the geneindex file.)(sample.list)
sample1
sample2
sample3

chromosome : File containing the list of chromosomes(chro.list)
A01
A02
A03

geneindex : Gene index file(geneindex.old.txt)
A01p00010.1_BnaDAR	-	-		-	-	-	-	-	-	-	-	-	-	-
A01p00020.1_BnaDAR	A01p000050.1_BnaEXP	BnaA01T0000700GG	-	A01g000006	-	-	BnaA01T0000800SL	-	-	BnaA01T0001100ZY	BnaA01T0001800ZS	ZS11A01G000060	ZY821000042	-
A01p00030.1_BnaDAR	A01p000060.1_BnaEXP	BnaA01T0000800GG	-	A01g000007	-	-	BnaA01T0000900SL	-	-	BnaA01T0001200ZY	BnaA01T0001900ZS	ZS11A01G000070	ZY821000044	-
A01p00034.1_BnaDAR	A01p000070.1_BnaEXP	BnaA01T0000900GG,BnaA01T0001000GG	-	A01g000008	-	-	BnaA01T0001000SL	-	-	BnaA01T0001300ZY	BnaA01T0002000ZS	ZS11A01G000080	ZY821000045	-
"""

def check_sample_file_exist(sample_file)->None:
    if not os.path.exists(sample_file):
        console.print(f"{sample_file} does not exist", style="red")
        exit(1)
    else:
        console.print(f"{sample_file} exists", style="green")

def check_chro_exist(chro_list)->None:
    if not os.path.exists(chro_list):
        console.print(f"{chro_list} does not exist", style="red")
        exit(1)
    else:
        console.print(f"{chro_list} exists", style="green")

def check_bed_exist(sample_file,bed_dir)->None:
    with open(sample_file) as f:
        samples = [line.strip() for line in f if line.strip()]
    for sample in samples:
        bed_file = os.path.join(bed_dir, f"{sample}.bed")
        if not os.path.exists(bed_file):
            console.print(f"Error: {bed_file} does not exist", style="red")
            exit(1)
        else:
            console.print(f"{bed_file} exists", style="green")

def split_chro(chro_list_file:str,sample_list_file:str,geneindex_file:str,bed_dir:str,output_dir:str,overWrite: bool=False)->None:
    with open(chro_list_file) as f:
        chros = [line.strip() for line in f if line.strip()]
    chroSet=set(chros)

    if overWrite:
        # 运行：覆盖已存在的文件
        pass
    else:
        # 若是存在,则跳过
        # 若是不存在,则运行
        flag=False
        for chro in chros:
            outfile1=os.path.join(output_dir,"chromosome",f"{chro}.num")
            outfile2=os.path.join(output_dir,"chromosome",f"{chro}.txt")
            if not os.path.exists(outfile1) or not os.path.exists(outfile2):
                console.print(f"{outfile1} or {outfile2} does not exist", style="green")
                flag=True
                break
        if not flag:
            console.print(f"{os.path.join(output_dir,'chromosome')}/*.num and *.txt already exist, skip split_chro step", style="green")
            return

    sample_list=[]
    with open(sample_list_file) as f:
        for line in f:
            sample=line.strip()
            if sample:
                sample_list.append(sample)
    
    geneindex=[]
    for i in open(geneindex_file):
        line=i.strip("\n")
        if line:
            geneindex.append(line.split(SEP))

    geneindex_T=list(zip(*geneindex))
    
    console.print()
    console.print(f"Get the best chromosome...", style="white")
    sample_chros=[]
    for sample_index,sample in enumerate(sample_list):
        bed_file=os.path.join(bed_dir,f"{sample}.bed")
        gene2chro={}
        console.print(f"parse {sample}", style="green")
        with open(bed_file) as f:
            for line in f:
                parts=line.strip().split("\t")
                gene=parts[GENE_INDEX]
                chro=parts[CHRO_INDEX]
                gene2chro[gene]=chro
        this_genelist=geneindex_T[sample_index]

        chro_list=[]
        for genes in this_genelist:
            if genes in NONE_GENE_SET:
                chro_list.append("-")
            elif "," in genes:
                genes_split=[g.strip() for g in genes.split(",")]
                a_chro=[]
                for gene in genes_split:
                    if gene in gene2chro:
                        this_chro=gene2chro[gene]
                        if this_chro in chroSet:
                            a_chro.append(this_chro)
                    else:
                        raise ValueError(f"{gene} not in {bed_file}")
                if a_chro:
                    chro_list.append(a_chro)
                else:
                    chro_list.append("-")
            else:
                if genes in gene2chro:
                    this_chro=gene2chro[genes]
                    if this_chro in chroSet:
                        chro_list.append(this_chro)
                    else:
                        chro_list.append("-")
                else:
                    raise ValueError(f"{genes} not in {bed_file}")

        sample_chros.append(chro_list)
    sample_chros_T=list(zip(*sample_chros))

    max_chros=[]
    for SG_chro_list in sample_chros_T:
        a=[]
        for i in SG_chro_list:
            if isinstance(i,list):
                a.extend(i)
            elif i!="-":
                a.append(i)
        if not a:
            max_chro="randomScaffold"
        else:
            max_chro=max(a,key=a.count)
        max_chros.append(max_chro)
    console.print()

    console.print(f"Get the location of genes in their respective chromosomes...", style="white")
    sample_starts=[]
    for sample_index,sample in enumerate(sample_list):
        bed_file=os.path.join(bed_dir,f"{sample}.bed")
        gene2start={}
        gene2chro={}
        console.print(f"parse {sample}", style="green")
        with open(bed_file) as f:
            for line in f:
                parts=line.strip().split("\t") 
                chro=parts[CHRO_INDEX]
                gene=parts[GENE_INDEX]
                start=parts[START_INDEX]
                gene2start[gene]=int(start)
                gene2chro[gene]=chro
        this_genelist=geneindex_T[sample_index]

        starts=[]
        for genesindex,genes in enumerate(this_genelist):
            this_bestchro=max_chros[genesindex]
            if genes in NONE_GENE_SET:
                starts.append("-")
            elif "," in genes:
                genes_split=[g.strip() for g in genes.split(",")]
                a_start=[]
                for gene in genes_split:
                    if gene in gene2chro and gene in gene2start:
                        this_chro=gene2chro[gene]
                        if this_chro==this_bestchro:
                            a_start.append(gene2start[gene])
                    else:
                        raise ValueError(f"{gene} not in {bed_file}")
                if a_start:
                    starts.append(min(a_start))
                else:
                    starts.append("-")
            else:
                if genes in gene2chro and genes in gene2start:
                    this_chro=gene2chro[genes]
                    if this_chro==this_bestchro:
                        starts.append(gene2start[genes])
                    else:
                        starts.append("-")
                else:
                    raise ValueError(f"{genes} not in {bed_file}")
        sample_starts.append(starts)
    sample_starts_T=list(zip(*sample_starts))
    console.print()

    # Add an index (sorting label)
    # If it is "-", then no label is added.
    console.print(f"Add index for genes...", style="white")
    chro2index={chro:index for index,chro in enumerate(chros)}
    chro2index["randomScaffold"]=len(chros)
    indexs=[[] for _ in range(len(chros))]

    for genesindex,best_chro in enumerate(max_chros):
        if best_chro=="randomScaffold":
            continue
        this_index=chro2index[best_chro]
        this_start=sample_starts_T[genesindex]
        indexs[this_index].append(this_start)

    x=[]
    for chro_index,chro in enumerate(chros):
        console.print(f"Get the sequence number of the gene in each chromosome:{chro}", style="green")
        if not starts:
            continue
        starts_T=list(zip(*indexs[chro_index]))

        z=[]
        for sample_index,sample in enumerate(sample_list):
            this_start=starts_T[sample_index]
            b=[]
            c=set()
            for genesindex,start in enumerate(this_start):
                if start=="-":
                    c.add(genesindex)
                else:
                    b.append((genesindex,start))
            b.sort(key=lambda x:x[1])
            genesindex2index={x[0]:i for i,x in enumerate(b)}
            result=[]
            for genesindex in range(len(this_start)):
                if genesindex in c:
                    result.append("-")
                else:
                    result.append(genesindex2index[genesindex])
            z.append(result)
        z_T=list(zip(*z))
        x.append(z_T)
    
    # Write in chro order
    console.print(f"Write the result to file...", style="white")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(os.path.join(output_dir,"chromosome")):
        os.makedirs(os.path.join(output_dir,"chromosome"))
    
    for chro in chros:
        with open(os.path.join(output_dir,"chromosome",f"{chro}.num"),"w") as f:
            chro_index=chro2index[chro]
            this_chro_indexs=x[chro_index]
            for i in this_chro_indexs:
                i=list(i)
                f.write("\t".join(map(str,i))+"\n")

    genes=[[] for _ in range(len(chros)+1)]
    for genesindex,best_chro in enumerate(max_chros):
        chro_index=chro2index[best_chro]
        genes[chro_index].append(geneindex[genesindex])
    
    for chro in chro2index.keys():
        chro_index=chro2index[chro]
        with open(os.path.join(output_dir,"chromosome",f"{chro}.txt"),"w") as f:
            for i in genes[chro_index]:
                f.write("\t".join(i)+"\n")


def sort_geneindex(chro_list_file:str,output_dir:str,thread:int=None,overWrite:bool=False)->None:
    with open(chro_list_file) as f:
        chros = [line.strip() for line in f if line.strip()]
    if overWrite:
        pass
    else:
        # Check if sorting results already exist.
        flag=False
        for chro in chros:
            chro_file_1=os.path.join(output_dir,"cluster",f"{chro}.cluster")
            chro_file_2=os.path.join(output_dir,"sort",f"{chro}.sorted.num")
            if os.path.exists(chro_file_1) or os.path.exists(chro_file_2):
                console.print(f"{chro_file_1} or {chro_file_2} already exist", style="green")
                flag=True
                break
        if flag:
            console.print()
            console.print(f'{os.path.join(output_dir,"cluster",f"*.cluster")} and {os.path.join(output_dir,"sort","*.sorted.num")} exist', style="green")
            console.print("Skip sort_geneindex step", style="green")
            return

    # Check if the input file exists.
    for chro in chros:
        chro_file_1=os.path.join(output_dir,"chromosome",f"{chro}.num")
        chro_file_2=os.path.join(output_dir,"chromosome",f"{chro}.txt")
        
        flag=False
        if not os.path.exists(chro_file_1):
            console.print(f"{chro_file_1} does not exist", style="red")
            flag=True
        if not os.path.exists(chro_file_2):
            console.print(f"{chro_file_2} does not exist", style="red")
            flag=True
        if flag:
            console.print("Please run split_chro first", style="red")
            raise FileNotFoundError(f"{chro_file_1} or {chro_file_2} does not exist")
    
    # Check if the cluster exists
    cluster_exe=os.path.join(os.path.dirname(__file__),"cluster")
    if not os.path.exists(cluster_exe):
        cluster_cpp=os.path.join(os.path.dirname(__file__),"cluster.cpp")
        if not os.path.exists(cluster_cpp):
            console.print(f"{cluster_cpp} does not exist", style="red")
            exit(1)
        else:
            console.print(f"Compile {cluster_cpp} first", style="red")
            console.print(f"g++ {cluster_cpp} -o {cluster_exe} -O3 -std=c++17", style="red")
            exit(1)
    
    # Check if clusterSort exists
    clusterSort_exe=os.path.join(os.path.dirname(__file__),"clusterSort")
    if not os.path.exists(clusterSort_exe):
        clusterSort_cpp=os.path.join(os.path.dirname(__file__),"clusterSort.cpp")
        if not os.path.exists(clusterSort_cpp):
            console.print(f"{clusterSort_cpp} does not exist", style="red")
            raise FileNotFoundError(f"{clusterSort_exe} does not exist")
        else:
            console.print(f"Compile {clusterSort_cpp} first", style="red")
            console.print(f"g++ {clusterSort_cpp} -o {clusterSort_exe} -O3 -std=c++17", style="red")
            raise FileNotFoundError(f"{clusterSort_exe} does not exist")

    # Check if Destruction-Fusion-Cycle.py exists.
    Destruction_Fusion_Cycle_script=os.path.join(os.path.dirname(__file__),"Destruction-Fusion-Cycle.py")
    if not os.path.exists(Destruction_Fusion_Cycle_script):
        console.print(f"{Destruction_Fusion_Cycle_script} does not exist", style="red")
        raise FileNotFoundError(f"{Destruction_Fusion_Cycle_script} does not exist")

    # Create output directory
    cluster_output_dir=os.path.join(output_dir,"cluster")
    if not os.path.exists(cluster_output_dir):
        os.makedirs(cluster_output_dir)
    sort_output_dir=os.path.join(output_dir,"sort")
    if not os.path.exists(sort_output_dir):
        os.makedirs(sort_output_dir)
    log_output_dir=os.path.join(output_dir,"log")
    if not os.path.exists(log_output_dir):
        os.makedirs(log_output_dir)

    # Multi-threaded execution of cluster and clusterSort
    def process_chro(chro):
        
        # Running cluster
        chro_file_1=os.path.join(output_dir,"chromosome",f"{chro}.num")
        cluster_file=os.path.join(cluster_output_dir,f"{chro}.cluster")
        cluster_log=os.path.join(log_output_dir,f"{chro}.cluster.log")
        console.print(f"Run cluster for {chro}", style="white")
        with open(cluster_log, "w") as f:
            run_success = subprocess.run([cluster_exe, "-i", chro_file_1, "-o", cluster_file], stdout=f, stderr=subprocess.STDOUT).returncode
        if run_success!=0:
            console.print(f"Error: {cmd} failed", style="red")
            return False

        # Run Destruction-Fusion-Cycle.py
        sort_prefix=os.path.join(sort_output_dir,f"{chro}")
        sort_log=os.path.join(log_output_dir,f"{chro}.sort.log")
        
        console.print()
        console.print(f"Run Destruction-Fusion-Cycle for {chro}", style="white")
        with open(sort_log, "w") as f:
            cmd_list=["python",Destruction_Fusion_Cycle_script,
                      "-i",cluster_file,
                      dfc_cmd,
                      "-o",sort_prefix]
            cmd=" ".join(cmd_list)
            run_success = subprocess.run(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT).returncode
        if run_success!=0:
            console.print(f"Error: {cmd} failed", style="red")
            return False
    
    # Set the number of threads
    if thread is None:
        thread=len(chros)
    else:
        thread=min(thread,len(chros))
    console.print(f"Use {thread} threads to sort geneindex", style="white")

    # Multiprocess processing
    with Pool(thread) as p:
        p.map(process_chro, chros)

def reSortGeneindex(chro_list_file:str,output_dir:str,thread: int=None)->None:
    
    resort_script=os.path.join(os.path.dirname(__file__),"reorder_geneindex_by_sorted_num.py")
    if not os.path.exists(resort_script):
        console.print(f"{resort_script} does not exist", style="red")
        raise FileNotFoundError(f"{resort_script} does not exist")
    
    with open(chro_list_file) as f:
        chros = [line.strip() for line in f if line.strip()]
    
    console.print()
    console.print(f"ReSort geneindex by sorted num file...", style="white")
    def process_chro(chro):
        chro_file_1=os.path.join(output_dir,"chromosome",f"{chro}.num")
        chro_file_2=os.path.join(output_dir,"chromosome",f"{chro}.txt")
        sort_file=os.path.join(output_dir,"sort",f"{chro}.sorted.num")
        if not os.path.exists(chro_file_1):
            console.print(f"{chro_file_1} does not exist", style="red")
            raise FileNotFoundError(f"{chro_file_1} does not exist")
        if not os.path.exists(chro_file_2):
            console.print(f"{chro_file_2} does not exist", style="red")
            raise FileNotFoundError(f"{chro_file_2} does not exist")
        if not os.path.exists(sort_file):
            console.print(f"{sort_file} does not exist", style="red")
            raise FileNotFoundError(f"{sort_file} does not exist")
        # Running the `reorder_geneindex_by_sorted_num.py` script will reorder the `geneindex` according to the order in which the `num` files are sorted.
        console.print(f"ReSort geneindex for {chro}", style="white")
        final_file=os.path.join(output_dir,"sort",f"{chro}.sorted.txt")
        cmd=f"python {resort_script} -n {chro_file_1} -s {sort_file} -g {chro_file_2} -o {final_file}"
        run_success = subprocess.run(cmd, shell=True).returncode
        if run_success!=0:
            console.print(f"Error: {cmd} failed", style="red")
            raise RuntimeError(f"Error: {cmd} failed")        
    with Pool(thread) as p:
        p.map(process_chro, chros)

# Integrate information from different chromosomes
def Integrate_information(chro_list_file: str,output_dir: str,sample_list_file: str) -> None:
    console.print()
    console.print(f"Integrate information from different chromosomes...", style="white")

    with open(chro_list_file) as f:
        chros = [line.strip() for line in f if line.strip()]
    
    # Check input file
    randomScaffold_file=os.path.join(output_dir,"chromosome","randomScaffold.txt")
    if not os.path.exists(randomScaffold_file):
        raise FileNotFoundError(f"{randomScaffold_file} does not exist")
    for chro in chros:
        sorted_file=os.path.join(output_dir,"sort",f"{chro}.sorted.txt")
        if not os.path.exists(sorted_file):
            raise FileNotFoundError(f"{sorted_file} does not exist")
    
    b=0
    d=[]
    for chro in chros:
        sorted_file=os.path.join(output_dir,"sort",f"{chro}.sorted.txt")
        with open(sorted_file) as f:
            for line in f:
                SGID=f"SG{b:07d}"  # SG0000001
                d.append(f"{SGID}\t{chro}\t{line}")
                b+=1
    with open(randomScaffold_file) as f:
        for line in f:
            SGID=f"SG{b:07d}"
            d.append(f"{SGID}\trandomScaffold\t{line}")
            b+=1

    # Read sample list file
    sample_list=[]
    with open(sample_list_file) as f:
        for line in f:
            sample=line.strip()
            if sample:
                sample_list.append(sample)
    header=["Orthogroup","chromosome",*sample_list]
    with open(os.path.join(output_dir,"ordoPan.txt"),"w") as f:
        f.write("\t".join(header)+"\n")
        for i in d:
            f.write(i)
    
if __name__ == "__main__":
    import argparse
    helptext = """Main script with sub-script parameters."""
    parser = argparse.ArgumentParser(
        description=helptext,
        formatter_class=argparse.RawTextHelpFormatter
    )
    # === Main script parameters ===
    parser.add_argument("-b", "--bed", type=str, required=True, help="Directory containing BED files for each sample")
    parser.add_argument("-s", "--sample", type=str, required=True, help="File containing the list of sample names")
    parser.add_argument("-c", "--chromosome", type=str, required=True, help="File containing the list of chromosomes")
    parser.add_argument("-g", "--geneindex", type=str, required=True, help="Gene index file")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output directory")
    parser.add_argument("-t", "--thread", type=int, default=None, help="Number of threads to use (default: all available)")
    parser.add_argument("--resume-step", choices=["start", "sort"], default="start", help=(
        "Resume pipeline from a given step:\n"
        "  start : run full pipeline (default)\n"
        "  sort  : resume after gene index sorting\n")
    )
    parser.add_argument("--overwrite", action="store_false", help="Overwrite existing files, default is True")

    # === Subscript Destruction-Fusion-Cycle parameters ===
    default_rate_list = [0.7, 0.8, 0.9, 0.95]
    parser.add_argument("-m", "--minimum-genome-length", type=int, default=100,
                        help="Minimum genome length (avoid inaccurate scores)")
    parser.add_argument("--lns-n", type=int, default=10, help="Parameter -n for clusterSort (default: 10)")
    parser.add_argument("--lns-n-iter", type=int, default=1, help="Parameter -n for clusterSort in destruction loop (default: 1)")
    parser.add_argument("-r", "--rate", type=float, nargs='+', default=default_rate_list,
                        help=f"Fusion thresholds (default: {' '.join(map(str, default_rate_list))})")
    parser.add_argument("--minimum-cluster-length", type=int, default=15,
                        help="Shortest cluster size corrected by finalAdjustment")
    parser.add_argument("-n", "--max-iter", type=int, default=5,
                        help="Max iterations for each threshold (default: 5)")

    args = parser.parse_args()

    # === Main script parameters ===
    bed_dir = args.bed
    sample_file = args.sample
    chro_list = args.chromosome
    geneindex = args.geneindex
    output_dir = args.output
    thread = args.thread
    resume_step = args.resume_step
    overWrite = args.overwrite

    # === Destruction-Fusion-Cycle.py Sub-script parameter collection ===
    dfc_params = {
        "--minimum-genome-length": args.minimum_genome_length,
        "--lns-n": args.lns_n,
        "--lns-n-iter": args.lns_n_iter,
        "--rate": ' '.join(map(str, args.rate)),
        "--minimum-cluster-length": args.minimum_cluster_length,
        "--max-iter": args.max_iter,
    }
    # Automatic splicing command parameters
    dfc_cmd = ' '.join([f"{k} {v}" for k, v in dfc_params.items()])
    
    startTime=time.time()
    if resume_step=="start":
        print("Run full pipeline")
        console.print("Check input files...", style="white")
        check_sample_file_exist(sample_file)
        check_chro_exist(chro_list)
        check_bed_exist(sample_file, bed_dir)
        split_chro(chro_list, sample_file, geneindex, bed_dir,output_dir,overWrite=overWrite)
        sort_geneindex(chro_list,output_dir,thread=thread,overWrite=overWrite)
        reSortGeneindex(chro_list,output_dir,thread=thread)
        Integrate_information(chro_list,output_dir,sample_file)
    elif resume_step=="sort":
        console.print("Resume after gene index sorting", style="white")
        reSortGeneindex(chro_list,output_dir,thread=thread)
        Integrate_information(chro_list,output_dir,sample_file)
    endTime=time.time()
    duration=endTime-startTime
    hours, rem = divmod(duration, 3600)
    minutes, seconds = divmod(rem, 60)
    console.print(f"Total time: {int(hours):02}:{int(minutes):02}:{int(seconds):02}", style="white")
    console.print("All done!", style="bold green")
