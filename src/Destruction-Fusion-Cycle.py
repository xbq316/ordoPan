#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: Destruction-Fusion-Cycle.py
# Created time: 2025/08/22 16:45:44
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python Destruction-Fusion-Cycle.py -h



import subprocess
import argparse
import os
import sys
import cluster
import time

# Global setting: stdout line buffering (flush on newline)
sys.stdout.reconfigure(line_buffering=True)

def format_duration(seconds: float) -> str:
    """Format elapsed time: seconds / minutes:seconds / hours:minutes:seconds"""
    if seconds < 60:
        return f"{seconds:.2f} seconds"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        sec = seconds % 60
        return f"{minutes} min {sec:.2f} sec"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        sec = seconds % 60
        return f"{hours} hours {minutes} min {sec:.2f} sec"

def run_command(cmd):
    """Run a command, check if it succeeds, and print elapsed time"""
    print("Running command:", " ".join(map(str, cmd)))
    start = time.time()
    
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    end = time.time()
    duration = end - start
    print("Elapsed time:", format_duration(duration))
    
    if result.returncode != 0:
        print("Command failed:", result.stderr, file=sys.stderr)
        sys.exit(1)
    
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description="Iteratively run clusterSort and finalAdjustment to sort geneindex.")
    parser.add_argument("-i", "--input", required=True, help="Original cluster file (output of cluster)")
    parser.add_argument("-s", "--sort", type=str, help="Existing sort file (output of clusterSort), skip first clusterSort if provided")
    parser.add_argument("--minimum-cluster-length",type=int,default=15,help="The shortest cluster size corrected by the finalAdjustment procedure.")
    parser.add_argument("-n", "--max-iter", type=int, default=5, help="Max iterations for each threshold (default: 5)")
    parser.add_argument("-r", "--rate", type=float,nargs='+', default=[0.7,0.8,0.9,0.95], 
                        help="Fusion thresholds (default: [0.7,0.8,0.9,0.95]); higher means stricter, thus fewer merges")
    parser.add_argument("-o", "--output", required=True,type=str, help="Output file prefix (two files: .sorted.num and .new.cluster)")
    parser.add_argument("-m", "--minimum-genome-length", default=100,help="Minimum genome length (used for scoring to avoid NA or inaccurate scores when some samples have too few genes)")
    parser.add_argument("--lns-n", type=int, default=10, help="Parameter -n for clusterSort (default: 10)")
    parser.add_argument("--lns-n-iter", type=int, default=1, help="Parameter -n for clusterSort in destruction loop (default: 1)")
    args = parser.parse_args()

    original_cluster = args.input
    outFile_prefix = args.output
    max_iter = args.max_iter
    rate_list = args.rate
    minimumGenomeLength= args.minimum_genome_length
    lns_n = args.lns_n
    lns_n_iter = args.lns_n_iter
    minimumClusterLength=args.minimum_cluster_length

    startTime= time.time()

    outfile = f"{outFile_prefix}.sorted.num"
    new_cluster_file = f"{outFile_prefix}.new.cluster"   # cluster file path

    script_dir = os.path.dirname(os.path.abspath(__file__))
    clusterSort_exe= os.path.join(script_dir, "clusterSort")
    finalAdjustment_exe= os.path.join(script_dir, "finalAdjustment")

    # Initial run
    rate = str(rate_list[0])
    if args.sort is None:
        run_command([clusterSort_exe, "-i", original_cluster, "-n",str(lns_n), "-o", outfile])
        
        # run finalAdjustment: --minimum-cluster-length -> -n 
        run_command([finalAdjustment_exe, "-c", original_cluster,"-n",str(minimumClusterLength),"-g", outfile, "-o", outfile])
        cluster.sort2cluster(originalClusterFile=original_cluster, rate=float(rate), sortedClusterFile=outfile, outClusterFile=new_cluster_file)
        init_sort_geneindex_file=outfile
    else:
        init_sort_geneindex_file=args.sort
        cluster.sort2cluster(originalClusterFile=original_cluster, rate=float(rate), sortedClusterFile=init_sort_geneindex_file, outClusterFile=new_cluster_file)
    
    this_geneindex=cluster.readNum(init_sort_geneindex_file)
    best_score=cluster.geneindex2rankScore(this_geneindex,int(minimumGenomeLength))
    geneindex_data=this_geneindex.copy()

    cluster_data = cluster.read_clusters(new_cluster_file)
    new_clusters_num = len(cluster_data)
    last_clusters_num = new_clusters_num
    best_cluster_data = cluster_data.copy()
    best_rate=rate

    print(f"Initial cluster count: {new_clusters_num}")

    # Iteratively run; stop if cluster count stabilizes
    for rate in rate_list:
        for i in range(1, max_iter + 1):
            print("*"*50)
            print(f"rate:{rate}-Iteration {i}, current cluster count: {new_clusters_num}")
            
            run_command([clusterSort_exe, "-i", new_cluster_file, "-n",str(lns_n_iter), "-o", outfile])
            run_command([finalAdjustment_exe, "-c", original_cluster, "-g", outfile, "-o", outfile])
            cluster.sort2cluster(originalClusterFile=original_cluster, rate=float(rate), sortedClusterFile=outfile, outClusterFile=new_cluster_file)
            
            
            # Read new cluster count
            cluster_data = cluster.read_clusters(new_cluster_file)
            new_clusters_num = len(cluster_data)

            # Stop if cluster count does not change
            if i > 1 and new_clusters_num == last_clusters_num:
                print("Cluster count unchanged, stopping iteration.")
                break
            
            last_clusters_num = new_clusters_num

            this_geneindex=cluster.readNum(outfile)
            this_score = cluster.geneindex2rankScore(this_geneindex,int(minimumGenomeLength))
            print(f"rate:{rate}-Current connection score: {this_score}")
            if this_score > best_score:
                best_score=this_score
                geneindex_data=this_geneindex.copy()
                best_cluster_data=cluster_data.copy()
                best_rate=rate


        # Run clusterSort one last time to get final sorted file
        run_command([clusterSort_exe, "-i", new_cluster_file, "-n",str(lns_n), "-o", outfile])
        run_command([finalAdjustment_exe, "-c", original_cluster, "-g", outfile, "-o", outfile])
        this_geneindex = cluster.readNum(outfile)
        this_score = cluster.geneindex2rankScore(this_geneindex,int(minimumGenomeLength))
        print(f"rate:{rate}-Current connection score: {this_score}")
        if this_score > best_score:
            best_score = this_score
            geneindex_data = this_geneindex.copy()
            best_cluster_data = cluster_data.copy()
            best_rate = rate
    
    print("*"*50)
    print(f"best rate:{best_rate}-kendallTau: {best_score}")
    geneindex_report=cluster.geneindex_permutation_report(geneindex_data)
    print(geneindex_report)
    print("*"*50)

    cluster.writeGeneIndexNum(geneindex_data, outfile)
    cluster.writeClusterList(best_cluster_data,new_cluster_file)

    print(f"Iteration completed. Final result saved to: {outfile}, best score(kendallTau): {best_score}, best clustering result saved to: {new_cluster_file}")
    endTime = time.time()

    duration = endTime - startTime
    hours = int(duration // 3600)
    minutes = int((duration % 3600) // 60)
    seconds = duration % 60
    print(f"Run time {hours}:{minutes}:{seconds} (hour:minute:second)")
    print("Done!")

if __name__ == "__main__":
    main()
