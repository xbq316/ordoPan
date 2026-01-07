


import sys
infile=sys.argv[1]

start_cluster=""
end_cluster=""
cluster_time=""
score=""
sort_time=""

sort_time_h=0
sort_time_m=0
sort_time_s=0

for i in open(infile):
    i=i.strip()
    if i.startswith("initial clustering... "):
        start_cluster=i.split(" ")[2]
    if i.startswith("Clustering completed. "):
        end_cluster=i.split(" ")[2]
    if i.startswith("Runtime: "):
        line=i.split(" ")
        h=line[1][:-1]
        m=line[2][:-1]
        s=line[3][:-1]
        cluster_time=f"{h}:{m}:{s}"
    if i.startswith("Total time taken: "):
        line=i.split(" ")
        h=line[3][:-1]
        m=line[4][:-1]
        s=line[5][:-1]
        cluster_time=f"{h}:{m}:{s}"
    if i.startswith("best rate:"):
        score=i.split(" ")[-1]
    if i.startswith("Run time "):
        sort_time=i.split(" ")[2]
    if i.startswith("耗时: "):
        if "小时" in i:
            sort_time_h+=float(i.split(" 小时")[0].split(" ")[-1])
        if "分" in i:
            sort_time_m+=float(i.split(" 分")[0].split(" ")[-1])
        if "秒" in i:
            sort_time_s+=float(i.split(" 秒")[0].split(" ")[-1])
    
if start_cluster and end_cluster and cluster_time:
    print(f"{start_cluster}\t{end_cluster}\t{cluster_time}")

if score and sort_time:
    print(f"{sort_time}\t{score}")
elif score and (sort_time_h or sort_time_m or sort_time_s):
    total_sort_time_s=sort_time_h*3600+sort_time_m*60+sort_time_s
    h=int(total_sort_time_s//3600)
    m=int((total_sort_time_s%3600)//60)
    s=total_sort_time_s%60
    print(f"{h}:{m}:{s}\t{score}")