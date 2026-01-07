
import sys
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description="Draw unlink bar chart")
parser.add_argument("-i","--input",required=True, type=str, help="Input sorted num file")
parser.add_argument("-o","--output",required=True ,type=str, help="Output PNG file")
args = parser.parse_args()
sorted_num_file = args.input
OUTpng = args.output



# 获取src目录路径
src_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(src_dir))

# 然后导入
import cluster


geneindex=cluster.readGeneIndex(sorted_num_file)
clusters=cluster.geneindex2clusters(geneindex)
speciesNum=clusters[0].speciesNum


b=[]
last_cluster=clusters[0]
for this_cluster in clusters[1:]:
    unlink_num=last_cluster.unlink_num(this_cluster)
    b.append(unlink_num)
    last_cluster=last_cluster+this_cluster
    
import matplotlib.pyplot as plt
#绘制柱状图
plt.bar(range(len(b)), b)
plt.savefig(OUTpng, dpi=300)

