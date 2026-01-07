

import sys

infile=sys.argv[1]
window_size=100

a=[]
for i in open(infile):
    i = i.strip().split("\t")
    a.append(i)

species_num=len(a[0])
gene_num=len(a)

b=[]
for i in range(0, gene_num, window_size):
    window = a[i:i+window_size]
    x=0
    for j in window:
        if j.count("-")==species_num-1:
            x+=1
    b.append(x)

import matplotlib.pyplot as plt
plt.plot(b)
plt.savefig(infile+".private.count.png", dpi=300)
