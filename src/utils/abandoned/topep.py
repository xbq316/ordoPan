


import sys
import os
from Bio import SeqIO
import warnings


genelistFile=sys.argv[1]
outfile=sys.argv[2]
pepdir="/public/home/bqxiao/data5/06.Bna-anc/14.Bna.orthofinder"

geneList=[]
for i in open(genelistFile):
    i=i.strip()
    geneList.append(i)

genename2pep={i:"" for i in geneList}

pepFileList=[]
for i in os.listdir(pepdir):
    if i.endswith(".pep"):
        pepFileList.append(i)
geneNum=len(pepFileList)


pep_abs_path=[]
for i in pepFileList:
    j=os.path.join(pepdir,i)
    pep_abs_path.append(j)

a=0
for fasta_file in pep_abs_path:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in genename2pep:
            genename2pep[record.id]=str(record.seq)
            a+=1
    if a==geneNum:
        break

with open(outfile,"w") as f:
    for i in geneList:
        if genename2pep[i]=="":
            warnings.warn(f"gene {i} not found in pep dir", category=UserWarning)
        else:
            f.write(f">{i}\n{genename2pep[i]}\n")
