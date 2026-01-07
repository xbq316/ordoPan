cd /mnt/d/project/泛基因组排序/ordoPan-2025.12.29/data

species_num=60
gene_num=5000
seq 1 100 | xargs -P 10 -I {} python ../ordopan.py imitate -sn ${species_num} -gn ${gene_num} -s {} -o ${species_num}.${gene_num}.{}

species_num=60
gene_num=1500
seq 1 100 | xargs -P 10 -I SEED python ../ordopan.py imitate -sn ${species_num} -gn ${gene_num} -s SEED -o ${species_num}.${gene_num}.SEED --Orthogroup
