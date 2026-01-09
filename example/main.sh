
# 若是没有泛基因组文件，则可以使用以下脚本构建
# python ../ordopan.py synpan -i data/maize -t 10 -o data/maize.synpan

python ../ordopan.py ordo -b data/maize -s data/maize/sample.list -c data/maize/chro.list -g data/maize.synpan/geneindex.final.txt -t 10 -o data/maize.ordo

