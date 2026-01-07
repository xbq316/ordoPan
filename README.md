- 编译程序&安装依赖
```bash
# g++ 应当满足 C++ 17 标准，下载 compilers 以避免 服务器 gblic 版本过低问题
conda install -c conda-forge compilers

# 本流程依赖于 diamond 进行序列比对
conda install diamond -c bioconda

# python版本应当 > 3.8
pip install -r requirements.txt
make
```

- 手册和版本
```bash
# 主命令
python ordopan.py -h
python ordopan.py -v

# 子命令手册
python ordopan.py ordo -h

```

---

# 一.准备输入文件

- 输入数据： file.pep file.bed sample.list chro.list

- **file.pep** 格式如下：
```text
>gene1
MATGENRTVQENLKKHLAVSVRNIQWSYGIFWSISASQPGV
SLLAKSASLLTVVCFPFLGGVLEIGTTEH
>gene2
MEKRRSPRSSKSPGTAT
MDSLRFLRILLVSASILLVSLLYTVTAGDIVHQDDLAPKKPGCENDFVLV
```


---

- **file.bed** 格式如下(支持be)：
```text
A01     831     1437    BnaA01g00010D   .   +
A01     1487    2436    BnaA01g00020D   .   +
A01     2665    5455    BnaA01g00030D   .   +
A01     8421    9623    BnaA01g00040D   .   +
```
---

- **chro.list** 格式如下(标注常染色体信息)：
```
A01
A02
A03
A04
A05
A06
```

```bash
cut -f1 example/Bna_Darmor_v10.0.bed | sort | uniq > example/chro.list
# 后续手动筛选常染色体
```

- sample.list 格式如下：
```txt
file1
file2
file3
file4
...

```

---

# 二.构建泛基因组

## 1.获取泛基因集
- 获取 geneindex 文件
```
python ordopan.py synpan -i input_dir -o output_dir -t thread

```

## 2.对泛基因组排序
```bash
python ordopan.py ordo -b bed_dir -s sample.list -c chro.list -g geneindex.final.txt -o output_dir -t thread_num

```

# 三.子脚本工具
## 1.对单个染色体进行排序
> 若是某个染色体的排序不理想，可使用以下脚本重新对该染色体进行排序
```bash
# 1.基本运行
python ordopan.py sort -i A01.cluster -o outdir/A01

# 2.在排序后的基础上进一步迭代排序
python ordopan.py sort -i A01.cluster -s A01.sorted.num -o outdir/A01

# 3.迭代次数设置
python ordopan.py sort -i A01.cluster  --lns-n 30 -o outdir/A01
# --lns-n 单位为万次

# 4.对lns结果的边进行切割，指定阈值列表
python ordopan.py sort -i A01.cluster -r 0.95 0.98 1 -o outdir/A01
```

## 2.对指定染色体的 geneindex 以上述排序为参考重新调整顺序
```bash
# 对所有染色体都进行调整，根据已经排好序的num文件(file.sorted.num)
python ordopan.py ordo -b bed_dir -s sample.list -c chro.list -g geneindex.final.txt -o output_dir -t thread_num --resume-step sort
# --resume-step sort 表示从排序后的文件开始运行

# 对指定染色体进行调整
python ordopan.py resort -n A01.num -s A01.sorted.num -g A01.txt -o A01.sorted.txt

```

## 3.生成模拟数据
```bash
python ordopan.py imitate -gn 5000 -c 1 2 3 -gn 5000 -sn 50  -o imitate
# -gn 5000 : 每条染色体5000个基因左右
# -c : 染色体列表
# -sn : 基因组数量
# -o : 输出目录
```
