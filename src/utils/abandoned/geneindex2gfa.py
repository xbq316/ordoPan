import networkx as nx
import argparse
import os
import sys, os
sys.path.append(os.path.dirname(__file__))
import cluster


class Edge:
    def __init__(self, name1: str,strand1: str,name2: str,strand2: str,overlap:int=0):
        if strand1 not in ["+","-"]:
            raise ValueError(f"Invalid strand1: {strand1}")
        if strand2 not in ["+","-"]:
            raise ValueError(f"Invalid strand2: {strand2}")
        if overlap<0:
            raise ValueError(f"Invalid overlap: {overlap}")
        self.name1 = name1
        self.strand1 = strand1
        self.name2 = name2
        self.strand2 = strand2
        self.overlap = overlap

    def __str__(self):
        return f"L\t{self.name1}\t{self.strand1}\t{self.name2}\t{self.strand2}\t{self.overlap}M"

    def coreInfo(self):
        return tuple([self.name1,self.strand1,self.name2,self.strand2])

class Edges:
    def __init__(self, edges: list[Edge]=[]):
        self.edges = edges
        self.edge_set=set()
        for edge in edges:
            edge_core_info=edge.coreInfo()
            if edge_core_info in self.edge_set:
                raise ValueError(f"Duplicate edge: {edge_core_info}")
            self.edge_set.add(edge_core_info)

    def add(self,edge:Edge):
        edge_core_info=edge.coreInfo()
        if edge_core_info in self.edge_set:
            pass
        else:
            self.edges.append(edge)
            self.edge_set.add(edge_core_info)

    def add_edges(self,edges: list[Edge]):
        for edge in edges:
            edge_core_info=edge.coreInfo()
            if edge_core_info in self.edge_set:
                raise ValueError(f"Duplicate edge: {edge_core_info}")
            self.edge_set.add(edge_core_info)
            self.edges.append(edge)
    def __iter__(self):
        return iter(self.edges)

def get_bestchromosome(geneindex: list[list[str]],bed_dir: str,sample_list: list[str])-> list[str]:
    geneindex=list(zip(*geneindex))
    chromosome_list=[]
    for species_index,sample_name in enumerate(sample_list):
        
        bed_list=cluster.read_bed6(os.path.join(bed_dir,f"{sample_name}.bed"))
        genename2chromosome={}
        for bed in bed_list:
            genename2chromosome[bed.name]=bed.chromosome

        for genes in geneindex[species_index]:
            if genes=="-":
                continue
            elif "," in genes:
                for gene in genes.split(","):
                    try:
                        chromosome=genename2chromosome[gene]
                        chromosome_list.append(chromosome)
                    except KeyError:
                        raise KeyError(f"gene {gene} in {sample_name} not found in bed file")
            else:
                try:
                    chromosome=genename2chromosome[genes]
                    chromosome_list.append(chromosome)
                except KeyError:
                    raise KeyError(f"gene {genes} in {sample_name} not found in bed file")

    bestchromosome = max(chromosome_list, key=chromosome_list.count)
    return bestchromosome


def get_strand_maxCount(geneindex: list[list[str]],bed_dir: str,sample_list: list[str],bestchromosome: str)-> list[str]:
    geneindex=list(zip(*geneindex))

    result=[]
    all_species_strand=[]
    for species_index,sample_name in enumerate(sample_list):
        
        bed_list=cluster.read_bed6(os.path.join(bed_dir,f"{sample_name}.bed"))
        genename2gene={}
        for gene in bed_list:
            genename2gene[gene.name]=gene

        this_species_strand=[]

        for genes in geneindex[species_index]:
            if genes=="-":
                this_species_strand.append(None)
                continue
            elif "," in genes:
                a=[]
                for gene in genes.split(","):
                    try:
                        this_species_chromosome=genename2gene[gene].chromosome
                        if this_species_chromosome!=bestchromosome:
                            continue
                        else:
                            strand=genename2gene[gene].strand
                            a.append(strand)
                    except KeyError:
                        raise KeyError(f"gene {gene} in {sample_name} not found in bed file")
                if a:
                    this_species_beststrand=max(set(a),key=a.count)
                    this_species_strand.append(this_species_beststrand)
                else:
                    this_species_strand.append(None)
            else:
                try:
                    this_species_chromosome=genename2gene[genes].chromosome
                    if this_species_chromosome!=bestchromosome:
                        this_species_strand.append(None)
                        continue
                    strand=genename2gene[genes].strand
                    this_species_strand.append(strand)
                except KeyError:
                    raise KeyError(f"gene {genes} in {sample_name} not found in bed file")
        all_species_strand.append(this_species_strand)
    all_species_strand=list(zip(*all_species_strand))

    for strand in all_species_strand:
        a=[]
        for s in strand:
            if s:
                a.append(s)
        if a:
            result.append(max(set(a),key=a.count))
        else:
            result.append(None)
    if None in result:
        return None
    else:
        return result

def get_edges_paths(geneindex: list[list[str]],bed_dir: str,sample_list: list[str],bestchromosome: str,)->tuple[Edges,dict[str,list[str]],dict[str,list[str]]]:
    geneindex=list(zip(*geneindex))
    
    edges=Edges([])
    paths={}
    genome={}

    for species_index,sample_name in enumerate(sample_list):
        
        bed_list=cluster.read_bed6(os.path.join(bed_dir,f"{sample_name}.bed"))
        genename2gene={}
        for gene in bed_list:
            genename2gene[gene.name]=gene

        this_species_genes=[]
        this_species_path=[]
        this_species_genome=[]

        for this_index,genes in enumerate(geneindex[species_index]):
            this_index+=1
            if genes=="-":
                continue
            elif "," in genes:
                for gene in genes.split(","):
                    try:
                        this_species_chromosome=genename2gene[gene].chromosome
                        if this_species_chromosome!=bestchromosome:
                            continue
                        else:
                            this_species_genes.append((gene,this_index))
                    except KeyError:
                        raise KeyError(f"gene {gene} in {sample_name} not found in bed file")
            else:
                try:
                    this_species_chromosome=genename2gene[genes].chromosome
                    if this_species_chromosome!=bestchromosome:
                        continue
                    else:
                        this_species_genes.append((genes,this_index))
                except KeyError:
                    raise KeyError(f"gene {genes} in {sample_name} not found in bed file")
        this_species_genes.sort(key=lambda x:genename2gene[x[0]].start)

        for i in range(len(this_species_genes)-1):
            genename1,genename2=this_species_genes[i][0],this_species_genes[i+1][0]
            nodeid1,nodeid2=this_species_genes[i][1],this_species_genes[i+1][1]
            strand1,strand2=genename2gene[genename1].strand,genename2gene[genename2].strand
            edges.add(Edge(nodeid1,strand1,nodeid2,strand2))
            this_species_path.append(f"{nodeid1}{strand1}")
            this_species_genome.append(genename1)
        this_species_path.append(f"{nodeid2}{strand2}")
        this_species_genome.append(genename2)
        paths[sample_name]=this_species_path
        genome[sample_name]=this_species_genome

    return edges,paths,genome



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", type=str, required=True, help="Input geneindex file")
    parser.add_argument("-s","--sample", type=str, required=True, help="Input sample file")
    parser.add_argument("-b","--bed", type=str, required=True, help="Input bed file directory")
    parser.add_argument("-o","--output", type=str, required=True, help="Output file for GFA")
    args = parser.parse_args()

    infile= args.input
    outfile= args.output
    
    beddir=args.bed
    sample_list=cluster.read_sample(args.sample)
    geneindex=cluster.read_geneindex(infile)
    bestchromosome=get_bestchromosome(geneindex,beddir,sample_list)
    ancestor_strand=get_strand_maxCount(geneindex,beddir,sample_list,bestchromosome)
    if ancestor_strand is None:
        raise ValueError("ancestor_strand is None, please check geneindex file")
    
    node_num=len(geneindex)
    node_id_list=list(range(1,node_num+1))
    node_name_list=[f"SG{i:07d}" for i in range(node_num)]

    edges,paths,genome=get_edges_paths(geneindex,beddir,sample_list,bestchromosome)

    with open(outfile,"w") as out:
        out.write(f"H\tVN:Z:1.1\n")
        # 泛基因组节点
        for this_node_id,this_node_name in zip(node_id_list,node_name_list):
            out.write(f"S\t{this_node_id}\t{this_node_name}\n")
        # 泛基因组边
        for i in range(node_num-1):
            this_node_id=node_id_list[i]
            this_strand=ancestor_strand[i]
            next_node_id=node_id_list[i+1]
            next_strand=ancestor_strand[i+1]
            # out.write(f"L\t{this_node_id}{this_strand}\t{next_node_id}{next_strand}\t*\n")
            edges.add(Edge(this_node_id,this_strand,next_node_id,next_strand))
        # 各基因组边和泛基因组
        for edge in edges:
            out.write(f"{edge}\n")
        # 各个基因组 path 数据
        for sample_name,this_path in paths.items():
            thispath=",".join(this_path)
            out.write(f"P\t{sample_name}\t{thispath}\t*\n")
        for sample_name,this_genome in genome.items():
            thisgenome=",".join(this_genome)
            out.write(f"#G\t{sample_name}\t{thisgenome}\t*\n")
        # 泛基因组边
        out.write("\t".join(["P","ANCESTOR",",".join([f"{index+1}{this_strand}" for index,this_strand in enumerate(ancestor_strand)]),"*\n"]))


if __name__=="__main__":
    
    main()
