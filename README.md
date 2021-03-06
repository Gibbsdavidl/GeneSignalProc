

# GeneSignalProc #

Application of graph signal processing to gene expression data.
Single sample gene set scoring.

0. Build a gene-gene functional graph from a collection of gene sets.
   a.) a sparse matrix is returned representing the adjacency matrix
   b.) a list of genes is returned, which are the row labels for the adjacency matrix
   c.) also saves a python pickle of igraph objects and gene ID dictionary
   d.) can use a gene lists to control the 'gene vocabulary' of the graph

1. Filter the gene expression data.
   a.) uses a matrix of gene expression, samples in rows, genes in cols
   b.) and a graph constructed from gene sets (src/see makeGraphs.py)
   c.) data is filtered to match the gene labels of the graph

2. Score samples for each gene set.
   a.) Scores are internally normalized by comparison to a background of
       subsampled gene sets (from gene set graph).

3. Compare the gene set scores
   a.) based on predictive capability using RandomForests
   b.) rank the sets based on cross validation


=======
To Run:
=======

python3.6 GeneSignalProc/src/main.py --help

<or>

python3.6 GeneSignalProc/src/main.py

        print("makegraphs args: ")
        print("    -m Mode")
        print("    -d data dir")
        print("    -o output dir")
        print("    -g genesets file as .gmt")
        print("    -a adjacency file if available")
        print("    -e gene list file if available")
        print("    -n number of subgraphs")
        print("    -x max size of subgraphs")
        print("    -t threshold for geneset overlaps")
        print("    -c number of cores")

<or>

python3.6 GeneSignalProc/src/main.py

        print("setscoring args: ")
        print("    -m Mode")
        print("    -d data dir")
        print("    -o output dir")
        print("    -r expression data, samples in rows, genes in columns, first col is sample name, first row is gene names")
        print("    -p phenotype file if available, same order as expression data")
        print("    -a adjacency file")
        print("    -s subgraph file")
        print("    -t threshold for geneset overlaps")
        print("    -x threshold for graph edges")
        print("    -e gene list")
        print("    -g genesets file as .gmt")
        print("    -f filter name")
        print("    -l number of scale levels")
        print("    -c num cores")

Examples:

Working with the set of Cancer Hallmark gene sets.

python3.6 GeneSignalProc/src/main.py
-m makegraphs
-d /gene_set_work/
-g h.all.v6.2.symbols.gmt
-c 4
-a h.all.v6.2.symbols.gmt_adjmat.tsv.gz
-e h.all.v6.2.symbols.gmt_genes.tsv.gz
-n 10
-x 7

python3.6 GeneSignalProc/src/main.py
-m setscoring
-d /gene_set_work/
-r ivy20t.tsv
-g h.all.v6.2.symbols.gmt
-a h.all.v6.2.symbols.gmt_adjmat.tsv.gz
-s h.all.v6.2.symbols.gmt_subgraphs.tsv.gz
-e h.all.v6.2.symbols.gmt_genes.tsv.gz
-p ivy20_pheno.tsv
-c 4
-l 10
-t 0.2
-x 0.05


