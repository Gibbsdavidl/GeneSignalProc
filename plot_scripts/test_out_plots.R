library(readr)
library(ggplot2)
library(igraph)
library(ggraph)


setmat <- read_tsv("setmatrix.tsv", col_names=F)
refset <- which(setmat[2,] == 1)-1

scrmat <- as.matrix(read_tsv("scorematrix.tsv", col_names=F))
rownames(scrmat) <- colnames(scrmat)
g <- igraph::graph_from_adjacency_matrix(adjmatrix=scrmat, mode='undirected', weighted=T, diag=F)

res0 <- read_tsv("analyout.tsv")
set1 <- (res0$genes)[20]
set1 <- str_sub(set1, start=2, end=str_length(set1)-1)
set1 <- sapply(str_split(set1, ',')[[1]], as.numeric)
set1 <- as.numeric(set1+1)

nodeCols <- rep("TF",120)
nodeCols[refset] <- "FN"
nodeCols[set1] <- "FP"
nodeCols[intersect(set1, refset)] <- "TP"

ggraph(g) + geom_node_point(aes(col=as.factor(nodeCols))) + geom_edge_fan()

V(g)$nodelabel <- as.character(1:120)
ggraph(g, 'igraph', algorithm = 'nicely') + geom_node_label(aes(label = nodelabel))+ geom_node_point(aes(col=as.factor(nodeCols))) + geom_edge_fan()
