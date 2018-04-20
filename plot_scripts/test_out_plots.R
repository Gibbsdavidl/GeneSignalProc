library(readr)
library(ggplot2)
library(igraph)

setmat <- read_tsv("setmatrix.tsv", col_names=F)
refset <- which(setmat[2,] == 1)-2

scrmat <- as.matrix(read_tsv("scorematrix.tsv", col_names=F))

g <- igraph::graph_from_adjacency_matrix(adjmatrix=scrmat, mode='undirected', weighted=T, diag=F)

res0 <- read_tsv("analyout.tsv")
set1 <- (res0$genes)[length(res0$genes)]
set1 <- str_sub(set1, start=2, end=str_length(set1)-1)
set1 <- sapply(str_split(set1, ',')[[1]], as.numeric)


nodeCols <- rep(1,120)
nodeCols[refset] <- 2
nodeCols[set1] <- 3
nodeCols[intersect(set1, refset)] <- 4

ggraph(g) + geom_node_point(aes(col=as.factor(nodeCols))) + geom_edge_fan()

