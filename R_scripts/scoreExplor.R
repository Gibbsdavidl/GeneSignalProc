
library(ggplot2)
library(superheat)
library(igraph)
library(ggraph)
library(dplyr)

trees <- read.table("denovo_trees.tsv", header=T)
t0 <- trees %>% filter(SampleID == 0 & TreeID == 0)
s0t1 <- trees %>% filter(SampleID == 0 & TreeID == 1)
s2t1 <- trees %>% filter(SampleID == 2 & TreeID == 1)
superheat(table(s0t1$Level, s0t1$GeneID), pretty.order.rows = F, pretty.order.cols = F)
superheat(table(s2t1$Level, s2t1$GeneID), pretty.order.rows = F, pretty.order.cols = F)

setscrs <- read.table("setscores.tsv")
pheno  <- read.table("phenotype.tsv")
data.frame(pheno, setscrs[,18:23])
data.frame(pheno, setscrs[,22])
table(pheno$V1, (setscrs[,6] > 0.9))

qplot(x=as.factor(pheno$V1), y=setscrs[,21], geom='boxplot')

scrmat <- as.matrix(read.table("scorematrix.tsv"))
rownames(scrmat) <- colnames(scrmat)
gsnet <- igraph::graph_from_adjacency_matrix(adjmatrix=scrmat, mode='undirected', weighted=T, diag=F)

gsmat <- read.table("setmatrix.tsv")
gsidx <- which(gsmat[2,] == 1)

expr <- read.table("exprdat.tsv", header=T)
superheat(expr[,-1], pretty.order.rows = T, pretty.order.cols = F, bottom.label.text.size=1, yr = pheno$V1)

p1 <- which(pheno == 1)
p0 <- which(pheno == 0)

s1 <- read.table("filtered_files/filtered_1.txt")
e1 <- as.numeric(expr[p1[1],-1])
x <- rbind(e1, s1)
superheat(x, pretty.order.rows = F, pretty.order.cols = F, bottom.label.text.size=1)

s0 <- read.table("filtered_files/filtered_3.txt")
e0 <- as.numeric(expr[p0[1],-1])
x <- rbind(e0, s0)
superheat(x, pretty.order.rows = F, pretty.order.cols = F, bottom.label.text.size=1)

V(gsnet)$nodelabel <- as.character(1:80)
ggraph(gsnet, 'igraph', algorithm = 'nicely') + geom_node_label(aes(label = nodelabel, fill=e1)) + geom_edge_fan()

all(e25 == as.numeric(x[1,]))  
#[1] TRUE

y0 <- as.numeric(x[p0[1],]) # and equals e25
ggraph(gsnet, 'igraph', algorithm = 'nicely') + geom_node_label(aes(label = nodelabel, fill=y0)) + geom_edge_fan()

# unfiltered
boxplot(list(as.numeric(x[1,gsidx]), as.numeric(x[1,-gsidx])))
# level 1
boxplot(list(as.numeric(x[2,gsidx]), as.numeric(x[2,-gsidx])))
# level 2
boxplot(list(as.numeric(x[3,gsidx]), as.numeric(x[3,-gsidx])))
# level 3
boxplot(list(as.numeric(x[4,gsidx]), as.numeric(x[4,-gsidx])))


filtList <- list()
for (li in sort(unique(s2t1$Level))) {
  print(li)
  filtList[[(li+1)]] <- s2t1[s2t1$Level == li,5]
}

boxplot(filtList)



filt <- read.table("filtered_files/filtered_49.txt",header=F)
superheat((filt), pretty.order.cols=F, pretty.order.rows=F, bottom.label.text.size=2)

delta <- matrix(data=0, ncol=80, nrow=9)
for (i in 1:9) {
    delta[i,] <- as.numeric(filt[i,] - filt[i+1,])
}
superheat((delta), pretty.order.cols=F, pretty.order.rows=F, bottom.label.text.size=2)

segm <- delta
for (i in 1:9) {
 dmed <- median(as.numeric(delta[i,]))
 dmad <- mad(as.numeric(delta[i,]))
 segm[i,] <- ifelse(abs(as.numeric(delta[i,])) > (dmed+dmad), 1, 0)
}

superheat((segm), pretty.order.cols=F, pretty.order.rows=F, bottom.label.text.size=2)