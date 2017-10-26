


# making simulated data

g <- sample_growing(500)

gcomp <- components(g)

induced_subgraph(g, v=which(gcomp$membership == 1))

gadj <- get.adjacency(g, type="both", sparse=F)

gadj[gadj == 1] <- rnorm(n=489)

gadj[gadj == 1] <- abs(rnorm(n=489))

gadj[upper.tri(gadj)] <- t(gadj)[upper.tri(gadj)]

write.table(gadj, file="adj_matrix.tsv", sep="\t", row.names=F, col.names=F, quote=F)

genenames <- sapply(1:500, function(i) paste(sample(letters, 5), collapse=""))

exprvals <- abs(rnorm(500))
write.table(data.frame(Genes=genenames, Expr=exprvals), file="sample1.tsv", sep="\t", row.names=F, quote=F)

exprvals <- abs(rnorm(500))
write.table(data.frame(Genes=genenames, Expr=exprvals), file="sample2.tsv", sep="\t", row.names=F, quote=F)

exprvals <- abs(rnorm(500))
write.table(data.frame(Genes=genenames, Expr=exprvals), file="sample3.tsv", sep="\t", row.names=F, quote=F)

exprvals <- abs(rnorm(500))
write.table(data.frame(Genes=genenames, Expr=exprvals), file="sample4.tsv", sep="\t", row.names=F, quote=F)

exprvals <- abs(rnorm(500))
write.table(data.frame(Genes=genenames, Expr=exprvals), file="sample5.tsv", sep="\t", row.names=F, quote=F)

filelist <- c("sample1.tsv", "sample2.tsv", "sample3.tsv", "sample4.tsv", "sample5.tsv")
write.table(filelist, file="filelist.tsv", sep="\t", row.names=F, col.names=F, quote=F)
