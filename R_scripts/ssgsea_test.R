# Title     : TODO
# Objective : TODO
# Created by: davidgibbs
# Created on: 10/1/18

library(GSVA)

S = matrix(abs(rnorm(10000)), ncol=1)
geneNames <- sapply(1:10000, function(i) paste(sample(letters, 10), collapse=''))
rownames(S) = geneNames
save(S, file='ssgsea_test_dat.rda')

T <- S
T <- as.data.frame(T)
T <- cbind(data.frame(GeneID=rownames(T)), T)
write.table(T, file='ssgsea_test_dat_py.tsv', sep='\t', row.names=F, col.names=F, quote=F)
write.table(G, file='ssgsea_test_dat_gs.tsv', row.names=F, col.names=F, quote=F)

Tt <- t(T[,2])
colnames(Tt) <- rownames(T)
Tt <- as.data.frame(Tt)
Tt <- cbind(data.frame(SameID='a'), Tt)
write.table(Tt, file='ssgsea_test_dat_py.tsv', sep='\t', row.names=F, col.names=T, quote=F)
write.table(G, file='ssgsea_test_dat_gs.tsv', row.names=F, col.names=F, quote=F)


load('ssgsea_test_dat.rda')
G <- head(rownames(s), n=100)

LG <- list(G)

N <- nrow(S)  # number of genes
P <- ncol(S)  # number of samples

R <- apply(S, 2, function(x,N) as.integer(rank(x)), N)  # not sure what this N does
O <- order(R[, 1], decreasing=TRUE)

rownames(R) <- rownames(S)
L = R[O,]   #################### gene expression ranks ordered ...

# head(L)
#apmwizydch iwctyqkngu jtnkvfcydx jpsnkrdvfi pneldjmrsf ehfjbkagxv
#     10000       9999       9998       9997       9996       9995

#apmwizydch
#  4.082186


gsva(expr=S, gset.idx.list=LG, method='ssgsea', tau=2, ssgsea.norm=F)
#2501.529


#[1,] 1606.651  # for tau = 1

PG <- function(G,L,i,w) {
    # i goes from 1 .. N in calling function
    # L is the expression rank, in order, max N, min 1
    # G is the gene set index
    # ci is the column of O .. will be 1

    val <- 0
    js <- 1:(i-1)
    jsInG <- js[js %in% G]
    norm <- sum(L[G]^w)

    if (length(jsInG) > 0) {

        for (j in jsInG) {
            val <- val + (L[j]^w / norm)
        }

        return(as.numeric(val))

    } else {
        return(0)
    }
}


PNG <- function(G,L,i) {

    # i goes from 1 .. N in calling function
    # L is the expression rank, in order, max N, min 1
    # G is the gene set index
    # ci is the column of O .. will be 1

    denom <- 1/(length(L) - length(G))   #
    js <- 1:(i-1)
    jsNotInG <- js[!js %in% G]
    res0 <- length(jsNotInG) * denom
    return(res0)
}


ES <- function(G, L, w) {
    sum(
        sapply(1:length(L), function(i){
            PG(G,L,i,w) - PNG(G,L,i)
        })
    )
}

G <- match(table=names(L), x=rownames(S)[1:100])

ES(G, L, 2)
[1] 2501.528


# python score:
#['rizskgubye', 'abevgukqxd', 'tfzmqlnrxh', 'igtnlxmupr', 'dribhqkwsp']

[2501.5286237144633]