#simsets.R

# need a set of sets, to find an ordering over all.
# related to contigs, but with repeated sets.

genes <- sapply(1:500, function(i) paste(sample(letters, size=5, replace=T),collapse=""))

sets <- list()
idx <- 1
while(sum(genes %in% unlist(sets)) < length(genes)) {
  sampsize <- abs(ceiling(rnorm(mean=30, sd=20, n=1)))
  print(sampsize)
  set1 <- sample(genes, size=sampsize, replace=T)
  sets[[idx]] <- set1
  idx <- idx+1
}
