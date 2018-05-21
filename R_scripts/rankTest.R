scores <- read.table("setscores.txt", sep='\t', header=F)
pheno <- read.table("phenotype.tsv")
analy <- read.table("analyout.tsv", sep='\t', header=T)

cbind(scores[,40:43], Pheno=pheno)

gs <- c(33, 35, 37, 20, 21, 22, 24, 25, 26, 27, 28, 29)

inputs <- read.table("filtered_2.txt", sep='\t')

r1 <- rank(inputs[1,])[gs]
r2 <- rank(inputs[2,])[gs]
r3 <- rank(inputs[3,])[gs]

subgraphs <- read.table("all_subgraphs.txt", stringsAsFactors=F)
commaCounts <- str_count(subgraphs[,1], ',')
idx <- which(commaCounts == 11)
sgs <- subgraphs[idx,]
sgsi <-as.numeric(str_split(sgs[1],',')[[1]])

s1 <- rank(inputs[1,])[sgsi]
s2 <- rank(inputs[2,])[sgsi]
s3 <- rank(inputs[3,])[sgsi]

data.frame(r1,r2,r3,s1,s2,s3)

data.frame(r1+r2+r3,s1+s2+s3)

sum(r1+r2+r3)

sgSums <- c()
for (si in sgs) {
sgsi <-as.numeric(str_split(si,',')[[1]])
s1 <- rank(inputs[1,])[sgsi]
s2 <- rank(inputs[2,])[sgsi]
s3 <- rank(inputs[3,])[sgsi]
sgSums <- c(sgSums, sum(s1+s2+s3))
 }
