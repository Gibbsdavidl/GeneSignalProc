# Title     : TODO
# Objective : TODO
# Created by: davidgibbs
# Created on: 5/1/18

library(ggplot2)

scores <- read.table("setscores.txt", sep='\t', header=F)
analy <- read.table("analyout.tsv", sep='\t', header=T)
pheno <- read.table("phenotype.tsv")


# last sets in analy will be the most predictive

qplot(y=scores[,42], x=as.factor(pheno$V1), geom="boxplot") + geom_point()


> cbind(scores[,40:43], Pheno=pheno)
     V40   V41   V42   V43 V1
1  0.985 1.000 1.000 1.000  1 *
2  0.985 1.000 1.000 1.000  1 *
3  1.000 1.000 1.000 1.000  1 *
4  0.945 0.730 0.655 0.655  0
5  1.000 1.000 1.000 1.000  1 *
6  0.030 0.305 0.530 0.530  0
7  0.735 0.705 0.970 0.970  0
8  1.000 1.000 1.000 1.000  1 *
9  0.975 0.630 0.265 0.265  0
10 0.760 0.275 0.420 0.420  0
11 0.985 1.000 1.000 1.000  1 *
12 1.000 1.000 1.000 1.000  1 *
13 0.370 0.075 0.175 0.175  0
14 0.475 0.015 0.090 0.090  0
15 0.735 0.495 0.290 0.290  0
16 0.985 1.000 1.000 1.000  1 *
17 0.985 1.000 1.000 1.000  1 *
18 0.700 0.465 0.325 0.325  0
19 0.955 0.985 1.000 1.000  1 *
20 0.645 0.250 0.380 0.380  0

