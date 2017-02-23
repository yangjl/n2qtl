### Jinliang Yang
### Feb 23th, 2017

install.packages("qtl")
library(qtl)

geno <- read.table("data/Genotypes_filtered_B73_final.txt", header=TRUE)


d <- read.cross("csvs", ".", "listeria_gen.csv", "listeria_phe.csv")

