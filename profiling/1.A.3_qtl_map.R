### Jinliang Yang
### Feb 23th, 2017

#install.packages("qtl")
library(qtl)

geno <- read.table("data/Genotypes_filtered_B73_final.txt", header=TRUE)
h <- geno[, 1:3]
names(h) <- c("id", "chr", "pos")
h$chr <- gsub("S|_.*", "", h$id)
h$pos <- gsub(".*_", "", h$id)
h$pos <- as.numeric(as.character(h$pos))/1000000


### Genotypic data
g <- cbind(h, apply(geno[, -1], 2, as.character))

g <- data.frame(lapply(g, function(x) {gsub("U", "-", x)}))

tg <- t(g)
row.names(tg)[2:3] <- c("", "")
write.table(tg, "cache/geno_B73.csv", sep=",", col.names=FALSE)


### phenotypic data for 15NT1
p <- read.csv("data/blup_15NT1.csv")
names(p)[1:2] <- c("id", "a15NT1")
p$id <- paste0("X", p$id)

write.table(p, "cache/pheno_a15NT1.csv", sep=",", row.names=FALSE, quote=FALSE)



#### Read in the data
library(qtl)
d <- read.cross("csvs", dir="cache", "geno_B73.csv", "pheno_a15NT1.csv")


summary(d)
plotPheno(d, pheno.col=2)


############################################################
# Single-QTL analysis
############################################################
d <- calc.genoprob(d, step=1)

out.em <- scanone(d, pheno.col=2)

summary(out.em)

summary(out.em, threshold=3)

plot(out.em)


####
out.hk <- scanone(d, pheno.col=2, method="hk")

plot(out.em, out.hk, col=c("blue", "red"))

plot(out.em, col="blue")
plot(out.hk, col="red", add=TRUE)

plot(out.hk - out.em, ylim=c(-0.3, 0.3), ylab="LOD(HK)-LOD(EM)")

sug <- sim.geno(sug, step=1, n.draws=64)
out.imp <- scanone(sug, method="imp")

plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"))

plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"), chr=c(7,15))

plot(out.imp - out.em, out.hk - out.em, col=c("green", "red"), ylim=c(-1,1))
