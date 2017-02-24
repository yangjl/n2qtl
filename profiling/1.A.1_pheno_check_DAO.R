### Jinliang Yang
### Feb 23th, 2017
### data checking


## Install some required packages if have not installed

#install.packages("lme4")
#install.packages(devtools)
#devtools::install_github("jyanglab/g3tools")

library("lme4")
library("g3tools")




pheno <- read.csv("data/1_QTL_AllFamilies_NoOutlier.csv")
dim(pheno)
#687  33

n2 <- subset(pheno, !is.na(X15NT1))

hist(n2$X15NT1, breaks=30)
hist(n2$X15NT1R, breaks=30)


fit <- get_BLUP(data = pheno, model = X15NT1 ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
         + (1 | Latitude), which.factor = "Genotype",
         outfile = "data/blup_15NT1.csv")
get_H2(fit, numerator="Genotype", denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))







fit2 <- get_BLUP(data = pheno, model = X15NT2 ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                + (1 | Latitude), which.factor = "Genotype",
                outfile = "data/blup_15NT2.csv")

get_H2(fit2, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))


fit3 <- get_BLUP(data = pheno, model = X15NT3 ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "Genotype",
                 outfile = "data/blup_15NT3.csv")

get_H2(fit3, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))




b <- read.csv("data/blup.csv")


fit1 <- get_BLUP(data = pheno, model = DT ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                + (1 | Latitude), which.factor = "Genotype",
                outfile = "data/blup_DT.csv")

get_H2(fit1, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit2 <- get_BLUP(data = pheno, model = DPS ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "Genotype",
                 outfile = "data/blup_DPS.csv")

get_H2(fit2, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit3 <- get_BLUP(data = pheno, model = AR ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "Genotype",
                 outfile = "data/blup_AR.csv")

get_H2(fit3, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit4 <- get_BLUP(data = myp, model = SC ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "Genotype",
                 outfile = "data/blup_SC.csv")

get_H2(fit4, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit5 <- get_BLUP(data = myp, model = PDM ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "Genotype",
                 outfile = "data/blup_PDM.csv")

get_H2(fit5, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit6 <- get_BLUP(data = myp, model = PTN ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "Genotype",
                 outfile = "data/blup_PTN.csv")

get_H2(fit6, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit7 <- get_BLUP(data = myp, model = PTN_percent ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "Genotype",
                 outfile = "data/blup_PTN_percent.csv")

get_H2(fit7, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit8 <- get_BLUP(data = myp, model = GDM ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "Genotype",
                 outfile = "data/blup_GDM.csv")

get_H2(fit8, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit9 <- get_BLUP(data = myp, model = GTN ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "Genotype",
                 outfile = "data/blup_GTN.csv")

get_H2(fit9, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit10 <- get_BLUP(data = myp, model = GTN_percent ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "Genotype",
                 outfile = "data/blup_GTN_percent.csv")

get_H2(fit10, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))
