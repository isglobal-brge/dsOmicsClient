rm(list=ls())
load("c:/Juan/CREAL/BayesianPrediction/Bayesian_clock/paper/data/GSE66351.Rdata")

set.seed(1234)
ss <- sampleNames(gse)
i <- sample(ss, 100, replace=FALSE)
gse1 <- gse[,i]
gse2 <- gse[, ss[!ss%in%i]]
dim(gse1)
dim(gse2)
