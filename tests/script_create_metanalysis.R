library(metaMA)

load("c:/juan/Creal/GitHub/brgedata/data/gse66351_1.rda")
load("c:/juan/Creal/GitHub/brgedata/data/gse66351_2.rda")
esets <- list(gse1, gse2)
classes <- list(gse1$diagnosis, gse2$diagnosis)

theScores <- pvalcombination(esets,classes)

ff <- function(x, ids){
  ans <- x[ ,"P.Value", drop=FALSE]
  rownames(ans) <- x$id
  out <- as.matrix(ans[ids,])
  out
}
idx <- Reduce(intersect, lapply(ans.limma, function(x) x$id))
pvalList <- lapply(ans.limma, ff)
out <- directpvalcombi(pvalList, c(100, 90))

group <- as.numeric(as.factor(gse1$diagnosis))-1
d.gse1 <- getdF(gse1, group)
d.adj.gse1 <- dstar(d.gse1, length(group))
var.d.adj.gse1 <- sigmad(d.adj.gse1, sum(group==0), sum(group==1))
head(d.adj.gse1)
