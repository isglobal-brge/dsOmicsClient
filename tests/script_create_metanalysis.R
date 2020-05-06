library(metap)

load("c:/juan/Creal/GitHub/brgedata/data/gse66351_1.rda")
load("c:/juan/Creal/GitHub/brgedata/data/gse66351_2.rda")





ff <- function(x, y){
  
  inner_join(x%>%select("id", "P.Value"), 
             y%>%select("id", "P.Value"), 
             by="id")
}


mm <- Reduce(ff, ans.limma)
colnames(mm)[2:3] <- paste("p", names(ans.limma), sep=".")
p.meta <- unlist(apply(mm[1:10,-1], 1, function(x) allmetap(x, method="sumlog")$p))
cbind(mm[1:10,], p.meta=p.meta)
