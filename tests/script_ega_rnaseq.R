library(resourcer)
library(DSLite)
library(dsOmicsClient)
library(dsBaseClient)
library(tidyverse)


csv <- newResource(name = "counts", 
                      url = "file:c:/Juan/CREAL/ATHLETE/DataSHIELD/EGA/RNAseq/RawGeneExpressionProfiles.csv", 
                      format = "csv")
pheno <- newResource(name = "pheno", 
                     url = "file:c:/Juan/CREAL/ATHLETE/DataSHIELD/EGA/RNAseq/rnaseq_rob_feno.csv", 
                     format = "csv")


pkgs <- c("dsBase", "resourcer", "dsOmics", "SummarizedExperiment", "GenomicRanges", "edgeR")
dslite.server <- newDSLiteServer(resources = list(csv=csv, pheno=pheno),
                                 config = DSLite::defaultDSConfiguration(include=pkgs))


builder <- DSI::newDSLoginBuilder()
builder$append(server = "server1", url = "dslite.server", resource = "csv", driver = "DSLiteDriver")
logindata.dslite.cnsim <- builder$build()

conns <- datashield.login(logindata.dslite.cnsim, assign=T, symbol = "csv")

datashield.assign.resource(conns, "csv", "csv")
datashield.assign.resource(conns, "pheno", "pheno")
datashield.assign.expr(conns = conns, symbol = "pheno", expr = quote(as.resource.data.frame(pheno)))
datashield.assign.expr(conns = conns, symbol = "rnaseq", expr = quote(as.resource.data.frame(csv)))

ds.dim('pheno')
ds.dim('rnaseq')

ds.createRSE('rnaseq', 'pheno', newobj.name = 'rse', annotCols = c("Symbol"))
ds.RNAseqPreproc('rse', group= 'FIR', newobj.name = 'rse.pre')

ds.dim('rse')
ds.dim('rse.pre')

ds.colnames('pheno')
ds.fvarLabels('rse.pre')

ans <- ds.limma( model = FIR ~ gender, 
                 Set="rse.pre", type.data="RNAseq", annotCols = c("Symbol"),
                 normalization = "quantile", robust = TRUE, sva=TRUE,
                 voomQualityWeights = TRUE)
o <- ans$server1 %>% filter(exp(abs(beta))>log2(1.5)) %>% arrange(desc(beta))
o





hist(ans$server1$P.Value, xlab="Raw p-value FIR effect",
     main="", las=1, cex.lab=1.5, cex.axis=1.2, col="gray")


datashield.errors()



ttFIR <- ans$server1
colnames(ttFIR)[3] <- "logFC"
colnames(ttFIR)[11] <- "Symbol"

library(plotrix)

FDRcutoff <- 0.05
fcCutoff <- 1.5
maskDE <- ttFIR$adj.P.Val < FDRcutoff & abs(ttFIR$logFC) > log2(fcCutoff)
DEgenes <- rownames(ttFIR)[maskDE]
nTopDEgenes <- min(10, length(DEgenes))
plot(ttFIR$logFC, -log10(ttFIR$P.Value), pch=".", cex=4, col=grey(0.75),
     main="", xlab=expression(paste(log[2], " Fold Change")), las=1, xlim=c(-1.2, 2),
     ylab=expression(paste(-log[10], " Raw P-value")), cex.lab=1.5, cex.axis=1.3)
abline(v=c(-log2(1.5), log2(1.5)), col=grey(0.5), lwd=1, lty=2)
abline(h=-log10(max(ttFIR$P.Value[ttFIR$adj.P.Val <= FDRcutoff])),
       col=grey(0.5), lwd=1, lty=2)
text(-log2(1.5), 0, "> 50%\nfold change", cex=0.7, pos=2)
text(log2(1.5), 0, "> 50%\nfold change", cex=0.7, pos=4)
points(ttFIR$logFC[match(DEgenes, rownames(ttFIR))],
       -log10(ttFIR$P.Value[match(DEgenes, rownames(ttFIR))]), pch=".", cex=4, col="red")
text(max(ttFIR$logFC)*0.90, -log10(max(ttFIR$P.Value[ttFIR$adj.P.Val <= FDRcutoff])),
     sprintf("%d%% FDR", 100*FDRcutoff), pos=1)
topDEgenesIdxByUpFC <- intersect(order(ttFIR$logFC, decreasing=TRUE), which(ttFIR$logFC > 0))
topDEgenesIdxByUpFC <- topDEgenesIdxByUpFC[rownames(ttFIR)[topDEgenesIdxByUpFC] %in% DEgenes]
topDEgenesIdxByDwFC <- intersect(order(ttFIR$logFC, decreasing=FALSE), which(ttFIR$logFC < 0))
topDEgenesIdxByDwFC <- topDEgenesIdxByDwFC[rownames(ttFIR)[topDEgenesIdxByDwFC] %in% DEgenes]
topDEgenesIdxByFC <- c(topDEgenesIdxByUpFC[1:nTopDEgenes], topDEgenesIdxByDwFC[1:nTopDEgenes])
topDEgenesIdxByFC <- topDEgenesIdxByFC[!is.na(topDEgenesIdxByFC)]
topDEgenesIdxByP <- order(ttFIR$adj.P.Val, decreasing=FALSE)
topDEgenesIdxByP <- topDEgenesIdxByP[rownames(ttFIR)[topDEgenesIdxByP] %in% DEgenes]
topDEgenesIdxByP <- setdiff(topDEgenesIdxByP, topDEgenesIdxByFC)
topDEgenesIdxByP <- topDEgenesIdxByP[1:min(length(topDEgenesIdxByP), nTopDEgenes)]
topDEgenesIdx <- c(topDEgenesIdxByFC, topDEgenesIdxByP)
## intelligent point label placement with the package plotrix
pos <- thigmophobe(ttFIR$logFC[topDEgenesIdx], -log10(ttFIR$P.Value[topDEgenesIdx]))
names(pos) <- sprintf("%s", ttFIR$Symbol[topDEgenesIdx])
## pos["GGH"] <- 4 ; pos["HAGH"] <- 4 ; pos["RBM38"] <- 2 ; pos["CDK1"] <- 1
## pos["LOC100008587"] <- 2 ; pos["SRGN"] <- 4 ; pos["MS4A3"] <- 1
text(ttFIR$logFC[topDEgenesIdx], -log10(ttFIR$P.Value[topDEgenesIdx]),
     labels=sprintf("%s", ttFIR$Symbol[topDEgenesIdx]), cex=0.7, pos=pos)
