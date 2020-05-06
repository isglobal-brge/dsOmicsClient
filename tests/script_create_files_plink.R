setwd("C:/Juan/CREAL/GitHub/dsOmicsClient")
setwd("c:/Juan/CREAL/GitHub/brgedata/inst/extdata/")

library(SeqArray)
seqBED2GDS("brge.bed", "brge.fam", "brge.bim", out.gdsfn = "brge.gds")
# load gds
genofile <- seqOpen("brge.gds")
# convert
SeqArray::seqGDS2VCF(genofile, "brge.vcf.gz")
seqClose(genofile)


dd<-read.delim("brge.txt")
head(dd)
dd$asthma <- dd$asthma-1
write.table(dd, file="brge.txt", sep="\t", row.names=FALSE)

dd<-read.delim("obesity.txt")
dd$asthma <- geno$fam$affected
write.table(dd, file="obesity.txt", row.names=FALSE, sep="\t")


dd<-read.delim("obesity.phe")
rownames(dd) <-dd$IID
head(dd)
geno <- snpStats::read.plink("obesity")
gg <- geno$genotypes
identical(rownames(gg), dd$IID)
dd2 <- dd[rownames(gg),]
identical(rownames(gg), dd2$IID)
head(geno$fam)
dd$FID <- geno$fam$pedigree

dd <- read.table("c:/Juan/CREAL/GitHub/brgedata/inst/extdata/obesity.phe", header=TRUE)
head(dd)
dd[is.na(dd)] <- -9

dd$gender <- as.numeric(dd$gender)
dd$smoke <-  abs(as.numeric(dd$smoke) - 4)
dd$smoke[is.na(dd$smoke)] <- -9
dd$FID <- as.character(dd$FID)
dd$IID <- as.character(dd$IID)


geno <- snpStats::read.plink("brge")
dd <- read.delim("brge.phe")


write.table(dd, file="c:/Juan/CREAL/GitHub/brgedata/inst/extdata/brge.phe", 
            row.names=FALSE, quote=FALSE)

