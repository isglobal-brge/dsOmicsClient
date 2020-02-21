rm(list=ls())
library(resourcer)
library(DSOpal)
library(dsBaseClient)
library(DSLite)
library(SNPRelate)

library(dsOmics)
library(dsOmicsClient)
# source("tests/getClientVCF.R")


################   GWAS data analysis

dslite.server <- newDSLiteServer(resources = list(
  vcf = newResource(name = "vcf", url = "https://raw.githubusercontent.com/isglobal-brge/brgeUtils/master/data/1KG_phase3_subset_chr22.vcf.gz", format = "VCF2GDS"),
  covars = newResource(name = "covars", url = "https://raw.githubusercontent.com/isglobal-brge/brgeUtils/master/data/covars_1000G.csv", format = "csv")
))


dslite.server <- newDSLiteServer(resources = list(
  vcf = newResource(name = "vcf", url = "https://raw.githubusercontent.com/isglobal-brge/brgedata/master/inst/extdata/obesity.vcf.gz", format = "VCF2GDS"),
  covars = newResource(name = "covars", url = "https://raw.githubusercontent.com/isglobal-brge/brgedata/master/inst/extdata/obesity.txt", format = "tsv")
))



resourcer::registerResourceResolver(GDSFileResourceResolver$new())


builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "dslite.server", resource = "vcf", driver = "DSLiteDriver")
logindata <- builder$build()
logindata
# login and assign resources
conns <- datashield.login(logins = logindata, assign = TRUE, symbol = "vcf.res")

datashield.assign.expr(conns, symbol = "gds", expr = quote(as.resource.object(vcf.res)))

# assign covars resource
datashield.assign.resource(conns, symbol = "covars.res", resource = list(study1 = "covars"))
datashield.assign.expr(conns, symbol = "covars", expr = quote(as.resource.data.frame(covars.res)))

ds.ls()
ds.class("gds")
ds.class("covars")


# inspect data in the first study (DSLite only!)
dslite.server$getSessionData(conns[[1]]@sid, "gds")
dslite.server$getSessionData(conns[[1]]@sid, "covars")


dslite.server$aggregateMethod("nsnpDS", "GWASTools::nsnp")
dslite.server$aggregateMethod("nscanDS", "GWASTools::nscan")
dslite.server$aggregateMethod("getVariable", "GWASTools::getVariable")


dslite.server$aggregateMethod("indexGdsnDS", "gdsfmt::index.gdsn")
dslite.server$aggregateMethod("readGdsnDS", "gdsfmt::read.gdsn")
dslite.server$aggregateMethod("snpgdsGetGenoDS", "SNPRelate::snpgdsGetGeno")
dslite.server$aggregateMethod("snpgdsSelectSNPDS", "SNPRelate::snpgdsSelectSNP")
dslite.server$aggregateMethod("snpgdsOpenDS", "SNPRelate::snpgdsOpen")
dslite.server$aggregateMethod("GWASDS", GWASDS)



dslite.server$assignMethod("selSNPDS", selSNPDS)

dslite.server$assignMethod("GenotypeDataDS", GenotypeDataDS)
dslite.server$assignMethod("selSNPDS", selSNPDS)

x <- 'gds'
xx <- 'covars'
columnId <- 1


cally <- paste0("GenotypeDataDS(", x, "," , xx, ",", columnId, ")")
cally
DSI::datashield.assign(conns,  symbol = "gdsData", as.symbol(cally))
ds.GenotypeData('gds', 'covars', 1, newobj.name = 'gds.Data' , conns)

snps.fit <- c("rs4970516", "rs2504786", "rs2474293", "rs11247693")
ds.glmSNP(snps.fit, model = obese ~ gender + age, genoData = 'gds.Data')

y <- "obese"
x <-  paste(c("age", "country"), collapse=",")
cally <- paste0("GWASDS(", 'gdsData', "," , deparse("obese"), ",", deparse(x), ")")
cally
DSI::datashield.aggregate(conns,  as.symbol(cally))

ds.GWAS('gdsData', model=obese~age+country)



cally <- paste0("nsnpDS(", gds, ")")
cally <- paste0("nscanDS(", gds, ")")
cally
DSI::datashield.aggregate(conns,  cally)
ds.dim("covars")

maf <- 0.05
missing.rate <- 0.95

cally <- paste0("readGdsnDS(indexGdsnDS(", gds, ", 'sample.id'))")
ids <- DSI::datashield.aggregate(conns,  cally)

cally <- paste0("readGdsnDS(indexGdsnDS(", gds, ", 'snp.rs.id'))")
snps <- unlist(DSI::datashield.aggregate(conns,  cally))
names(snps) <- 1:length(snps)
cally <- paste0("snpgdsSelectSNPDS(", gds, ",maf=", maf, 
                ",missing.rate=", missing.rate, ")")



a<-ids$study1[1:3]
sample.id <- paste(a,collapse=",")
cally <- paste0("snpgdsSelectSNPDS(", gds,  ",", deparse(sample.id) , ")")
cally
DSI::datashield.aggregate(conns,  cally)



snps.id <- unlist(DSI::datashield.aggregate(connections, cally))
snps.fit <- snps[snps.id]
names(snps.fit) <- 1:length(snps.fit)

#out <- ds.glmSNP(model="group~sex", gds='gds', covars='covars', connections=conns)
out1 <- ds.glmSNP("rs115743375", model="group~sex", gds='gds', covars='covars', connections=conns)
out2 <- ds.glmSNP("rs115743375", model="group~1", gds='gds', covars='covars', connections=conns)

cally <- "readGdsnDS(indexGdsnDS(gds, 'sample.id'))"
ids <- datashield.aggregate(conns,  cally)
if (!identical(ids,vars[,1]))
 stop("There VCF ids are not identical to those in the first column of 'covars' ")


i <- which(snp==snps)
cally <- paste0("snpgdsGetGenoDS(gds, snp.id=",1,")")
datashield.aggregate(conns, as.symbol(cally))




cally <- "selSNPDS(gds, i=1, covars, c('group', 'sex'))"
datashield.assign(conns, 'dat', as.symbol(cally))
ds.ls()


dslite.server$closeSession("conns")

# new package GENESIS

SeqArray::seqBED2GDS("c:/juan/CREAL/GitHub/brgedata/inst/extdata/obesity.bed",
                     "c:/juan/CREAL/GitHub/brgedata/inst/extdata/obesity.fam",
                     "c:/juan/CREAL/GitHub/brgedata/inst/extdata/obesity.bim",
                     out="c:/tmp/obesit.gds")
SeqArray::seqGDS2VCF("c:/tmp/obesit.gds", "c:/tmp/obesity.vcf.gz")


SNPRelate::snpgdsVCF2GDS("c:/tmp/1KG_phase3_subset_chr22.vcf.gz",
                          out.fn="c:/tmp/1KG_phase3_subset_chr22.gds")
gds.f <- SNPRelate::snpgdsOpen("c:/tmp/1KG_phase3_subset_chr22.gds", readonly = FALSE)
sample.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds.f, "sample.id"))
dd <- read.delim("c:/tmp/1000G_info.txt", as.is=TRUE)
covars <- dd[dd$Sample%in%sample.id,]
set.seed(1234)
covars$casco <- sample(c(0,1), nrow(covars), replace=TRUE) 
readr::write_csv(covars, "c:/juan/CREAL/GitHub/brgeUtils/data/covars_1000G.csv")



SNPRelate::snpgdsVCF2GDS("c:/tmp/obesity.vcf.gz",
                         out.fn="c:/tmp/obesity.gds")
covars <- read.delim("c:/juan/CREAL/GitHub/brgedata/inst/extdata/obesity.txt", as.is=TRUE)
rownames(covars) <- covars$id
covars <- covars[getScanID(g),]

library(GWASTools)
covars <- read.delim("c:/juan/CREAL/GitHub/brgedata/inst/extdata/obesity.txt", as.is=TRUE)  
names(covars)[1] <- "scanID"
scanAnnot <- ScanAnnotationDataFrame(covars)
g <- GdsGenotypeReader("c:/tmp/obesity.gds")
geno <- GenotypeData(g, scanAnnot = scanAnnot)
close(g)

getGenotype(geno, char=TRUE)

getGenotypeSelection(g, snp=1)




scan.exclude=NULL
res <- assocRegression(geno,
                       outcome="obese",
                       model.type="logistic",
                       covar=c("age", "gender"),
                       scan.exclude=scan.exclude)


nullmod <- GENESIS::fitNullModel(genoData, outcome = "obese", 
                        covars = c("age", "gender"), 
                        family = binomial)
genoIterator <- GenotypeBlockIterator(genoData, snpBlock=5000)
assoc <- GENESIS::assocTestSingle(genoIterator, null.model = nullmod)
assoc$rs<-getVariable(genoData, "snp.rs.id")[ans$variant.id]
ans <- assoc %>% tibble::as_tibble() %>%
  arrange(Score.pval) %>% select(variant.id, rs, everything())

                    

# .... R
gds.f<-SNPRelate::snpgdsOpen("c:/tmp/gds.f.gds", readonly = FALSE)
covars <- read.delim("c:/juan/CREAL/GitHub/brgeUtils/data/covars_1000G.txt", sep="")

sample.id <- read.gdsn(index.gdsn(gds.f, "sample.id"))
snps.id <- read.gdsn(index.gdsn(gds.f, "snp.id"))
snps <- read.gdsn(index.gdsn(gds.f, "snp.rs.id"))
g <- as.numeric(snpgdsGetGeno(gds.f, snp.id = snps.id[1]))
g



datashield.logout(conns)



##############  First try with Yannick's function

source("tests/getClientVCF.R")


# register resolver so that magic happens

resourcer::registerResourceResolver(GDSFileResourceResolver$new())



# use it (method and snpfirstdim are optional)

res <- resourcer::newResource(url = "https://raw.githubusercontent.com/isglobal-brge/brgeUtils/master/data/1000G.vcf", format = "VCF2GDS")

client <- resourcer::newResourceClient(res)


gds.f <- client$getValue()

sample.id <- read.gdsn(index.gdsn(gds.f, "sample.id"))
snps.id <- read.gdsn(index.gdsn(gds.f, "snp.id"))
snps <- read.gdsn(index.gdsn(gds.f, "snp.rs.id"))
g <- as.numeric(snpgdsGetGeno(gds.f, snp.id = snps.id[1]))
g

client$close()

