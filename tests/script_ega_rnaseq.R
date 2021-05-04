library(resourcer)
library(DSLite)
library(dsOmics)
library(dsOmicsClient)
library(dsBaseClient)
library(SummarizedExperiment)


csv <- newResource(name = "counts", 
                      url = "file:c:/Juan/CREAL/ATHLETE/DataSHIELD/EGA/RNAseq/RawGeneExpressionProfiles.csv", 
                      format = "csv")
pheno <- newResource(name = "pheno", 
                     url = "file:c:/Juan/CREAL/ATHLETE/DataSHIELD/EGA/RNAseq/rnaseq_rob_feno.csv", 
                     format = "csv")


pkgs <- c("dsBase", "resourcer", "dsOmics", "SummarizedExperiment", "GenomicRanges", "edgeR")
dslite.server <- newDSLiteServer(resources = list(csv=csv, pheno=pheno),
                                 config = DSLite::defaultDSConfiguration(include=pkgs))

dslite.server$aggregateMethod("filteredByExprDS", "filteredByExprDS")


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

ds.createRSE('rnaseq', 'pheno', newobj.name = 'rse')
ds.filterByExpr('rse')

ds.dim('rse')
ds.dim('rse.filter')

ds.limma( model = FIR ~ Sex, 
          Set="rse.filter", type.data="RNAseq")
datashield.errors()
