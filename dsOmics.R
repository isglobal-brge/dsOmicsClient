
## ----requiredRPackages---------------------------------------------------
library(resourcer)
library(DSI)
library(DSOpal)
library(dsBaseClient)
library(dsOmicsClient)


#
# Example 1: Methylation data analysis (ExpressionSet)
#

load("c:/juan/CREAL/GitHub/dsOmicsClient/data/GSE66351.Rdata")
gse66351.sel


## ----login_assign--------------------------------------------------------
builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "test.GSE66351", driver = "OpalDriver")
builder$append(server = "study2", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "test.GSE80970", driver = "OpalDriver")

logindata <- builder$build()

conns <- DSI::datashield.login(logins = logindata, assign = TRUE, 
                               symbol = "res")


## ----show_assign---------------------------------------------------------
ds.ls()


## ----coerce_df-----------------------------------------------------------
datashield.assign.expr(conns, symbol = "methyl_df", 
                       expr = quote(as.resource.data.frame(res)))
ds.class("methyl_df")


## ----inspect-------------------------------------------------------------
ds.summary("methyl_df$casecon")
ds.summary("methyl_df$cg07363416")


## ----glm-----------------------------------------------------------------
ds.glm(cg07363416 ~ casecon + Sex, data="methyl_df",
       family="binomial")


## ----assign_es-----------------------------------------------------------
datashield.assign.expr(conns, symbol = "methy", 
                       expr = quote(as.resource.object(res)))
ds.class("methy")


## ----show_featureNames---------------------------------------------------
fn <- ds.featureNames("methy")
lapply(fn, head)


## ----show_phenoNames-----------------------------------------------------
ds.varLabels("methy")


## ----one_cpg-------------------------------------------------------------
ans <- ds.lmFeature(feature="cg07363416", 
                    model=casecon~Sex, 
                    eSet="methy",
                    datasources=conns)
ans


## ----multiple_cpg, eval=FALSE--------------------------------------------
ans <- ds.lmFeature(model = casecon~Sex, 
                    eSet = "methy",
                    datasources = conns,
                    mc.cores = 20)


## ----limma_methy---------------------------------------------------------
ans.limma <- ds.limma(model = ~ casecon + Sex,
                      Set = "methy", 
                      datasources = conns)
ans.limma


## ----one_cpg_cellCount, error=TRUE---------------------------------------
ans.cell <- ds.lmFeature(feature = "cg07363416", 
                    model = casecon ~ Sex, 
                    eSet = "methy", 
                    cellCountsAdjust = TRUE,
                    datasources = conns)


## ----all_cpg_sva, error=TRUE---------------------------------------------
ans.sva <- ds.limma(model = casecon ~ Sex, 
                    Set = "methy",
                    sva = TRUE)


## ----close_ds------------------------------------------------------------
datashield.logout(conns)


#
# Example 2: RangedSummarizedExperiment (TCGA - Recount)
#



## ----pipeline_gene_expr--------------------------------------------------
builder <- newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "test.tcga_liver", driver = "OpalDriver")

logindata <- builder$build()

conns <- datashield.login(logins = logindata, assign = TRUE, 
                          symbol = "res")


## ----get_rse-------------------------------------------------------------
datashield.assign.expr(conns, symbol = "rse", 
                       expr = quote(as.resource.object(res)))
ds.class("rse")


## ----dim_rse-------------------------------------------------------------
ds.dim("rse")


## ----name_feature_rse----------------------------------------------------
name.features <- ds.featureNames("rse")
lapply(name.features, head)


## ----name_covar_rse------------------------------------------------------
name.vars <- ds.featureData("rse")
lapply(name.vars, head, n=15)


## ----table_gender--------------------------------------------------------
ds.table1D("rse$gdc_cases.demographic.gender")


## ----voom_gender, eval=FALSE---------------------------------------------
ans.gender <- ds.limma(model =  ~ gdc_cases.demographic.gender, 
                   Set = "rse", type.data = "RNAseq", 
                   sva = FALSE)


## ----close_ds2-----------------------------------------------------------
datashield.logout(conns)

#
# Example 3: GWAS data analysis
#

covars <- read.delim("c:/juan/CREAL/GitHub/brgedata/inst/extdata/obesity.txt", as.is=TRUE)  
names(covars)[1] <- "scanID"
head(covars)
scanAnnot <- ScanAnnotationDataFrame(covars)
g <- GdsGenotypeReader("c:/tmp/obesity.gds")
g
geno <- GenotypeData(g, scanAnnot = scanAnnot)
geno
close(g)



## ----add_resources_vcf---------------------------------------------------
builder <- newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org",
               user = "dsuser", password = "password",
               resource = "test.obesity_vcf", driver = "OpalDriver")
logindata <- builder$build()

conns <- datashield.login(logins = logindata, assign = TRUE,
                          symbol = "res")


## ----assign_vcf----------------------------------------------------------

datashield.assign.resource(conns, symbol = "vcf.res", 
                           resource = list(study1 = "test.obesity_vcf"))
datashield.assign.expr(conns, symbol = "gds", 
                       expr = quote(as.resource.object(vcf.res)))


datashield.assign.resource(conns, symbol = "covars.res", 
                           resource = list(study1 = "test.obesity"))
datashield.assign.expr(conns, symbol = "covars", 
                       expr = quote(as.resource.data.frame(covars.res)))


## ----ls_vcf--------------------------------------------------------------
ds.ls()


## ----show_covars---------------------------------------------------------
ds.colnames("covars")


## ----show_group----------------------------------------------------------
ds.table1D("covars$obese")


## ----createGenoData------------------------------------------------------
ds.GenotypeData(x='gds', covars = 'covars', columnId = 1, newobj.name = 'gds.Data')


## ----snp_analysis--------------------------------------------------------
ds.glmSNP(snps.fit = "rs11247693", model = obese ~ gender + age, genoData='gds.Data')


## ----GWAS----------------------------------------------------------------
ds.GWAS('gds.Data', model=obese~age+country)


## ----close_conns3--------------------------------------------------------
datashield.logout(conns)

