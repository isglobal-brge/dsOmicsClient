
## ----requiredRPackages---------------------------------------------------
library(resourcer)
library(DSI)
library(DSOpal)
library(dsBaseClient)
library(dsOmicsClient)


## ----login_assign--------------------------------------------------------
builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "test.GSE66351_1", driver = "OpalDriver")
builder$append(server = "study2", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "test.GSE66351_2", driver = "OpalDriver")

logindata <- builder$build()

conns <- DSI::datashield.login(logins = logindata, assign = TRUE, 
                               symbol = "res")


## ----show_assign---------------------------------------------------------
ds.ls()


## ----coerce_df-----------------------------------------------------------
datashield.assign.expr(conns, symbol = "methyl_df", 
                       expr = quote(as.resource.data.frame(res)))
ds.class("methyl_df")


## ----dim-----------------------------------------------------------------
ds.dim("methyl_df")


## ----colnamesEset--------------------------------------------------------
varNames <- ds.colnames("methyl_df")

# CpGs
lapply(varNames, head)

# covariables
lapply(varNames, tail, n=20)


## ----inspect-------------------------------------------------------------
ds.table1D("methyl_df$diagnosis")
ds.table1D("methyl_df$Sex")


## ----glm-----------------------------------------------------------------
ds.glm(cg07363416 ~ diagnosis + Sex, data="methyl_df",
       family="binomial")


## ----close_glm, echo=FALSE-----------------------------------------------
datashield.logout(conns)


## ----login_assign_eSet---------------------------------------------------
builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "test.GSE66351_1", driver = "OpalDriver")
builder$append(server = "study2", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "test.GSE66351_2", driver = "OpalDriver")

logindata <- builder$build()

conns <- DSI::datashield.login(logins = logindata, assign = TRUE, 
                               symbol = "res")


# Assign to the original R class (e.g ExpressionSet)
datashield.assign.expr(conns, symbol = "methy", 
                       expr = quote(as.resource.object(res)))



## ----assign_es-----------------------------------------------------------
ds.class("methy")


## ----show_featureNames---------------------------------------------------
fn <- ds.featureNames("methy")
lapply(fn, head)


## ----show_phenoNames-----------------------------------------------------
ds.varLabels("methy")


## ----one_cpg-------------------------------------------------------------
ans <- ds.lmFeature(feature = "cg07363416", 
                    model = ~ diagnosis + Sex, 
                    Set = "methy",
                    datasources = conns)
ans


## ----multiple_cpg, eval=FALSE--------------------------------------------
## ans <- ds.lmFeature(model = ~ diagnosis + Sex,
##                     Set = "methy",
##                     datasources = conns,
##                     mc.cores = 20)


## ----limma_methy---------------------------------------------------------
ans.limma <- ds.limma(model = ~ diagnosis + Sex,
                      Set = "methy", 
                      datasources = conns)
ans.limma


## ----show_annot_cols-----------------------------------------------------
ds.fvarLabels("methy")


## ----remove_ans_limma, echo=FALSE----------------------------------------
ds.rm("ans.limma")


## ----limma_methy_annot---------------------------------------------------
ans.limma.annot <- ds.limma(model = ~ diagnosis + Sex,
                            Set = "methy", 
                            annotCols = c("CHR", "UCSC_RefGene_Name"),
                            datasources = conns)
ans.limma.annot


## ----remove_ans_limma_annot, echo=FALSE----------------------------------
ds.rm("ans.limma.annot")


## ----one_cpg_cellCount, eval=FALSE---------------------------------------
## ans.cell <- ds.lmFeature(feature = "cg07363416",
##                     model = ~ diagnosis + Sex,
##                     Set = "methy",
##                     cellCountsAdjust = TRUE,
##                     datasources = conns)


## ----login_assign_eSet_new, echo=FALSE-----------------------------------
datashield.logout(conns)
builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "test.GSE66351_1", driver = "OpalDriver")
builder$append(server = "study2", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "test.GSE66351_2", driver = "OpalDriver")

logindata <- builder$build()

conns <- DSI::datashield.login(logins = logindata, assign = TRUE, 
                               symbol = "res")


# Assign to the original R class (e.g ExpressionSet)
datashield.assign.expr(conns, symbol = "methy", 
                       expr = quote(as.resource.object(res)))



## ----all_cpg_sva---------------------------------------------------------
ans.sva <- ds.limma(model = ~ diagnosis + Sex, 
                    Set = "methy",
                    sva = TRUE, annotCols = c("CHR", "UCSC_RefGene_Name"))
ans.sva


## ----close_ds------------------------------------------------------------
datashield.logout(conns)


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


## ----voom_gender---------------------------------------------------------
ans.gender <- ds.limma(model =  ~ gdc_cases.demographic.gender, 
                   Set = "rse", type.data = "RNAseq", 
                   sva = FALSE)


## ----close_ds2-----------------------------------------------------------
datashield.logout(conns)


## ----resourceVCF, echo=FALSE, fig.cap="Description of how a VCF file can be added to the opal resources", out.height= '5%', fig.align='center'----
knitr::include_graphics("fig/opal_resource_VCF.png")


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

