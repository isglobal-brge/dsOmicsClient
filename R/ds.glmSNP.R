##' Performing a logistic regression analysis on pooled data from multiple studies for every SNP
##' 
##' Function fits a generalized linear model to gentic data for each SNP
##' in the data sets considered, using user specified outcome and phenotypic variabes as covariates
##' Outputs a matrix containing a beta value, standard error and p-value for each SNP
##' @title Logistic regression analysis of pooled data for each SNP site in study
##' @param snps.fit an optional parameter input as a character vector of
##' SNPs (rs numbers) that should be analysed. If missing all SNPs are analysed
##' @param model list of phenotypic variables to use as covariates in the regression analysis in the form: 
##' "outcome ~ covar1 + covar2 +  ... + covarN"
##' @param genoData name of the DataSHIELD object to which the genotype (snpMatrix) and phenotypic data (data.frame) has been assigned
##' @param datasources Opal object or list of opal objects denoting the opal server(s) information
##' @param mc.cores optional parameter that allows the user to specify the number of CPU cores to use during 
##' parallel processing. Argument can only be > 1 when the function is run on a linux machine
##' @export
##' @examples
##' 

ds.glmSNP <- function(snps.fit=NULL, model, genoData, datasources=NULL,
                      type.p.adj='fdr', mc.cores = 1){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  mt <- as.formula(model)
  vars <- all.vars(mt)
  
# check vars are present in covars
  
  
  # Get snp rs's
  cally <- paste0("getVariable(", genoData, ", 'snp.rs.id')")
  snps <- unlist(DSI::datashield.aggregate(datasources,  cally))
  
  # Setting the SNPs to loop over if no SNPs are specified.
  # They are based on those that pass QC
  if(is.null(snps.fit)){
    snps.fit <- snps
  }

  
  ans <- t(as.data.frame(parallel::mclapply(snps.fit, glmSNP, 
                                            snps=snps, 
                                            vars=vars,
                                            genoData=genoData, 
                                            datasources=datasources,
                                            mc.cores = mc.cores))) 
  
  if (nrow(ans) > 1) {
  #Ordering matrix by ascending p.value
  ans <- ans[order(ans[, 'p-value']), ]
  
  #Calculating adjusted p-values
  ans <- cbind(ans, p.adj = p.adjust(ans[,'p-value'],  method=type.p.adj))
  }
  
  
  class(ans) <- c("dsGlmSNP", class(ans))
  
  return(ans)
  
}
