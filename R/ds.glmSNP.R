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
##' @param data name of the DataSHIELD object to which the genotype (snpMatrix) and phenotypic data (data.frame) has been assigned
##' @param datasources Opal object or list of opal objects denoting the opal server(s) information
##' @param mc.cores optional parameter that allows the user to specify the number of CPU cores to use during 
##' parallel processing. Argument can only be > 1 when the function is run on a linux machine
##' @export
##' @examples
##' 

ds.glmSNP <- function(snps.fit=NULL, model, gds, covars, connections,
                      maf=0.05, missing.rate=0.9, 
                      type.p.adj='fdr', mc.cores = 1){
  
  
  mt <- as.formula(model)
  vars <- all.vars(mt)
  
  # Check whether data are in the same order
  cally <- paste0("readGdsnDS(indexGdsnDS(", gds, ", 'sample.id'))")
  ids <- datashield.aggregate(connections,  cally)
#  if (!identical(ids, vars[,1]))
#    stop("There VCF ids are not identical to those in the first column of 'covars' ")
  
  
  # Get snp rs's
  cally <- paste0("readGdsnDS(indexGdsnDS(", gds, ", 'snp.rs.id'))")
  snps <- unlist(datashield.aggregate(connections,  cally))
  
  # Setting the SNPs to loop over if no SNPs are specified.
  # They are based on those that pass QC
  if(is.null(snps.fit)){
    cally <- paste0("snpgdsSelectSNPDS(", gds, ",maf=", maf, 
                    ",missing.rate=", missing.rate, ")")
    
    snps.id <- unlist(datashield.aggregate(connections, cally))
    snps.fit <- snps[snps.id]
  }
  
  
  ans <- t(as.data.frame(parallel::mclapply(snps.fit, glmSNPDS, 
                                            snps=snps, 
                                            covars=covars,
                                            vars=vars,
                                            gds=gds, 
                                            connections=connections,
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
