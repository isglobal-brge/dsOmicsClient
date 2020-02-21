#'
#' @title Fits a GLM to assess association between the outcome (binary variable) and a given SNP
#' @description To be supplied
#' @param snp
#' @param snps
#' @param genoData
#' @param vars
#' @param datasources
#' @param family
#' @return a vector with effect estimates, standard error and associated p-value
#' 
#' @export
#' @author Gonzalez, JR.


glmSNP <- function(snp, snps, genoData, vars, datasources, family="binomial"){
  i <- which(snp==snps)
  cally <- paste0("selSNPDS(", genoData, ", i=", i,")")
  DSI::datashield.assign(datasources, 'dat', as.symbol(cally))
  
  y <- vars[1]
  dep <- paste(c("snp", vars[-1]), collapse="+")
  mm <- as.formula(paste(y, "~", dep))
  mod <- try(ds.glm(mm, family=family, data='dat', 
                    datasources=datasources, viewIter = FALSE), TRUE)
  if (inherits(mod, "try-error")) {
    metrics <- data.frame(a=rep(NA,3))
    row.names(metrics) <- c("Estimate", "Std. Error", "p-value")  
  }
  else if (mod$errorMessage!="No errors"){
    metrics <- data.frame(a=rep(NA,3))
    row.names(metrics) <- c("Estimate", "Std. Error", "p-value")  
  }
  else {
    metrics <- as.data.frame(mod$coefficients[2, c(1,2,4)])
  }
  colnames(metrics) <- snp
  return(metrics)
}
