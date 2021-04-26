#'
#' @title Fits a GLM to assess association between the outcome (binary variable) and a given SNP
#' @description To be supplied
#'
#' @param snp
#' @param snps
#' @param genoData
#' @param vars
#' @param datasources
#' @param strata \code{character} Categorical variable to perform a stratified glm
#' @param family
#'
#' @return a vector with effect estimates, standard error and associated p-value
#' 
#' @author Gonzalez, JR.


glmSNP <- function(snp, snps, genoData, vars, strata, datasources, family="binomial"){
  if(!(snp %in% snps)){stop('Selected SNP [',snp,'] Not present in the genoData object')}
  i <- which(snp==snps)
  y <- vars[1]
  if(any(vars %in% strata)){
    stop("Stratification variable can't be in the model")
  }
  cally <- paste0("selSNPDS(", genoData, ", i=", i, ")")
  DSI::datashield.assign(datasources, 'dat', as.symbol(cally))
  dep <- paste(c('snp', vars[-1]), collapse="+")
  mm <- as.formula(paste(y, "~", dep))
  if(!is.null(strata)){
    levels <- ds.asFactor(paste0('dat$', strata))$all.unique.levels
    results <- lapply(levels, function(x){
      cally <- paste0("selSNPDS(", genoData, ", i=", i, ", '", strata, "', '", x, "')")
      DSI::datashield.assign(datasources, 'dat', as.symbol(cally))
      .model(mm, family, datasources, snp, x)
    })
    results <- do.call(cbind, results)
  } else{
    results <- .model(mm, family, datasources, snp, NULL)
  }
}

.model <- function(mm, family, datasources, snp, strat_var = NULL){
  mod <- try(ds.glm(mm, family=family, data='dat', 
                    datasources=datasources, viewIter = FALSE), TRUE)
  if (inherits(mod, "try-error")) {
    metrics <- data.frame(a=rep(NA,4))
    row.names(metrics) <- c("Estimate", "Std. Error", "p-value", "n")
    warning(datashield.errors())
  }
  else if (mod$errorMessage!="No errors"){
    metrics <- data.frame(a=rep(NA,4))
    row.names(metrics) <- c("Estimate", "Std. Error", "p-value", "n")
    warning(mod$errorMessage)
  }
  else {
    metrics <- as.data.frame(c(mod$coefficients[2, c(1,2,4)], mod$Nvalid))
    row.names(metrics) <- c("Estimate", "Std. Error", "p-value", "n")
  }
  if(!is.null(strat_var)){
    colnames(metrics) <- paste(snp,strat_var,sep="_")
  } else{
    colnames(metrics) <- snp
  }
  return(metrics)
}
