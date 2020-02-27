##' Performing a linear regression analysis on pooled data from multiple studies for every feature
##' 
##' The function fits a generalized linear model of a ExpressionSet for each feature (gene, CpG site, ...) 
##' in the data sets considered, using user specified condition and covariates
##' Outputs a matrix containing a beta value, standard error and p-value for each feature
##' @title Linear regression analysis of pooled data for each CpG site in study
##' @param features an optional parameter input as a vector of integer values which indicates the indices of specific 
##' features (e.g. genes, CpGs, ...) that should be analysed. If missing all features are analysed
##' @param model formula indicating the condition (left side) and other covariates to be adjusted for 
##' (i.e. condition ~ covar1 + ... + covar2). The fitted model is: feature ~ condition + covar1 + ... + covarN
##' @param Set name of the DataSHIELD object to which the ExpresionSet or RangedSummarizedExperiment has been assigned
##' @param type.p.adj multiple comparison correction method. Default 'fdr' 
##' @param cellCountsAdjust logical value which indicates whether models should be
##' adjusted for cell counts that are estimated using 'meffil.estimate.cell.counts.from.betas'
##' function from \code{meffil} package.
##' @param mc.cores optional parameter that allows the user to specify the number of CPU cores to use during 
##' parallel processing. Argument can only be > 1 when the function is run on a linux machine
##' models should be adjusted for the estimated cell counts by including the variables in the models.
##' NOTE: This assumes that the Opal pheno tables for every study include the necessary estimated cell count data 
##' originally computed when running the createOpalFiles function
##' ##' @param datasources ....
##' 
##' @export
##' @examples
##' 

ds.lmFeature <- function(features=NULL, model, Set, 
                         type.p.adj='fdr', cellCountsAdjust = FALSE,
                         mc.cores = 1, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  # Compute cell-types if cellCountsAdjust argument has been specified
  if(isTRUE(cellCountsAdjust)){
    cally <- paste0("cellCountsDS(", Set, ")")
    DSI::datashield.assign(datasources, 'cell.counts', as.symbol(cally))
    check.cell <- ds.dim("cell.counts")
    if (any(sapply(check.cell, is.null)))
      stop("There is any problem with cell-type estimation in at least one study")
  }
  else {
    DSI::datashield.assign(datasources, "cell.counts", as.symbol(NA))
  }
  
  
  mt <- as.formula(model)
  vars <- all.vars(mt)
  
  # Setting the features to loop over as the total number of 
  # features in the studies if no features are specified
  if(is.null(features)){
    cally <- paste0("featureNamesDS(", Set, ")")
    ff <- DSI::datashield.aggregate(datasources, as.symbol(cally))    
    features <- Reduce(intersect, ff)
  }
  
  ans <- t(as.data.frame(parallel::mclapply(features, lmFeature, 
                                            vars=vars,
                                            Set=Set,
                                            cellCountsAdjust=cellCountsAdjust,
                                            datasources=datasources, 
                                            mc.cores = mc.cores))) 

  if (nrow(ans) > 1) {
    #Ordering matrix by ascending p.value
    p.val <- ans[, 'p-value']
    ans <- ans[order(p.val), ]
  
    #Calculating adjusted p-values
    ans <- cbind(ans, p.adj = p.adjust(p.val,  method=type.p.adj))
  }
  
  class(ans) <- c("dsLmFeature", class(ans))
  
  return(ans)
  
}
