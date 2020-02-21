##' Performing differential expression analysis using limma 
##' 
##' The function  ...
##' Outputs a matrix containing ...
##' @title Differential expression using limma
##' @param model formula indicating the condition (left side) and other covariates to be adjusted for 
##' (i.e. condition ~ covar1 + ... + covar2). The fitted model is: feature ~ condition + covar1 + ... + covarN
##' @param Set name of the DataSHIELD object to which the ExpresionSet or RangedSummarizedExperiment has been assigned
##' @param type.data optional parameter that allows the user to specify the number of CPU cores to use during 
##' @param sva logical value 
##' @param fNames ...
##' @param datasources ....
##' 
##' @export
##' @examples
##' 

ds.limma <- function(model, Set, type.data="microarray", 
                     sva=FALSE, fNames=NULL, 
                     datasources=NULL){
  
  type <- charmatch(type.data, c("microarray", "RNAseq"))
  if(is.na(type))
    stop("type.data must be 'microarray' or 'RNAseq'")
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  mt <- all.vars(model)
  variable_names <- mt[1] 
  if (length(mt)>1)
    covariable_names <- paste(mt[-1], collapse=",")
  else
    covariable_names <- NULL
  cally <- paste0("limmaDS(", Set, ",", deparse(variable_names), ",", 
                  deparse(covariable_names), ",", type, ",", sva, ",", deparse(fNames), ")")
  ans <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  class(ans) <- c("dsLimma", class(ans))
  
  return(ans)
  
}
