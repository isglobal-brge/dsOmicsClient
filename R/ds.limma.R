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
##' @param connections ....
##' 
##' @export
##' @examples
##' 

ds.limma <- function(model, Set, type.data="microarray", 
                         sva=FALSE, connections=NULL){
  
  if (is.null(connections)) {
    connections <- datashield.connections_find()
  }
  
  mt <- as.formula(model)
  vars <- all.vars(mt)
  
  
  # Compute SVA if sva argument has been specified
  if(isTRUE(sva)){
    cally <- paste0("svaDS(", deparse(vars), ",", Set, ")")
    datashield.assign(connections, 'sva', as.symbol(cally))
    check.sva <- ds.dim("sva")
    if (any(sapply(check.sva, is.null)))
      stop("There is any problem with SVA estimation in at least one study")
  }
  else {
    datashield.assign(connections, "sva", as.symbol(NA))
  }
  
  
  # vars[1] is to avoid non-disclosive access since only
  # the design part is necessary (DS does not allow '~ A + B')
  mod <- paste(vars[1], " ~", paste(vars, collapse="+"))
  
  cally <- paste0("limmaDS(", mod, ",", Set, ",", sva, ")")
  ans <- datashield.aggregate(connections, as.symbol(cally))
  
  class(ans) <- c("dsLimma", class(ans))
  
  return(ans)
  
}
