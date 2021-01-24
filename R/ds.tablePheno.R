##' @title Retrieve the number of counts at each factor levels for metadata covariates in eSets and RSEs.
##' @description This function is similar to base function \code{table} 
##' 
##' @param object name of the DataSHIELD object to which the ExpresionSet or RangedSummarizedExperiment has been assigned
##' @param datasources ....
##' 
##' @return a numeric value
##' 
##' @export
##' @examples
##' 

ds.tablePheno <- function(object, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  

  cally <- paste0("tablePhenoDS(", object, ")")
  ans <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  class(ans) <- c("dstablePheno", class(ans))
  
  return(ans)
  
}
