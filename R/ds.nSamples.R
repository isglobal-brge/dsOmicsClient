##' @title Retrieve the number of samples for eSets and RSEs.
##' @description This function is similar to Bioconductor function \code{dim} 
##' 
##' @param object name of the DataSHIELD object to which the ExpresionSet or RangedSummarizedExperiment has been assigned
##' @param datasources ....
##' 
##' @return a numeric value
##' 
##' @export
##' @examples
##' 

ds.nSamples <- function(object, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  

  cally <- paste0("nSamplesDS(", object, ")")
  ans <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  class(ans) <- c("dsnSamples", class(ans))
  
  return(ans)
  
}
