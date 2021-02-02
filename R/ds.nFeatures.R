##' @title Retrieve the number of features for eSets and RSEs.
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

ds.nFeatures <- function(object, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  

  cally <- paste0("nFeaturesDS(", object, ")")
  ans <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  class(ans) <- c("dsnFeatures", class(ans))
  
  return(ans)
  
}
