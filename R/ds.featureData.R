##' @title Retrieve information on features from eSets and RSEs.
##' @description This function is similar to Bioconductor function \code{featureData} 
##' 
##' @param object name of the DataSHIELD object to which the ExpresionSet has been assigned
##' @param datasources ....
##' 
##' @return \code{featureNames} returns a (usually long!) character vector uniquely identifying each feature.
##' 
##' @export
##' @examples
##' 

ds.featureData <- function(object, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  
  cally <- paste0("featureDataDS(", object, ")")
  ans <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  class(ans) <- c("dsFeatureData", class(ans))
  
  return(ans)
  
}
