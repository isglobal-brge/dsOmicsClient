##' @title Retrieve information on features from eSets.
##' @description This function is similar to Bioconductor function \code{featureData} 
##' 
##' @param object name of the DataSHIELD object to which the ExpresionSet has been assigned
##' @param connections ....
##' 
##' @return \code{featureNames} returns a (usually long!) character vector uniquely identifying each feature.
##' 
##' @export
##' @examples
##' 

ds.featureData <- function(object, connections=NULL){
  
  if (is.null(connections)) {
    connections <- datashield.connections_find()
  }
  
  
  cally <- paste0("featureDataDS(", object, ")")
  ans <- datashield.aggregate(connections, as.symbol(cally))
  
  class(ans) <- c("dsFeatureData", class(ans))
  
  return(ans)
  
}
