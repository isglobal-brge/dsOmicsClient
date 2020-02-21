##' @title Retrieve feature names from eSets
##' @description This function is similar to Bioconductor function \code{featureNames} 
##' 
##' @param object name of the DataSHIELD object to which the ExpresionSet has been assigned
##' @param datasources ....
##' 
##' @return \code{featureNames} returns a (usually long!) character vector uniquely identifying each feature.
##' 
##' @export
##' @examples
##' 

ds.featureNames <- function(object, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  

  cally <- paste0("featureNamesDS(", object, ")")
  ans <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  class(ans) <- c("dsFeatureNames", class(ans))
  
  return(ans)
  
}
