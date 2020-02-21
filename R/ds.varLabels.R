##' @title Retrieve phenotypic data (e.g. covariates) from eSets.
##' @description This function is similar to Bioconductor function \code{varLabels} 
##' 
##' @param object name of the DataSHIELD object to which the ExpresionSet has been assigned
##' @param datasources ....
##' 
##' @return a character vector of measured variables.
##' 
##' @export
##' @examples
##' 

ds.varLabels <- function(object, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  
  cally <- paste0("varLabelsDS(", object, ")")
  ans <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  class(ans) <- c("dsvarLabels", class(ans))
  
  return(ans)
  
}
