##' @title Retrieve phenotypic data (e.g. covariates) from eSets.
##' @description This function is similar to Bioconductor function \code{varLabels} 
##' 
##' @param object name of the DataSHIELD object to which the ExpresionSet has been assigned
##' @param connections ....
##' 
##' @return a character vector of measured variables.
##' 
##' @export
##' @examples
##' 

ds.varLabels <- function(object, connections=NULL){
  
  if (is.null(connections)) {
    connections <- datashield.connections_find()
  }
  
  
  cally <- paste0("varLabelsDS(", object, ")")
  ans <- datashield.aggregate(connections, as.symbol(cally))
  
  class(ans) <- c("dsvarLabels", class(ans))
  
  return(ans)
  
}
