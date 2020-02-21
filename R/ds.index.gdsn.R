##' @title Retrieve phenotypic data (e.g. covariates) from eSets.
##' @description This function is similar to function \code{index.gdsn}  in the 'gdsfmt' package 
##' 
##' @param x an object of class
##' @param datasources ....
##' 
##' @return a character vector of measured variables.
##' 
##' @export
##' @examples
##' 

ds.index.gdsn <- function(x=NULL, datasources=NULL){
  if (is.null(datasources)) {
    datasources <- datashield.connections_find()
  }
  if (is.null(x)) {
      stop("Please provide the name of the input object!", 
           call. = FALSE)
  }
  defined <- isDefined(datasources, x)
  cally <- paste0("class(", x, ")")
  output <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  return(output)
}
  
  
