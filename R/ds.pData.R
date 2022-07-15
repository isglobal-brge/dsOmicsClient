#' @title Retrieve information on experimental phenotypes recorded in eSet and ExpressionSet-derived classes
#' 
#' @description These generic functions access the phenotypic data (e.g., covariates) 
#' and meta-data (e.g., descriptions of covariates) associated with an experiment.
#'
#' @param x \code{character} Name of an object (on the study server), 
#' possibly derived from \link{eSet-class} or \link{AnnotatedDataFrame}
#' @param newobj.name \code{character} (default \code{NULL}) If \code{NULL}, the created object will take the 
#' name \code{"new_pData"}
#' @param datasources a list of \code{\link{DSConnection-class}} (default \code{NULL}) objects obtained after login
#'
#' @return This function does not have an output. It creates (or overwrites) a data.frame on the study server.
#' @export

ds.pData <- function(x, newobj.name = NULL, datasources = NULL){
  
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }
  
  if(is.null(newobj.name)){
    newobj.name <- "new_pData"
  }
  
  dsBaseClient:::isDefined(datasources, x)
  cls <- dsBaseClient:::checkClass(datasources, x)
  if(!any((cls %in% c("ExpressionSet")))){
    stop("'x' is not an 'ExpressionSet'")
  }
  
  cally <- paste0("pDataDS(", x, ")")
  
  DSI::datashield.assign.expr(datasources, newobj.name, as.symbol(cally))
  
}