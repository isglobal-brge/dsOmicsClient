#' @title Filter Genes By Expression Level
#' 
#' @description Determine which genes have sufficiently large counts to be retained in a statistical analysis. This is a similar function to \code{edgeR::filterByExpr}
#'
#' @param object \code{character} Name of the \code{eSet}, \code{RangedSummarizedExperiment} on the server
#' @param newobj.name \code{character} (default \code{NULL}) Name to be assigned 
#' on the server to the RangedSummarizedExperiment. If \code{NULL} the created \code{RangedSummarizedExperiment} will 
#' be assigned to a variable named \code{'rse.filter'}
#' @param datasources a list of \code{\link{DSConnection-class}} (default \code{NULL}) objects obtained after login
#'
#' @return This function does not have an output. It creates (or overwrites) a data frame on the study server.
#' @export

ds.filterByExpr <- function(object,  newobj.name = NULL, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if(is.null(newobj.name)){
    newobj.name <- 'rse.filter'
  }
  
  dsBaseClient:::isDefined(datasources, object)
  cls <- dsBaseClient:::checkClass(datasources, object)
  if(!cls %in% c("RangedSummarizedExperiment", "SummarizedExperiment", "ExpressionSet")){
    stop("[",object,"] is not a 'data.frame'")
  }

  cally <- paste0("filterByExprDS(", object, " )")
  DSI::datashield.assign.expr(datasources, newobj.name, cally)
}