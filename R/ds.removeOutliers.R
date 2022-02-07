#' @title Filter potential CpG outliers
#' 
#' @description Trimming methylation beta values to remove potential outliers.
#' Trimming scheme is as follows: trim values beyond the lower and upper outer fences. These are defined by 
#' the `pct` argument.
#'
#' @param object \code{character} Name of an \code{ExpressionSet} or \code{SummarizedExperiment / RangedSummarizedExperiment}
#' on the server side to which perform the filtering
#' @param pct \code{numeric} (default \code{0.125}) Tail and head quantiles that will be considered as 
#' outliers; for example: a \code{pct=0.125} will use the function 
#' \code{matrixStats::rowQuantiles(probs = c(0.125, 1-0.125))} to detect outliers.
#' @param newobj.name \code{character} (default \code{NULL}) If provided, the filtered object will be assigned
#' to a variable on the server with this name. If \code{NULL}, the input variable will be overwritten.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login.
#'
#' @return This function does not have an output. It creates (or overwrites) an \code{ExpressionSet} on the study server.
#' @export
#'
#' @examples
ds.removeOutliers <- function(object, pct = 0.125, newobj.name = NULL, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if (is.null(newobj.name)) {
    newobj.name <- object
  }
  
  cally <- paste0("removeOutliersDS(", object, ",", pct, ")")
  DSI::datashield.assign(datasources,  symbol = newobj.name, as.symbol(cally))
  
}