#' @title Filter potential CpG outliers
#' 
#' @description Trimming methylation beta values to remove potential outliers.
#' Trimming scheme is as follows – trim values beyond the lower and upper outer #fences. These are defined by:
#' Values < 25th percentile minus 3 times IQRANDValues > 75th percentile plus 3 times IQR
#' (IQR = interquartile range)
#'
#' @param object \code{character} Name of an \code{ExpressionSet} or \code{SummarizedExperiment / RangedSummarizedExperiment}
#' on the server side to which perform the filtering
#' @param newobj.name \code{character} (default \code{NULL}) If provided, the filtered object will be assigned
#' to a variable on the server with this name. If \code{NULL}, the input variable will be overwritten.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login.
#'
#' @return
#' @export
#'
#' @examples
ds.removeOutliers <- function(object, newobj.name = NULL, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if (is.null(newobj.name)) {
    newobj.name <- object
  }
  
  cally <- paste0("removeOutliersDS(", object, ")")
  DSI::datashield.assign(datasources,  symbol = newobj.name, as.symbol(cally))
  
}