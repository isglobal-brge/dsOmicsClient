#' @title Principal Component Analysis (PCA) on SNP genotype data
#' 
#' @description PCA for genotype data on the study server
#'
#' @param gds \code{character} Name of the \code{GDS} object
#' @param prune \code{bool} (default \code{TRUE}) \code{TRUE} to prune the GDS file using \code{SNPRelate::snpgdsLDpruning}
#' @param ld.threshold Threshold for the pruning (see \code{\link{snpgdsLDpruning}})
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#'
#' @return \code{ggplot} object with PCA plot (x axis First principal component; y axis Second principal component)
#' @export

ds.PCASNPS <- function(gds, prune = TRUE, ld.threshold = 0.2, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0("PCASNPSDS(", gds, ", ", as.character(prune), ", ", ld.threshold, ")")
  res <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  class(res) <- c(class(res), "PCASNPS")
  return(list(res = res, set = gds))
}