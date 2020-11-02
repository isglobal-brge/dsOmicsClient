#' Title
#'
#' @param gds 
#' @param ld.threshold 
#' @param datasources 
#'
#' @return
#' @export
#'
#' @examples
ds.PCASNPS <- function(gds, prune = TRUE, ld.threshold = 0.2, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0("PCASNPSDS(", gds, ", ", as.character(prune), ", ", ld.threshold, ")")
  res <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
}