#' Title
#'
#' @param x 
#' @param datasources 
#'
#' @return
#' @export
#'
#' @examples
ds.genoDimensions <- function(x, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0('genoDimensionsDS(', x, ')')
  DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  
}