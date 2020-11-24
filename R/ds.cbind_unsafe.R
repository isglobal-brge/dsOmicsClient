#' Title
#'
#' @param obj1 
#' @param obj2 
#' @param name 
#' @param datasources 
#'
#' @return
#' @export
#'
#' @examples
ds.cbind_unsafe <- function(obj1, obj2, name, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0("cbind_unsafeDS(", obj1, ", ", obj2, ")")
  DSI::datashield.assign.expr(datasources, name, as.symbol(cally))
  
}