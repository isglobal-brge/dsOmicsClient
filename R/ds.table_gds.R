#' Title
#'
#' @param gds 
#' @param variable 
#' @param datasources 
#'
#' @return
#' @export
#'
#' @examples
ds.table_gds <- function(gds, variable, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0("extractPhenoFromGDSDS(", gds, ")")
  DSI::datashield.assign.expr(datasources, "gds_feno_table_aux", cally)
  
  ds.table(paste0("gds_feno_table_aux$", variable))
  
}
