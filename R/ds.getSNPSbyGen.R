#' Title
#'
#' @param gdsFile 
#' @param genes 
#' @param name 
#' @param datasources 
#'
#' @return
#' @export
#'
#' @examples
ds.getSNPSbyGen <- function(gdsFile, genes, name = NULL, datasources = NULL){
 
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if(is.null(name)){
    name <- gdsFile
  }
  
  cally <- paste0("getSNPSbyGenDS(", gdsFile, ", old_assign = FALSE, '", paste0(genes, collapse = "','"), "')")
  DSI::datashield.assign.expr(datasources, name, as.symbol(cally))
  
  if(name != gdsFile){
    cally <- paste0("getSNPSbyGenDS(", gdsFile, ", old_assign = TRUE)")
    DSI::datashield.assign.expr(datasources, gdsFile, as.symbol(cally))
  }
  
}