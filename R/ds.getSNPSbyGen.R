#' @title Subset GDS with gene(s)
#' 
#' @description Subset a server side GDS with the SNPs from selected gene(s)
#'
#' @param gdsFile \code{character} Name of the GDS on the server
#' @param genes \code{character vector} Name of the gene(s) to select (HGCN nomenclature)
#' @param name \code{character} (default \code{NULL}) Name of the subsetted GDS, if \code{NULL} 
#' entry file will be over-written.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#'
#' @return This function does not have an output. It creates (or overwrites) a GDS object on the study server.
#' @export

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