#' @title Contingency tables of phenotypes
#' 
#' @description Wrapper of the \code{dsBaseClient::ds.table} to directly use the phenotype data stored on a 
#' \code{GenotypeData} object on the server
#'
#' @param gds \code{character} Name of the \code{GenotypeData} object on the server
#' @param variable \code{character} Phenotype variable to get the contingency table
#' @param datasources a list of \code{DSConnection-class} objects obtained after login
#'
#' @return Results of the \code{dsBaseClient::ds.table}
#' @export

ds.table_gds <- function(gds, variable, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0("extractPhenoFromGDSDS(", gds, ")")
  DSI::datashield.assign.expr(datasources, "gds_feno_table_aux", cally)
  
  ds.table(paste0("gds_feno_table_aux$", variable))
  
}
