#' @title Retrieve SNPs from  
#'
#' @param genoData \code{character} Object of class \code{GenotypeData} or 
#' \code{GdsGenotypeReader} on the server
#' @param datasources a list of \code{\link{DSConnection-class}} (default \code{NULL}) objects obtained after login
#'
#' @return Character vector with the SNP names
#' @export

ds.getSNPs <- function(genoData, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0("getVariable(", genoData, ", 'snp.rs.id')")
  snps <- unlist(DSI::datashield.aggregate(datasources,  cally))
  return(snps)
}