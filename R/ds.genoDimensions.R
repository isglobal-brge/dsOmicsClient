#' @title Get main dimensions of Genotype data
#' 
#' @description Get the number of SNPs, number of scans and number of chromosomes on the genotype file
#'
#' @param x \code{character} Name of the \code{GenotypeData} or \code{GdsGenotypeReader} object on the server
#' @param datasources a list of \code{DSConnection-class} objects obtained after login
#'
#' @return \code{list} with results
#' @export

ds.genoDimensions <- function(x, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0('genoDimensionsDS(', x, ')')
  DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  
}