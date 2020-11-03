#' @title Get names of chromosomes
#' 
#' @description Get the names of the chromosomes of a GenotypeData object on the server side
#'
#' @param genoData \code{character} Name of the \code{\link{GenotypeData}} object on the server
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#'
#' @return A \code{list} with: \cr
#' -autosomes: integer codes for the autosomes \cr
#' -Xchr: integer code for the X chromosome \cr
#' -pseudoautosomalRegionXY: integer code for the pseudoautosomal region \cr
#' -Ychr: integer code for the Y chromosome \cr
#' -mitochondrial: integer code for mitochondrial SNPs \cr
#' 
#' @export

ds.getChromosomeNames <- function(genoData, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0("getChromosomeNamesDS(", genoData, ")")
  res <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  return(res)
  
}