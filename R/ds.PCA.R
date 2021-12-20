#' Title
#'
#' @param genoData 
#' @param standardize 
#' @param snpBlock 
#' @param datasources 
#'
#' @return
#' @export
#'
#' @examples
ds.PCA <- function(genoData, standardize = TRUE, snpBlock = 20000L, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  # TODO comprobar que els geno files segueixen el mateix ordre?? rollo que les column i files estan alineades
  # TODO Add dsBaseClient:::isDefined(datasources, object)
  # and dsBaseClient:::checkClass(datasources, object) for the genoData objects
  
  if(standardize){
    standard_data <- standardizeGenoData(genoData, snpBlock, datasources)
  }
  
  browser()
  # TODO aqui començar a afegir tot el codi que tinc al servidor de isglobal 
  # i aplicar els coeficients de standaritzacio si fa falta al llegir el geno data
  # block by block!!!
  
  
  # TODO remove normalized genotype data on the servers
  
}
