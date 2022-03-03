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
ds.PCA <- function(genoData, snp_subset = TRUE, standardize = TRUE, snpBlock = 20000L, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if(snp_subset){
    if(length(genoData) > 1){
      lapply(1:length(genoData), function(x){
        DSI::datashield.assign.expr(datasources, paste0("subsetGenoData_", x), 
                                    paste0("subsetGenoDS(c('", genoData[[x]], "'), 'ethnic_snps')"))
      })
      genoData <- paste0("subsetGenoData_", 1:length(genoData))
    } else {
      DSI::datashield.assign.expr(datasources, paste0("subsetGenoData_", 1), 
                                  paste0("subsetGenoDS(c('", genoData[[1]], "'), 'ethnic_snps')"))
      genoData <- "subsetGenoData_1"
    }
  }

  # TODO comprobar que els geno files segueixen el mateix ordre?? rollo que les column i files estan alineades
  # TODO Add dsBaseClient:::isDefined(datasources, object)
  # and dsBaseClient:::checkClass(datasources, object) for the genoData objects
  # TODO que es pugui agafar nomes una seleccion de SNPs (correu isglobal amb dolors)
  
  # TODO Linkage diseq filter!
  
  if(standardize){
    standard_data <- standardizeGenoData(genoData, datasources)
    toAssign <- paste0("'",paste(unlist(standard_data[,1]), collapse = ","), "'")
    DSI::datashield.assign.expr(datasources, "pca_rs", toAssign)
    toAssign <- paste0("'",paste(unlist(standard_data[,2]), collapse = ","), "'")
    DSI::datashield.assign.expr(datasources, "pca_means", toAssign)
    toAssign <- paste0("'",paste(unlist(standard_data[,3]), collapse = ","), "'")
    DSI::datashield.assign.expr(datasources, "pca_sd_hw", toAssign)
    
    cally <- paste0("PCADS(c(", paste(genoData, collapse = ", ")
                    ,"), pca_rs, pca_means, pca_sd_hw, ", snpBlock, ")")
  } else {
    cally <- paste0("PCADS(c(", paste(genoData, collapse = ", ")
                    ,"), NULL, NULL, NULL, ", snpBlock, ")")
  }
  
  res <- DSI::datashield.aggregate(datasources, cally)
  
  results <- svd(Reduce(rbind, res))$u
  
  return(list(res = list(data.frame(results)), set = genoData))

}
