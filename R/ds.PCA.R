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
ds.PCA <- function(genoData, snp_subset = TRUE, standardize = TRUE, snpBlock = 20000L, ncomp = 5L, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  # Aux variable
  genoData_og <- genoData
  classes <- ds.class(genoData_og, datasources)
  if(!all(unlist(classes) == "GenotypeData")){
    stop("[genoData] objects must all be `GenotypeData` objects")
  }
  
  # Temp name
  tfile <- substring(tempfile(pattern = "", tmpdir = ""), 2)
  
  if(snp_subset){
    if(length(genoData) > 1){
      assigneds <- unlist(lapply(1:length(genoData), function(x){
        tryCatch({
          DSI::datashield.assign.expr(datasources, paste0("subsetGenoData_", x, tfile), 
                                      paste0("subsetGenoDS(c('", genoData[[x]], "'), 'ethnic_snps')"))
          x
        }, error = function(w){
          message('[',genoData[[x]], '] Geno file could not be subsetted because of disclosue risk')
        })
      }))
      genoData <- paste0("subsetGenoData_", assigneds, tfile)
    } else {
      DSI::datashield.assign.expr(datasources, paste0("subsetGenoData_", 1, tfile), 
                                  paste0("subsetGenoDS(c('", genoData[[1]], "'), 'ethnic_snps')"))
      genoData <- paste0("subsetGenoData_", 1, tfile)
    }
  }

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
  
  partial_svd <- t(Reduce("cbind", res))
  total_svd <- svd(partial_svd)$v
  if(ncol(total_svd) < ncomp) {
    message('[ncomp] is too large, fixed to: ', ncol(total_svd))
    ncomp <- ncol(total_svd)
  }
  
  DSI::datashield.assign.expr(datasources, paste0("genoPCA_results", tfile), 
                              paste0("geno_pca_pooled_addPCDS(c(", paste(genoData, collapse = ", "), 
                                     "), c(", paste(unlist(total_svd[,1:ncomp]), collapse = ", "), 
                                     "), ", ncomp,")"))
  
  lapply(genoData_og, function(x){
    DSI::datashield.assign.expr(datasources, x, paste0("geno_pca_pooled_addPC2GenoDS(",
                                                       x, ", ", paste0("genoPCA_results", tfile), ")"))
  })
  
  
  return(list(pca_res = paste0("genoPCA_results", tfile), geno = genoData))

}
