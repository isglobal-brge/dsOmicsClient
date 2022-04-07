#' @title Principal Component Analysis (PCA) on SNP genotype data
#' 
#' @description PCA for genotype data on the study server
#' 
#' @details Pooled method implemented using block method ("Parallel Algorithms for the Singular Value Decomposition." Berry et al. 2005). 
#' The \code{snp_subset} option uses gene regions that have been linked to ethnic groupings, it is suggested to use this option to optimize 
#' the computing time and get noise-less principal components. 
#'
#' @param genoData \code{character} Name of the \code{GenotypeData} object on the server
#' @param snp_subset \code{bool} (default \code{TRUE}) Use only SNPs that have been linked to ethnic groups
#' @param standardize \code{bool} (default \code{TRUE}) Standardize the genotype data before performing the analysis
#' @param ncomp \code{integer} (default \code{5L}) Number of principal components to compute
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login
#'
#' @return Returns \code{list} with the names of the objects created on the server side that will be used if the principal 
#' components are to be plotted
#' @export

ds.PCA <- function(genoData, snp_subset = TRUE, standardize = TRUE, ncomp = 5L, datasources = NULL){
  
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
                    ,"), pca_rs, pca_means, pca_sd_hw", ")")
  } else {
    cally <- paste0("PCADS(c(", paste(genoData, collapse = ", ")
                    ,"), NULL, NULL, NULL", ")")
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
