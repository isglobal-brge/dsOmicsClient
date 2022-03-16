#' @title Allelic frequency
#'
#' @description Calculates the frequency of the A allele on a server side GenotypeData object.
#' 
#' @details (From alleleFrequency documentation): Counts male heterozygotes on the X and Y chromosomes
#'  as missing values, and any genotype for females on the Y chromosome as missing values. A "sex" 
#'  variable must be present in the scan annotation slot of genoData.
#'  Samples with missing sex are included in the allele counts for "all" and "MAF" for autosomes, 
#'  but not for sex chromosomes.
#'
#' @param genoData \code{character} Name of the \code{\link{GenotypeData}} object on the server
#' @param type \code{character} (default \code{"combined"}) Type of analysis, if (\code{"split"}), the  
#' frequencies will be calculated for each study server. If  (\code{"combined"}) 
#' a pooled methodology will be performed.
#' @param method \code{character} (default \code{"fast"}) If \code{"fast"} an optimized fast method will 
#' be used; if \code{"slow"} the function \code{GWASTools::alleleFrequency} will be used, this is notably 
#' slower.
#' @param snpBlock \code{integer} (default \code{5000L}) number of SNPs to read at each iteration, 
#' tune this parameter to improve performance. This argument is only used when \code{method = "fast"}
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#'
#' @return A \code{data frame} with a row for each SNP. Columns "M" for males, "F" for females, and "all" 
#' for all scans give frequencies of the A allele. Sample size for males, females, and all is returned as 
#' "n.M", "n.F", and "n", respectively. "MAF" is the minor allele frequency over all scans.
#' @export

ds.alleleFrequency <- function(genoData, type = "combined", method = "fast",
                               snpBlock = 5000L, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if(length(genoData) > 1){
    res <- lapply(genoData, function(x){
      if(method == "fast"){
        cally <- paste0("fastAlleleFrequencyDS(", x, ", ", snpBlock, ")")
        return(DSI::datashield.aggregate(datasources, as.symbol(cally)))
      } else {
        cally <- paste0("alleleFrequencyDS(", x, ")")
        return(DSI::datashield.aggregate(datasources, as.symbol(cally)))
      }
    })
    
    res <- do.call(c, res)
    keys <- unique(names(res))
    alFreq <- NULL
    for(i in keys){
      ii <- which(names(res) == i)
      alFreq[[i]] <- do.call(rbind, res[ii])
    }
    
  } else {
    if(method == "fast"){
      cally <- paste0("fastAlleleFrequencyDS(", genoData, ", ", snpBlock, ")")
      alFreq <- DSI::datashield.aggregate(datasources, as.symbol(cally))
    } else {
      cally <- paste0("alleleFrequencyDS(", genoData, ")")
      alFreq <- DSI::datashield.aggregate(datasources, as.symbol(cally))
    }
  }
  
  if(type == "combined"){
    aux_df <- lapply(alFreq, function(x){
      return(data.frame(rs = x$rs, n = x$n, part_MAF = x$n * x$MAF))
    })
    # Check consistency of rs IDs
    common_rs <- Reduce(intersect, lapply(aux_df, function(x){x$rs}))
    if(length(common_rs) == 0){
      stop("No common SNPs between the study servers.")
    }
    # Use only common SNPs for each server
    aux_df <- lapply(aux_df, function(x){
      x[x$rs %in% common_rs,]
    })

    # Get sum of n*MAF
    sum_part_MAF <- Reduce("+", lapply(aux_df, function(x){x$part_MAF}))
    
    # Get sum of n
    sum_n <- Reduce("+", lapply(aux_df, function(x){x$n}))
    
    # Assemble results
    results <- data.frame(rs = aux_df[[1]]$rs, n = sum_n, part_MAF = sum_part_MAF)
    
    # Get pooled MAF
    results$pooled_MAF <- results$part_MAF / results$n

    # Assemble final results to output
    alFreq <- tibble::as_tibble(results[,c(1,2,4)])
    
  } else if (!(type %in% c("combined", "split"))){
    stop("Allowed values for argument 'type' are ['combined' or 'split']")
  } 
  
  class(alFreq) <- c("dsalleleFrequency", class(alFreq))
  return(alFreq)
  
}
