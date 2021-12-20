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
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#'
#' @return A \code{data frame} with a row for each SNP. Columns "M" for males, "F" for females, and "all" 
#' for all scans give frequencies of the A allele. Sample size for males, females, and all is returned as 
#' "n.M", "n.F", and "n", respectively. "MAF" is the minor allele frequency over all scans.
#' @export

ds.alleleFrequency <- function(genoData, type = c("combined", "split"), datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if(length(genoData) > 1){
    res <- lapply(genoData, function(x){
      cally <- paste0("alleleFrequencyDS(", x, ")")
      return(DSI::datashield.aggregate(datasources, as.symbol(cally)))
    })
    
    res <- do.call(c, res)
    keys <- unique(names(res))
    alFreq <- NULL
    for(i in keys){
      ii <- which(names(res) == i)
      alFreq[[i]] <- do.call(rbind, res[ii])
    }
    
  } else {
    cally <- paste0("alleleFrequencyDS(", genoData, ")")
    alFreq <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  }
  
  if(type == "combined"){
    aux_df <- lapply(alFreq, function(x){
      return(data.frame(rs = x$rs, n = x$n, part_MAF = x$n * x$MAF))
    })
    
    # Check consistency of rs IDs
    if(!all(Vectorize(identical, "x")(lapply(aux_df, function(x){x$rs}), aux_df[[1]]$rs))){
      stop("The study servers have different rs ID, not consistent.")
    }
    
    # Get sum of n*MAF
    sum_part_MAF <- Reduce("+", lapply(aux_df, function(x){x$part_MAF}))
    
    # Get sum of n
    sum_n <- Reduce("+", lapply(aux_df, function(x){x$n}))
    
    # Assemble results
    results <- data.frame(rs = aux_df[[1]]$rs, n = sum_n, part_MAF = sum_part_MAF)
    
    # Get pooled MAF
    results$pooled_MAF <- results$part_MAF / results$n

    # Assemble final results to output
    alFreq <- as_tibble(results[,c(1,2,4)])
    
  } else if (!(type %in% c("combined", "split"))){
    stop("Allowed values for argument 'type' are ['combined' or 'split']")
  } 
  
  class(alFreq) <- c("dsalleleFrequency", class(alFreq))
  return(alFreq)
  
}
