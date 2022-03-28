#' Title
#'
#' @param genoData 
#' @param datasources 
#'
#' @return
#'
#' @examples
standardizeGenoData <- function(genoData, datasources){

  # TODO this is an internal function! no need to search datasources etc.

  # Get alternate allele frequencies for each SNP
  allel_freqs <- ds.alleleFrequency(genoData, type = "combined", fast = TRUE, datasources = datasources)

  # Calculate mean using alternate allele freq. Mean = 2 * MAF
  means <- allel_freqs$pooled_MAF * 2
  
  # Calculate standard error under HW equilibrium. SD_hw = sqrt(2 * MAF * (1 - MAF))
  sd_hw <- sqrt(means * (1 - allel_freqs$pooled_MAF))
  
  # Return data frame with rs ID means and sd_hw to be applied on other functions when reading
  # the actual geno data. The formula to standardize is (Xij - means) / sd_hw. Where Xij 
  # corresponds to the geno data in 0, 1, 2 coding.
  return(tibble(rs = allel_freqs$rs, means = means, sd_hw = sd_hw))
  
}