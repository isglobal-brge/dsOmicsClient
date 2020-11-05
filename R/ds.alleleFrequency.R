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
#' @param sexcol \code{character} (default \code{"sex"}) Name of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object
#' @param male \code{character} (default \code{"M"}) Name of the male indicator of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object. (Note that it is case sensitive so it's not the same \code{male} than \code{Male})
#' @param female \code{character} (default \code{"M"}) Name of the female indicator of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object. (Note that it is case sensitive so it's not the same \code{female} than \code{Female})
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#'
#' @return A \code{data frame} with a row for each SNP. Columns "M" for males, "F" for females, and "all" 
#' for all scans give frequencies of the A allele. Sample size for males, females, and all is returned as 
#' "n.M", "n.F", and "n", respectively. "MAF" is the minor allele frequency over all scans.
#' @export

ds.alleleFrequency <- function(genoData, sexcol = "sex", male = "M", female = "F", datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0("alleleFrequencyDS(", genoData, ", '", sexcol, "', '", male, "', '", female, "')")
  alFreq <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  class(alFreq) <- c("dsalleleFrequency", class(alFreq))
  return(alFreq)
  
}