#' @title Hardy-Weinberg Equilibrium testing
#' 
#' @description This function performs exact Hardy-Weinberg Equilibrium testing (using Fisher's Test) 
#' over a selection of SNPs. It also counts genotype, calculates allele frequencies, 
#' and calculates inbreeding coefficients.
#'
#' @details (from exactHEW documentation): HWE calculations are performed with the \link{HWExact} function in the \link{GWASExactHW} package.
#' For the X chromosome, only female samples will be used in all calculations (since males are excluded from HWE
#'  testing on this chromosome). The X chromosome may not be included in a block with SNPs from other chromosomes.
#'   If the SNP selection includes the X chromosome, the scan annotation of genoData should include a "sex" column.
#' 
#' Y and M and chromsome SNPs are not permitted in the SNP selection, since the HWE test is not valid for these SNPs.
#' 
#' If permute=TRUE, alleles will be randomly shuffled before the HWE calculations. Running permutation can yield 
#' the expected distribution of p-values and corresponding confidence intervals.
#'
#' @param genoData \code{character} Name of the \code{\link{GenotypeData}} object on the server
#' @param sexcol \code{character} (default \code{"sex"}) Name of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object
#' @param male \code{character} (default \code{"M"}) Name of the male indicator of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object. (Note that it is case sensitive so it's not the same \code{male} than \code{Male})
#' @param female \code{character} (default \code{"M"}) Name of the female indicator of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object. (Note that it is case sensitive so it's not the same \code{female} than \code{Female})
#' @param chromosome \code{character} Chromosome to study. \code{"all"} to study all available chromosomes. Use
#' \code{ds.getChromosomeNames(genoData)} to retrieve the name encodings of the chromosomes.
#' @param geno.counts \code{bool} (default \code{TRUE}) if \code{TRUE}, genotype counts are returned in the output data.frame
#' @param block.size \code{numeric} (default \code{5000}) number of SNPs to read in at once
#' @param permute \code{bool} (default \code{FALSE}) logical indicator for whether to permute alleles before calculations
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#'
#' @return A \code{data frame} with the following columns: \cr
#' - \code{snpID}: the snpIDs \cr
#' - \code{chr}: chromosome SNPs are on \cr
#' 
#' If \code{geno.counts=TRUE}: \cr
#' - \code{nAA}: number of AA genotypes in samples \cr
#' - \code{nAB}: number of AB genotypes in samples \cr
#' - \code{nBB}: number of BB genotypes in samples \cr
#' - \code{MAF}: minor allele frequency \cr
#' - \code{minor.allele}: which allele ("A" or "B") is the minor allele \cr
#' - \code{f}: the inbreeding coefficient \cr
#' - \code{pval}: exact Hardy-Weinberg Equilibrium (using Fisher's Test) p-value. 
#' pval will be NA for monomorphic SNPs (MAF=0). \cr
#' 
#' @export

ds.exactHWE <- function(genoData, sexcol = "sex", male = "M", female = "F", chromosome = "all", geno.counts = TRUE,
                        block.size = 5000, permute = FALSE, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  cally <- paste0("exactHWEDS(", genoData, ", '", sexcol, "', '", male, "', '", female,
                  "', geno.counts = ", if(!is.null(geno.counts)){geno.counts}else{"NULL"}, 
                  ", chromosome = '", chromosome, "', block.size = ", block.size, ", permute = ",
                  permute, ")"
                  )
  ans <- datashield.aggregate(datasources, as.symbol(cally))
  
  class(ans) <- c("dsexactHWE", class(ans))
  return(ans)
  
}