% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ds.exactHWE.R
\name{ds.exactHWE}
\alias{ds.exactHWE}
\title{Hardy-Weinberg Equilibrium testing}
\usage{
ds.exactHWE(
  genoData,
  chromosome = "all",
  geno.counts = TRUE,
  block.size = 5000,
  permute = FALSE,
  controls = FALSE,
  controls_column = NULL,
  datasources = NULL
)
}
\arguments{
\item{genoData}{\code{character} Name of the \code{\link{GenotypeData}} object on the server}

\item{chromosome}{\code{character} Chromosome to study. \code{"all"} to study all available chromosomes. Use
\code{ds.getChromosomeNames(genoData)} to retrieve the name encodings of the chromosomes.}

\item{geno.counts}{\code{bool} (default \code{TRUE}) if \code{TRUE}, genotype counts are returned in the output data.frame}

\item{block.size}{\code{numeric} (default \code{5000}) number of SNPs to read in at once}

\item{permute}{\code{bool} (default \code{FALSE}) logical indicator for whether to permute alleles before calculations}

\item{controls}{\code{bool} (default \code{FALSE}) logical to calculate the HWE test only on the controls}

\item{controls_column}{\code{character} (default \code{NULL}) name of the case/controls column of the covariates 
used to create the GenotypeData object. Only used if \code{controls = TRUE}}

\item{datasources}{a list of \code{\link{DSConnection-class}} objects obtained after login.}
}
\value{
A \code{data frame} with the following columns: \cr
- \code{snpID}: the snpIDs \cr
- \code{chr}: chromosome SNPs are on \cr

If \code{geno.counts=TRUE}: \cr
- \code{nAA}: number of AA genotypes in samples \cr
- \code{nAB}: number of AB genotypes in samples \cr
- \code{nBB}: number of BB genotypes in samples \cr
- \code{MAF}: minor allele frequency \cr
- \code{minor.allele}: which allele ("A" or "B") is the minor allele \cr
- \code{f}: the inbreeding coefficient \cr
- \code{pval}: exact Hardy-Weinberg Equilibrium (using Fisher's Test) p-value. 
pval will be NA for monomorphic SNPs (MAF=0). \cr
}
\description{
This function performs exact Hardy-Weinberg Equilibrium testing (using Fisher's Test) 
over a selection of SNPs. It also counts genotype, calculates allele frequencies, 
and calculates inbreeding coefficients.
}
\details{
(from exactHEW documentation): HWE calculations are performed with the \link{HWExact} function in the \link{GWASExactHW} package.
For the X chromosome, only female samples will be used in all calculations (since males are excluded from HWE
 testing on this chromosome). The X chromosome may not be included in a block with SNPs from other chromosomes.
  If the SNP selection includes the X chromosome, the scan annotation of genoData should include a "sex" column.

Y and M and chromsome SNPs are not permitted in the SNP selection, since the HWE test is not valid for these SNPs.

If permute=TRUE, alleles will be randomly shuffled before the HWE calculations. Running permutation can yield 
the expected distribution of p-values and corresponding confidence intervals.
}
