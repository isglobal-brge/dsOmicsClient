#' @title Logistic regression analysis of pooled data for each SNP site in study
#' 
#' @description Fits a generalized linear model to gentic data for each SNP
#' in the data sets considered, using user specified outcome and phenotypic variabes as covariates
#' Outputs a matrix containing a beta value, standard error and p-value for each SNP
#'
#' @param snps.fit an optional parameter input as a character vector of
#' SNPs (rs numbers) that should be analysed. If missing all SNPs are analysed
#' @param model list of phenotypic variables to use as covariates in the regression analysis in the form:
#' "outcome ~ covar1 + covar2 +  ... + covarN"
#' @param genoData name of the DataSHIELD object to which the genotype (snpMatrix) and phenotypic data (data.frame) has been assigned
#' @param datasources Opal object or list of opal objects denoting the opal server(s) information
#' @param mc.cores optional parameter that allows the user to specify the number of CPU cores to use
#' @param type.p.adj
#' @param family
#' @param strata \code{character} Categorical variable to perform a stratified glm
#' 
#' @export
#' @examples
#'

ds.glmSNP <- function(snps.fit=NULL, model, genoData, datasources=NULL,
                      type.p.adj='fdr', mc.cores = 1, family = 'binomial', strata = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  mt <- as.formula(model)
  vars <- all.vars(mt)
  
# check vars are present in covars
  
  
  # Get snp rs's
  cally <- paste0("getVariable(", genoData, ", 'snp.rs.id')")
  snps <- unlist(DSI::datashield.aggregate(datasources,  cally))
  
  # Setting the SNPs to loop over if no SNPs are specified.
  # They are based on those that pass QC
  if(is.null(snps.fit)){
    snps.fit <- snps
  }
  if(is(snps.fit, "GenomicRanges")){
    cally <- paste0("getVariable(", genoData, ", 'snp.position')")
    positions <- unlist(DSI::datashield.aggregate(datasources,  cally))
    cally <- paste0("getVariable(", genoData, ", 'snp.chromosome')")
    chroms <- unlist(DSI::datashield.aggregate(datasources,  cally))
    snps_annoted <- data.frame(snps, chroms, positions)
    snps_annoted <- snps_annoted[snps_annoted$chrom %in% GenomeInfoDb::seqnames(snps.fit)@values,]
    snps.fit <- snps_annoted[between(snps_annoted$positions, IRanges::ranges(snps.fit)@start, 
                                     IRanges::ranges(snps.fit)@start + IRanges::ranges(snps.fit)@width),]$snps
    if(length(snps.fit) == 0){stop('No SNPs found on the range')}
    message('[', length(snps.fit), '] SNPs found on the range')
  }
  
  ans <- t(as.data.frame(parallel::mclapply(snps.fit, glmSNP, 
                                            snps=snps, 
                                            vars=vars,
                                            strata=strata,
                                            genoData=genoData, 
                                            datasources=datasources,
                                            family=family,
                                            mc.cores = mc.cores)))
  if (nrow(ans) > 1 & length(snps.fit) == 1) {
  #Ordering matrix by ascending p.value
  ans <- ans[order(ans[, 'p-value']), ]
  }
  #Calculating adjusted p-values
  ans <- cbind(ans, p.adj = p.adjust(ans[,'p-value'],  method=type.p.adj))
  
  class(ans) <- c("dsGlmSNP", class(ans))
  
  return(ans)
  
}
