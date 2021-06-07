#' @title Get Ploygenic Risk Score
#' @subtitle Using server side VCF resources and retrieving the risk 
#' factors from the PGS catalog (https://www.pgscatalog.org/)
#' 
#' @details This function resolves a list of resources subsetting them by the 
#' SNPs of risk, this does not ensure that all the SNPs of risk will be found on the 
#' data. From all the found SNPs of risk, if an individual has less than 'snp_threshold' (percetage)
#' of SNPs with data, it will be dropped (SNP with no data is marked on the VCF as ./.). If an individual 
#' passes this threshold filter but still has SNPs with no data, those SNPs will be counted on the 
#' polygenic risk score as non-risk-alleles, to take this infomation into account, the number of SNPs 
#' with data for each individual is returned as 'n_snps'.
#'
#' @param resources \code{list} of all the VCF resources with biallelic genotype information. It is advised to 
#' have one VCF resource per chromosome, a big VCF file with all the information is always slower 
#' to use.
#' @param pgs_id \code{character} ID of the PGS catalog to be used to calculate the polygenic risk score. 
#' Polygenic Score ID & Name from https://www.pgscatalog.org/browse/scores/
#' @param snp_threshold \code{numeric} (default \code{80}) Threshold to drop individuals. See details for 
#' further information.
#' @param datasources a list of \code{\link{DSConnection-class}} (default \code{NULL}) objects obtained after login
#'
#' @return \code{data.frame} were the rownames are the individuals. The columns found are: \cr
#' prs: Polygenic risk score per individual \cr
#' prs_nw: Polygenic risk score without weights (weight 1 for each risk allele) \cr
#' n_snps: Number of SNPs with information for each individual
#' @export
#'

ds.PRS <- function(resources, pgs_id, snp_threshold = 80, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  # Get PGS data to later calculate genetic risk score
  ROI <- .retrievePGS(pgs_id)
  # Assign all the resources to GDSs to be studied later
  lapply(names(resources), function(x){
    if(colnames(ROI) == "rsID"){
      DSI::datashield.assign.expr(conns = datasources, symbol = paste0(x, "_gds"), 
                                  expr = as.symbol(paste0(
                                    "as.resource.object(x = ", resources[[x]],
                                    ", snps = c('", paste(ROI$rsID, collapse = "', '"), "'))"
                                  )))
    } else {
      DSI::datashield.assign.expr(conns = datasources, symbol = paste0(x, "_gds"), 
                                  expr = as.symbol(paste0(
                                    "as.resource.object(x = ", resources[[x]], ", seq = c('", 
                                    paste(ROI$chr_name, collapse = "', '"), "'), start.loc = c(",
                                    paste(ROI$start, collapse = ", "), "), end.loc = c(",
                                    paste(ROI$end, collapse = ", "), "))"
                                  )))
    }
  })
  cally <- paste0("PRSDS(c(",
                  paste0(paste0(names(resources), "_gds"), collapse = ", "),
                  "), '", pgs_id, "', ", snp_threshold, ")")
  DSI::datashield.aggregate(datasources, as.symbol(cally))
}

#' @title Internal function: Get PGS catalog table of polygenic risks
#'
#' @param pgs_id \code{character} ID of the PGS catalog to be used to calculate the polygenic risk score. 
#' Polygenic Score ID & Name from https://www.pgscatalog.org/browse/scores/
#'
#' @return \code{data.frame} with the columns: \cr
#' - If chr_name and poition are found: \cr
#' + start \cr
#' + end \cr
#' + effect_allele \cr
#' + effect_weight \cr
#' + weight_type (if present on the catalog) \cr
#' - If rsID is found: \cr
#' + rsID \cr
#' + effect_allele \cr
#' + effect_weight \cr
#' + weight_type (if present on the catalog) \cr
#'

.retrievePGS <- function(pgs_id){
  pgs <- httr::GET(paste0("https://www.pgscatalog.org/rest/score/", pgs_id))
  pgs_text <- httr::content(pgs, "text")
  pgs_scoring_file <- jsonlite::fromJSON(pgs_text, flatten = TRUE)$ftp_scoring_file
  if(is.null(pgs_scoring_file)){stop('[', pgs_id, '] Not found in www.pgscatalog.org')}
  
  # Download scoring file
  destination <- tempfile(fileext = ".txt.gz")
  download.file(pgs_scoring_file, destination)
  
  # Unzip file
  unziped_file <- tempfile(fileext = ".txt")
  R.utils::gunzip(destination, unziped_file)
  
  # Read file into R
  # Logic https://www.pgscatalog.org/downloads/
  scorings <- read.delim(unziped_file, comment.char = "#")
  if(c("chr_name", "chr_position") %in% colnames(scorings)){
    data <- data.frame(chr_name = scorings$chr_name,
                       start = scorings$chr_position,
                       end = scorings$chr_position,
                       # reference_allele = scorings$reference_allele,
                       effect_allele = scorings$effect_allele,
                       effect_weight = scorings$effect_weight)
    if("weight_type" %in% colnames(scorings)){
      data <- data %>% tibble::add_column(weight_type = scorings$weight_type)}
    return(data)
  } else {
    data <- data.frame(rsID = scorings$rsID,
                       # reference_allele = scorings$reference_allele,
                       effect_allele = scorings$effect_allele,
                       effect_weight = scorings$effect_weight)
    if("weight_type" %in% colnames(scorings)){
      data <- data %>% tibble::add_column(weight_type = scorings$weight_type)}
    return(data)
  }
}
