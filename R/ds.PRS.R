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
#' When using a user provided ROI table instead of a PGS catalog ID to calculate the PRS, it is important to note that 
#' the provided data.frame has to have a very strict structure regarding column names (order is not relevant). Please 
#' follow one of this two schemas: \cr
#' - Schema 1 (provide SNP positions): \cr
#' + "chr_name", "chr_position", "effect_allele", "reference_allele", "effect_weight" \cr
#' \cr
#' - Schema 2 (provide SNP id's): \cr
#' + "rsID", "effect_allele", "reference_allele", "effect_weight" \cr
#' \cr
#' It is important to note that this "effect_weight" corresponds to the beta value of the SNP (log(OR)).
#' 
#' As a rule of thumb, it is advised to use when possible the Schema 1 (provide SNP positions), as the implementation 
#' to subset the VCF files is miles faster. 
#'
#' @param resources \code{list} of all the VCF resources with biallelic genotype information. It is advised to 
#' have one VCF resource per chromosome, a big VCF file with all the information is always slower 
#' to use.
#' @param pgs_id \code{character} (default \code{NULL}) ID of the PGS catalog to be used to calculate the polygenic risk score. 
#' Polygenic Score ID & Name from https://www.pgscatalog.org/browse/scores/
#' @param ROI \code{data.frame} (default \code{NULL}) Table containing the genomic region of interest to calculate 
#' the polygenic risk score. See the details for information of the structure this table has to have.
#' @param snp_threshold \code{numeric} (default \code{80}) Threshold to drop individuals. See details for 
#' further information.
#' @param datasources a list of \code{\link{DSConnection-class}} (default \code{NULL}) objects obtained after login
#'
#' @return \code{data.frame} were the rownames are the individuals. The columns found are: \cr
#' prs: Polygenic risk score per individual \cr
#' prs_nw: Polygenic risk score without weights (weight 1 for each risk allele) \cr
#' p_prs_nw: Risk probability using prs_nw and SNPassoc::pscore \cr
#' n_snps: Number of SNPs with information for each individual
#' @export
#'

ds.PRS <- function(resources, pgs_id = NULL, ROI = NULL, snp_threshold = 80, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  # Logic to check if ROI or pgs_id is supplied. If both are supplied, the ROI wil be used.
  if(!is.null(ROI) & !is.null(pgs_id)){
    warning('Both [ROI] and [pgs_id] supplied, only the [ROI] will be used to calculate the PRS.')
    pgs_id <- NULL
  }
  else if(!is.null(ROI)){
    # Check the structure of the supplied ROI is as required (check details of roxygen), plus remove
    # unneeded columns
    if(all(c("chr_name", "chr_position", "effect_allele", "reference_allele", "effect_weight") %in% colnames(ROI))){
      ROI <- ROI %>% dplyr::select(chr_name, chr_position, effect_allele, reference_allele, effect_weight)
      ROI <- .recodeROI(ROI)
    } else if (all(c("rsID", "effect_allele", "reference_allele", "effect_weight") %in% colnames(ROI))){
      ROI <- ROI %>% dplyr::select(rsID, effect_allele, reference_allele, effect_weight)
      ROI <- .recodeROI(ROI)
    } else {
      stop('The supplied [ROI] table is not structure as required. Please read the @details of ?dsOmicsClient::ds.PRS')
    }
  }
  else if(is.null(ROI) & !is.null(pgs_id)){
    # Get PGS data to later calculate genetic risk score
    ROI <- .retrievePGS(pgs_id)
  }
  else if(is.null(ROI) & is.null(pgs_id)){
    error('Supply a [ROI] or [pgs_id] to calculate the PRS.')
  }
  # Assign all the resources to GDSs to be studied later
  assigned_resources <- lapply(names(resources), function(x){
    if(colnames(ROI) == "rsID"){
      tryCatch({
        DSI::datashield.assign.expr(conns = datasources, symbol = paste0(x, "_gds"), 
                                    expr = as.symbol(paste0(
                                      "as.resource.object(x = ", resources[[x]],
                                      ", snps = c('", paste(ROI$rsID, collapse = "', '"), "'))"
                                    )))
        paste0(x, "_gds")
      }, error = function(w){
        NULL
      })
      
    } else {
      tryCatch({
        DSI::datashield.assign.expr(conns = datasources, symbol = paste0(x, "_gds"), 
                                    expr = as.symbol(paste0(
                                      "as.resource.object(x = ", resources[[x]], ", seq = c('", 
                                      paste(ROI$chr_name, collapse = "', '"), "'), start.loc = c(",
                                      paste(ROI$start, collapse = ", "), "), end.loc = c(",
                                      paste(ROI$end, collapse = ", "), "))"
                                    )))
        paste0(x, "_gds")
      }, error = function(w){
        NULL
      })
    }
  })
  # Different cally builds depending on needing to send to the server the ROI with ID or positions
  if(all(c("chr_name", "chr_position") %in% colnames(ROI))){
    ROI_type <- "chr_name"
  } else if ("rsID" %in% colnames(ROI)){
    ROI_type <- "rsID"
  }
  browser()
  cally <- paste0("PRSDS(c(",
                  paste0(assigned_resources, collapse = ", "),
                  "), NULL, ", snp_threshold, ", '", paste(unlist(ROI), collapse = "', '"),
                  "', ", if(ROI_type == "rsID"){3}else{5}, ")")
  
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
  } else {
    data <- data.frame(rsID = scorings$rsID,
                       # reference_allele = scorings$reference_allele,
                       effect_allele = scorings$effect_allele,
                       effect_weight = scorings$effect_weight)
    if("weight_type" %in% colnames(scorings)){
      data <- data %>% tibble::add_column(weight_type = scorings$weight_type)}
  }
  # If weight_type is present and is equal to "OR" or "HR", convert the effect_weight
  # to log(effect_weight) to get the beta
  # This is done row by row in case not all rows have the same weight_type
  if(!is.null(data$weight_type)){
    for(i in seq(1, nrow(data))){
      if(c("OR", "HR") %in% data$weight_type[i]){
        data$effect_weight[i] <- log(data$effect_weight[i])
      }
    }
    data$weight_type <- NULL
  }
  
  return(data)
  
}



.recodeROI <- function(scorings){
  if(c("chr_name", "chr_position") %in% colnames(scorings)){
    data <- data.frame(chr_name = scorings$chr_name,
                       start = scorings$chr_position,
                       end = scorings$chr_position,
                       # reference_allele = scorings$reference_allele,
                       effect_allele = scorings$effect_allele,
                       effect_weight = scorings$effect_weight)
    return(data)
  } else {
    data <- data.frame(rsID = scorings$rsID,
                       # reference_allele = scorings$reference_allele,
                       effect_allele = scorings$effect_allele,
                       effect_weight = scorings$effect_weight)
    return(data)
  }
}