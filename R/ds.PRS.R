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
#' When using a user provided prs_table table instead of a PGS catalog ID to calculate the PRS, it is important to note that 
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
#' Since the actual results of the PRS is sensitive information, the results are not returned to the client, 
#' however they can be merged into a table on the server(s). The main use of that is to add the PRS results 
#' to a phenotypes table and assess relationships between PRS scores and the phenotypes. This merge is performed 
#' via the individuals ID, specified on the argument (\code{table_id_column}); the table is specified using 
#' the argument \code{table}. When merging the results to a table, by default the column names will be: \cr
#' - When using \code{pgs_id}: \cr
#' + prs_\code{pgs_id} \cr
#' + prs_nw_\code{pgs_id} \cr
#' - When using \code{prs_table}: \cr
#' + prs_prs_custom_results \cr
#' + prs_nw_prs_custom_results \cr
#' 
#' If another designation is desired, make use of the \code{table_prs_name} argument, which by default is 
#' \code{NULL}. Note that this parameter only changes the tail 
#' of the names, the columns added (2) will begin by \code{prs_} and \code{prs_nw_}. This columns correspond to the 
#' actual PRS calculated and the PRS without weights (or PRS where all weights equal 1).
#'
#' @param resources \code{list} of all the VCF resources with biallelic genotype information. It is advised to 
#' have one VCF resource per chromosome, a big VCF file with all the information is always slower 
#' to use.
#' @param pgs_id \code{character} (default \code{NULL}) ID of the PGS catalog to be used to calculate the polygenic risk score. 
#' Polygenic Score ID & Name from https://www.pgscatalog.org/browse/scores/
#' @param prs_table \code{data.frame} (default \code{NULL}) Table containing the genomic region of interest to calculate 
#' the polygenic risk score. See the details for information of the structure this table has to have.\
#' @param table \code{character} (default \code{NULL}) If not \code{NULL}, it is the name of the table (on the server(s)) 
#' that will be used to merge the PRS results (typically a phenotypes table).
#' @param table_id_column \code{character} (default \code{NULL}) Argument only used when the \code{table} argument 
#' is supplied, it corresponds to the column name of the \code{table} that contains the individual IDs to perform 
#' the merge.
#' @param table_prs_name \code{character} (default \code{NULL}) If not \code{NULL} it's the name that will be 
#' used to design the column names added to \code{table}. Read the details for further information.
#' @param snp_threshold \code{numeric} (default \code{90}) Threshold to drop individuals. See details for 
#' further information.
#' @param datasources a list of \code{\link{DSConnection-class}} (default \code{NULL}) objects obtained after login
#'
#' @return The function has no client return. The results are stored on the server(s)
#' @export
#'

ds.PRS <- function(resources, pgs_id = NULL, prs_table = NULL, table = NULL, table_id_column = NULL,
                   table_prs_name = NULL, snp_threshold = 90, snp_assoc = FALSE, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  # Logic to check if prs_table or pgs_id is supplied. If both are supplied, the prs_table wil be used.
  if(!is.null(prs_table) & !is.null(pgs_id)){
    warning('Both [prs_table] and [pgs_id] supplied, only the [prs_table] will be used to calculate the PRS.')
    pgs_id <- NULL
  }
  else if(!is.null(prs_table)){
    # Check the structure of the supplied prs_table is as required (check details of roxygen), plus remove
    # unneeded columns
    if(all(c("chr_name", "chr_position", "effect_allele", "reference_allele", "effect_weight") %in% colnames(prs_table))){
      prs_table <- prs_table %>% dplyr::select(chr_name, chr_position, effect_allele, reference_allele, effect_weight)
      prs_table <- .recodeprs_table(prs_table)
    } else if (all(c("rsID", "effect_allele", "reference_allele", "effect_weight") %in% colnames(prs_table))){
      prs_table <- prs_table %>% dplyr::select(rsID, effect_allele, reference_allele, effect_weight)
      prs_table <- .recodeprs_table(prs_table)
    } else {
      stop('The supplied [prs_table] table is not structure as required. Please read the @details of ?dsOmicsClient::ds.PRS')
    }
  }

  else if(is.null(prs_table) & is.null(pgs_id)){
    error('Supply a [prs_table] or [pgs_id] to calculate the PRS.')
  }

  # Different cally builds depending on needing to send to the server the prs_table with ID or positions
  if(!is.null(prs_table)){
    if(all(c("chr_name", "start") %in% colnames(prs_table))){
      prs_table_type <- "chr_name"
    } else if ("rsID" %in% colnames(prs_table)){
      prs_table_type <- "rsID"
    }
  }

  # Build cally
  if(!is.null(pgs_id)){
    cally <- paste0("PRSDS(resources = c(", paste0(resources, collapse = ", "),
                    "), snp_threshold = ", snp_threshold, ", pgs_id = '", pgs_id, 
                    "', snp_assoc = ", snp_assoc, ")")
  } else {
    cally <- paste0("PRSDS(c(",
                    paste0(resources, collapse = ", "),
                    "), ", snp_threshold, ", pgs_id = NULL", ", snp_assoc = ", snp_assoc, ", '", paste(unlist(prs_table), collapse = "', '"),
                    "', ", if(prs_table_type == "rsID"){3}else{5}, ")")
  }
  
  DSI::datashield.assign.expr(datasources, "prs_results", cally)
  
  # If table is NULL the results will not be added (join by ID) to a table on the server
  if(!is.null(table)){
    if(is.null(table_id_column)){
      stop('Please provide the argument "table_id_column" to identify the column that contains the
           IDs to merge.')
    }
    if(is.null(table_prs_name)){
      if(!is.null(pgs_id)){
        warning("Default PRS custom name used ['", pgs_id, "']")
        table_prs_name <- pgs_id
      }else{
        warning("Default PRS custom name used ['prs_custom_results']")
        table_prs_name <- 'prs_custom_results'
      }
    }
    cally <- paste0("PRSDS_aux(prs_results, '", table_prs_name, "', ", table, ", '", table_id_column, "')")
    DSI::datashield.assign.expr(datasources, table, cally)
  }
}


.recodeprs_table <- function(scorings){
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
