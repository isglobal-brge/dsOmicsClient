#' Title
#'
#' @param resources 
#' @param pgs_id 
#' @param datasources 
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param pgs_id 
#'
#' @return
#'
#' @examples
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
