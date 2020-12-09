#' @title Get GA4GH URL
#'
#' @param url \code{character} URL of the server
#' @param id \code{character} ID of the file
#' @param chr \code{character} Chromosome to filter (encoding is server dependant)
#' @param start \code{numeric} Start position to filter
#' @param end \code{numeric} End position to filter
#' @param type \code{character} Type of file (VCF/BAM)
#'
#' @return
#' URL of ga4gh VCF/BAM file
#' @export

get_ga4gh_url <- function(url, id, chr, start, end, type){
  if(type == "VCF"){
    endpoint <- paste0(url, "/variants/", id, "?format=VCF&referenceName=", chr, "&start=", format(start, scientific = FALSE),
          "&end=", format(end, scientific = FALSE))
  }
  else if(type == "BAM"){
    endpoint <- paste0(url, "/reads/", id, "?format=BAM&referenceName=", chr, "&start=", format(start, scientific = FALSE),
           "&end=", format(end, scientific = FALSE))
  }
  return(endpoint)
}
