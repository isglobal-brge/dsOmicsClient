#' @title Get EGA URL
#'
#' @param id \code{character} ID of the file
#' @param chr \code{character} Chromosome to filter (encoding is server dependant)
#' @param start \code{numeric} Start position to filter
#' @param end \code{numeric} End position to filter
#'
#' @return
#' URL of ga4gh VCF/BAM file
#' @export

get_EGA_url <- function(id, chr, start, end){
  
  endpoint <- paste0("https://ega.ebi.ac.uk:8052/elixir/tickets/tickets/", id, 
                     "?referenceName=", chr, "&start=", format(start, scientific = FALSE),
                     "&end=", format(end, scientific = FALSE))
  
  return(endpoint)
}
