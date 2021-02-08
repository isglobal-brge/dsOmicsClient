#' @title Interface to run SNPtest commands on a ssh connection (with SNPtest 2.5.2 installed)
#'
#' @param client \code{character} Name of the ssh resource on the server
#' @param snptest.arguments \code{character} SNPtest arguments, the leading "snptest" and trailing 
#' "-o output.out" do not need to be provided, they are handled on the server side.
#' @param datasources a list of \code{\link{DSConnection-class}} (default \code{NULL}) objects obtained after login
#'
#' @return List containing: \cr
#' - Results: Table of results (typical SNPtest output file) \cr
#' - Output: Console output of the ssh query
#' 
#' @export
#' 

ds.snptest <- function(client, snptest.arguments, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  # Remove double spaces
  snptest.arguments <- gsub("^ *|(?<= ) | *$", "", snptest.arguments, perl = TRUE)
  
  args <- strsplit(snptest.arguments, " ")[[1]]
  
  cally <- paste0("snptestDS(", client, ", '", paste0(args, collapse = "', '"), "')")
 
  ans <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  class(ans) <- c("dssnptest", class(ans))
  
  return(ans)
  
}
