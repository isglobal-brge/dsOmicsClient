#' @title Genome-wide association analysis (GWAS) using PLINK 
#' 
#' @description Performs GWAS using PLINK using shell command lines
#' 
#' 
#' @param client ...
#' @param plink.command ...
#' @param datasourcers ...
#' 
#' @return ...
#' 
#' @author Gonzalez, JR.
#'
#' @export 
#' 

ds.PLINK <- function(client, plink.arguments, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  args <- tolist(plink.arguments)
  
  cally <- paste0("plinkDS(", client, ", ", "'", 
                  gsub(" ''=", " ", paste(names(args), args, collapse = "', '", sep = "'='")), "')")

  ans <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  
  class(ans) <- c("dsPLINK", class(ans))
  
  return(ans)
  
}
