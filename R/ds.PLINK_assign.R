#' @title Genome-wide association analysis (GWAS) using PLINK 
#' 
#' @description Performs GWAS using PLINK using shell command lines and stores the
#' results on the servers.
#' 
#' 
#' @param client ...
#' @param plink.command ...
#' @param newobj name of the object on the server-side which the resuls will be assigned to.
#' @param datasourcers ...
#' 
#' @return ...
#' 
#' @author Gonzalez JR, Wheater SM,
#'
#' @export 
#' 

ds.PLINK_assign <- function(client, plink.arguments, newobj = "plink_newobj", datasources=NULL){
  
    if (is.null(datasources)) {
        datasources <- DSI::datashield.connections_find()
    }
  
    args <- tolist(plink.arguments)
  
    cally <- paste0("plinkDS.assign(", client, ", ", "'", 
                    gsub(" ''=", " ", paste(names(args), args, collapse = "', '", sep = "'='")), "')")

    DSI::datashield.assign(datasources, newobj, as.symbol(cally))
}
