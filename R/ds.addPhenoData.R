#' @title Add Phenotype data to ExpressionSet
#' 
#' @description Add new phenotype data contained on a data.frame to an ExpressionSet. The ExpressionSet
#' may or may not already have phenotype data. If the data.frame contains a phenotype already present 
#' on the ExpressionSet, the server function will throw an exception. The pehnotypes data.frame has to 
#' contain an ID column which does not need to contain exactly the same individuals as the ExpressionSet, 
#' only the matching individuals will be updated, no new individuals will be introduced or removed from the 
#' ExpressionSet
#'
#' @param x \code{character} Name of the ExpressionSet on the study server
#' @param pheno \code{character} Name of the data.frame with the new phenotypes on the 
#' study server
#' @param identifier \code{character} (default \code{"ID"}) Name of the ID column on the phenotypes data.frame
#' @param alternate_eset_id \code{character} (default \code{NULL}) Alternate ID of the eSet pheno data, by default the rownames 
#' of the eSet pheno data act as ID, use this argument if the ID to merge the individuals is on a column of the pheno data. 
#' Input NULL for the standard behaviour of using the rownames of the pheno data as ID.
#' @param newobj.name \code{character} (default \code{NULL}) If \code{NULL}, the original ExpressionSet will be overwritten,
#' otherwise the new ExpressionSet will be assigned to a variable named after this argument
#' @param complete_cases \code{bool} (default \code{TRUE}) If \code{TRUE} only the matching individuals 
#' between the ExpressionSet and the phenotypes table will be included on the resulting ExpressionSet. If 
#' \code{FALSE} all the individuals on the input ExpressionSet will be on the output ExpressionSet
#' @param datasources a list of \code{\link{DSConnection-class}} (default \code{NULL}) objects obtained after login
#'
#' @return This function does not have an output. It creates (or overwrites) a data frame on the study server.
#' @export

ds.addPhenoData2eSet <- function(x, pheno, identifier = "ID", alternate_eset_id = NULL,
                            newobj.name = NULL, complete_cases = TRUE, datasources = NULL){
  
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }
  
  if(is.null(newobj.name)){
    newobj.name <- x
  }
  
  dsBaseClient:::isDefined(datasources, x)
  dsBaseClient:::isDefined(datasources, pheno)
  cls <- dsBaseClient:::checkClass(datasources, x)
  cls2 <- dsBaseClient:::checkClass(datasources, pheno)
  if(!any((cls %in% c("ExpressionSet")))){
    stop("'x' is not an 'ExpressionSet'")
  }
  if(!any((cls2 %in% c("data.frame")))){
    stop("The 'pheno' is not a 'data.frame'")
  }
  
  cally <- paste0("addPhenoDataDS(", x, ", ", pheno, ", '", identifier, "', ", 
                  if(is.null(alternate_eset_id)){"NULL, "}else{paste0("'",alternate_eset_id,"', ")},
                  complete_cases, ")")
  DSI::datashield.assign.expr(datasources, newobj.name, as.symbol(cally))
  
}
