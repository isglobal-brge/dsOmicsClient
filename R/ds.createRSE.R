#' @title Create a RangedSummarizedExperiment from a RNAseq count table and 
#' a phenotypes table
#' 
#' @details The RNAseq count table has to have the following structure: \cr
#' - First column: Entrez identificator \cr
#' - All the remaining columns: One col. per individual + annotation columns (no special
#' order required, they can even be mixed) \cr
#' The phenotypes table has to have the following structure: \cr
#' - First column: ID of the individuals. To match with the column names 
#' of the RNAseq table \cr
#' - All the remaining columns: One col. per phenotype
#'
#' @param rnaseq \code{character} Name of the RNAseq table on the server
#' @param phe \code{character} Name of the phenotypes table on the server
#' @param newobj.name \code{character} (default \code{NULL}) Name to be assigned 
#' on the server to the RangedSummarizedExperiment. If \code{NULL} the created RSE will 
#' be assigned to a variable named \code{'RSE'}
#' @param annotCols \code{character vector} (default \code{NULL}) Columns of the counts table to be used as 
#' annotation when creating the SummarizedExperiment. If \code{NULL}, biomaRt will be used to find 
#' annotation using the first column of the counts table (EntrezID)
#' @param datasources a list of \code{\link{DSConnection-class}} (default \code{NULL}) objects obtained after login
#'
#' @return This function does not have an output. It creates (or overwrites) a data frame on the study server.
#' @export

ds.createRSE <- function(rnaseq, phe, annotCols = NULL, newobj.name = NULL, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if(is.null(newobj.name)){
    newobj.name <- 'RSE'
  }
  
  dsBaseClient:::isDefined(datasources, rnaseq)
  dsBaseClient:::isDefined(datasources, phe)
  cls <- dsBaseClient:::checkClass(datasources, rnaseq)
  cls2 <- dsBaseClient:::checkClass(datasources, phe)
  if(!any((cls %in% c("data.frame")))){
    stop("[",rnaseq,"] is not a 'data.frame'")
  }
  if(!any((cls2 %in% c("data.frame")))){
    stop("[",phe,"] is not a 'data.frame'")
  }
  cally <- paste0("createRSEDS(", rnaseq, ", ", phe, ", ", 
                  if(is.null(annotCols)){"NULL"}else{paste0("'",paste(annotCols, collapse = "','"),"'")}, ")")
  DSI::datashield.assign.expr(datasources, newobj.name, cally)
  
}