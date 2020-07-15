#'
#' @title Server-side Differential Gene Expression analysis using DESeq2
#' @description  This function performs a non-disclosive
#' Differential Gene Expression Analysis using \code{DESeq2} package
#' from Bioconductor.
#' 
#' @details Implementation of Bioconductor's \code{DESeq2} in DataSHIELD 
#' infraestructure. 
#' 
#' This function is similar to Bioconductor's 
#' \code{DESeq2} package \code{DESeq} function.
#' @param model formula indicating the condition and other covariates to be ajusted
#' @param set a 
#' @param test 
#' @param contrast the comparison to extract from the object to build the results table 
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#' If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link{datashield.connections_default}}.
#' @return Data frame containing the results of the Differential Gene Expression analysis 
#' @author L. Abarrategui, for DataSHIELD development team
#' @export
#'

ds.edgeR <- function(model, set, test = "QLF", dispersion = "both",
                     levels = "design", coef = 2,
                      contrast = NULL, datasources=NULL)
{   
  if (is.null(datasources))
  {
    datasources <-  DSI::datashield.connections_find()
  }
  
  if(is.null(model)){
    
    stop(" Please provide a valid model formula", call.=FALSE)
  }
  model.terms<-attributes(stats::terms(stats::as.formula(model)))
  variable_names <- model.terms$term.labels
  variable_names <- paste(variable_names, collapse=",")
  intercept <- model.terms$intercept
  
  
  test <- charmatch(test, c("QLF", "LRT"))
  if(is.na(test))
  {
    stop("Function argument 'test' has to be either 'QLF' or  'LRT'", call.=FALSE) 
  }
    
  
  
  dispersion <- charmatch(dispersion, c("both", "common","tagwise"))
  
  if(is.na(dispersion))
  {
    stop("Function argument 'fitType' has to be either 'both','common' or 'tagwise'", call.=FALSE)  
  }
    
  
  if (levels[1] != "design")
  {
    levels <- paste(levels, collapse=",")
  }
  
  if(!is.null(contrast))
  {
    contrast <- paste(contrast, collapse=",")
  }
  # call the server side function
  calltext <- call("edgeRDS",set, variable_names, intercept, dispersion, contrast, levels, test, coef)
  output <- datashield.aggregate(datasources, calltext)
  return(output)
}
#ds.edgeR