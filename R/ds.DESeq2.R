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
#' @param set a \code{RangedSummarizedExperiment} object.  
#' @param test \code{"Wald"} or \code{"LRT"}.
#' @param fitType \code{"parametric"}, \code{"local"}, or \code{"mean"}.
#' @param sfType 	\code{"ratio"}, \code{"poscounts"}, or \code{"iterate"}.
#' @param reduced a reduced formula for \code{test="LRT"}
#' @param contrast the comparison to extract from the object to build the results table 
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#' If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link{datashield.connections_default}}.
#' @return Data frame containing the results of the Differential Gene Expression analysis 
#' @author L. Abarrategui, for DataSHIELD development team
#' @export
#'

ds.DESeq2 <- function(model, set, test = "Wald",
                    fitType = "parametric", sfType = "ratio",
                    reduced = NULL, 
                    contrast = NULL, datasources=NULL)
{   
    if (is.null(datasources))
    {
        datasources <-  DSI::datashield.connections_find()
    }
    
    if(is.null(model)){
      stop(" Please provide a valid model formula", call.=FALSE)
    }
    
    vars <- all.vars(model)
    vars<- paste(vars, collapse=",")
    
    if(test != "Wald" & test != "LRT")
      {
      stop("Function argument 'test' has to be either 'Wald' or  'LRT'", call.=FALSE)
    }   
    
    if(fitType != "parametric" & fitType != "local" & fitType != "mean")
    {
      stop("Function argument 'fitType' has to be either 'parametric','local' or 'mean'", call.=FALSE)
    }  
    
    if(sfType != "ratio" & sfType != "poscounts" & sfType != "iterate")
    {
      stop("Function argument 'sfType' has to be either 'ratio', 'poscounts', 'iterate'", call.=FALSE)
    }  
    
    if(!is.null(reduced))
    {
      reduced.term <- attributes(terms(as.formula(reduced)))$term.labels
      
      if (length(reduced.term)==0){
        reduced<-attributes(stats::terms(stats::as.formula(reduced)))$intercept
        reduced<-as.character(reduced)
        
      }else{
        reduced<-all.vars(reduced)
        reduced<-paste(reduced, collapse=",")
      }
    }
    
    if(!is.null(contrast))
      {
      contrast <- paste(contrast, collapse=",")
    }
    # call the server side function
    calltext <- call("DESeq2DS",vars, set, test, fitType, sfType, reduced, contrast)
    output <- datashield.aggregate(datasources, calltext)
    return(output)
}
#ds.DESeq2