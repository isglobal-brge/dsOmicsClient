#'
#' @title Server-side Differential Expression Analysis using edgeR
#' @description  This function performs a non-disclosive
#' Differential Expression Analysis using \code{edgeR} package functions
#' from Bioconductor.
#' 
#' @details Differential Expression Analysis of RNA-seq data based on  
#' genewise negative binomial generalized linear models using either 
#' likelihood ratio or quasi-likelihood F-tests. 
#' 
#' The steps implemented by DataSHIELD \code{ds.edgeR}  
#' client-side and \code{edgeRDS} server-side function is the following:\cr
#' (1) Create \code{DGEList} object\cr
#' (2) Filter genes using Expression Level. Implemented by \code{filterByExpr}
#' function from \code{edgeR} package.\cr
#' (3) Calculate the Normalization Factors using \code{calcNormFactors} function.
#' The \code{normalization} parameter can be set as follows:
#'     \itemize{
#'    \item{\code{"TMM"}}{: uses a trimmed mean of M-values between each pair of samples.} 
#'    \item{\code{"TMMwsp"}}{: TMM with singleton pairing.
#'    A variant of \code{TMM} that performs better  with a high proportion of zeros data} 
#'    \item{\code{"RLE"}}{: relative log expression. 
#'    The median library is calculated from the geometric mean of all columns 
#'    and the median ratio of each sample to the median library is taken as the scale factor.} 
#'    \item{\code{"upperquartile"}}{: the scale factors are calculated from the 75\% quantile of the 
#'    counts for each library, after removing genes that are zero in all libraries.} 
#'    \item{\code{"none"}}{: the normalization factors are set to 1.} 
#'     }
#' Note that normalization is only necessary for sample-specific effects.\cr
#' (4) Estimate the dispersion. 
#' The \code{dispersion} parameter can be set as follows:\cr
#'    \itemize{
#'    \item{\code{"both"}}{: estimate common and tagwise dispersions in one run (\code{estimateDisp}).} 
#'    \item{\code{"common"}}{: estimate common dispersion (\code{estimateCommonDisp}).}
#'    \item{\code{"tagwise"}}{: estimate tagwise dispersions (\code{estimateTagwiseDisp}). 
#'     Note that before needs to 
#'     estimate common dispersion. This is done automatically  when \code{tagwise} option
#'     is selected.}
#'     }
#' (5) Differential Expression Analysis. \code{test} can be set as follows:
#'     \itemize{
#'       \item{\code{"QLF"}}{: perform quasi-likelihood F-tests (\code{glmQLFit} and \code{glmQLFTest}).
#'       Highly recommended for differential expression analysis of RNA-seq.} 
#'       \item{\code{"LRT"}}{: perform likelihood ratio tests (\code{glmFit} and \code{glmLRT}).
#'       Useful in single cell RNA-seq and datasets with no replicates.}
#'     }
#' 
#' 
#' @param model formula indicating the condition and other covariates to be adjusted.
#' @param set \code{SummarizedExperiment} object.
#' @param test a character string specifying the test to carry out. 
#' \code{test} parameter can be set as \code{"QLF"} or \code{LRT}. Default \code{"QLF"}. 
#' For more information see \strong{Details}. 
#' @param dispersion a character string specifying the type of dispersion to estimate. 
#' This can be set as \code{both}, \code{common} or \code{tagwise}. 
#' For more information see \strong{Details}. Default \code{"both"}. 
#' @param normalization a character string specifying the normalization method to be used. 
#' This can be set as \code{"TMM"},\code{"TMMwsp"},\code{"RLE"},\code{"upperquartile"} or 
#' \code{"none"}. Default \code{"TMM"}. For more information see \strong{Details}.
#' @param levels a character or factor vector specifying the names of the parameters to be
#' used in \code{contrast}. Also, the design matrix can be specified. Default \code{"design"}.
#' @param coef an integer or character specifying the coefficients of the linear model
#' to be tested equal to zero. Default \code{2}. 
#' @param contrast the comparison to extract from the object to build the results table.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#' If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link{datashield.connections_default}}.
#' @return \code{ds.edgeR} returns to the client-side a data frame containing differential expression
#' results for the top genes sorted by adjusted p-value. 
#' @examples  
#' \dontrun{
#'  #required packages
#' 
#'  library(DSI)
#'  library(DSOpal)
#'  library(dsBaseClient)
#'  library(dsOmicsClient)
#' 
#'  # Connecting to the Opal servers
#'  builder <- DSI::newDSLoginBuilder()
#'  builder$append(server = "study1", url = "https://opal-demo.obiba.org/", 
#'                 user = "dsuser", password = "password", 
#'                 resource = "RSRC.tcga_liver", driver = "OpalDriver")
#'
#'  logindata <- builder$build()
#'
#'  conns <- datashield.login(logins = logindata, 
#'                            assign = TRUE, 
#'                            symbol = "res")
#'                            
#'  #coerce the resource to a RangedSummarizedExperiment
#'
#'  datashield.assign.expr(conns = conns, 
#'                         symbol = "rse",
#'                         expr = quote(as.resource.object(res)))
#'                         
#'  #Differential Expression analysis
#'  
#'  ##Default settings
#'  ds.edgeR(model =~ gdc_cases.demographic.gender,set = "rse")
#' 
#'  # Clear the Datashield R sessions and logout  
#'   
#'   DSI::datashield.logout(conns) 
#'   }
#' 
#' @author L. Abarrategui, for DataSHIELD development team
#' @export
#'

ds.edgeR <- function(model, set, test = "QLF", dispersion = "both",
                     normalization = "TMM", 
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
  
  if(normalization != "TMM" & normalization != "TMMwsp" & normalization != "RLE" & normalization != "upperquartile" & normalization != "none")
  {
    stop("Function argument 'normalization' has to be either 'TMM','TMMwsp','RLE','upperquartile' or 'none'", call.=FALSE)
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
  calltext <- call("edgeRDS",set, variable_names, intercept, dispersion, normalization, contrast, levels, test, coef)
  output <- datashield.aggregate(datasources, calltext)
  return(output)
}
#ds.edgeR
