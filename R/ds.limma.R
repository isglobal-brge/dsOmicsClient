#'
##' @title Server-side Differential Gene Expression analysis using limma
##'
##' @description  This function performs a non-disclosive
##' Differential Gene Expression Analysis using \code{limma} package from Bioconductor.
#'
##' @details Implementation of Bioconductor's \code{limma} in DataSHIELD using \code{MEAL} package
#'
##' @param model formula indicating the condition (left side) and other covariates to be adjusted for
##' (i.e. condition ~ covar1 + ... + covar2). The fitted model is: feature ~ condition + covar1 + ... + covarN
##' @param Set name of the DataSHIELD object to which the ExpresionSet or RangedSummarizedExperiment has been assigned
##' @param type.data optional parameter that allows the user to specify the number of CPU cores to use during
##' @param sva logical value
##' @param annotCols the column names of the annotation available in the ExpresionSet or RangedSummarizedExperiment (see fData() function)
##' @param method String indicating the method used in the regression: "ls" or
#' "robust". (Default: "ls")
#' @param robust Logical value indicating whether robust method is applied in the eBayes function of limma. Default is FALSE.
#' @param normalization String indicating the normalize method used when using voom for RNAseq data
#' (see normalized.method argument in limma::vomm for possible values)
#' #' @param voomQualityWeights Logical value indicating whether limma::voomWithQualityWeights should be used instead of
#' limma::voom. Default is FALSE and hence the pipeline uses limma::voom to transform RNAseq data.
#'
##' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login.
#' If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link{datashield.connections_default}}.
##'
##' @export
##' @examples
##'

ds.limma <- function(model, Set, type.data="microarray",
                     contrasts = NULL, levels = "design", coef = 2,
                     sva=FALSE, annotCols=NULL, method = "ls", robust = FALSE,
                     normalization = "none", voomQualityWeights = FALSE,
                     big = FALSE, sort.by = "none", datasources=NULL){

  type <- charmatch(type.data, c("microarray", "RNAseq"))
  if(is.na(type))
    stop("type.data must be 'microarray' or 'RNAseq'")

  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  if(is.null(model)){
    stop(" Please provide a valid model formula", call.=FALSE)
  }
  method <- charmatch(method, c("ls", "robust"))
  if(is.null(method)){
    stop("method must be 'ls' or 'robust'")
  }

  mt <- all.vars(model)
  variable_names <- mt[1]
  if (length(mt)>1)
    covariable_names <- paste(mt[-1], collapse=",")
  else
    covariable_names <- NULL

  if (!is.null(annotCols)){
    annotCols <- paste(annotCols, collapse=",")
  }

  if (levels[1] != "design")
    {
    levels <- paste(levels, collapse=",")
  }
  if (!is.null(contrasts))
  {
    contrasts<-paste(contrasts, collapse=",")
  }
  
  calltext <- call("limmaDS", Set, variable_names,
                   covariable_names, type, contrasts,
                   levels, coef, sva, annotCols, method,
                   robust, normalization, voomQualityWeights, big, sort.by)
  datashield.assign.expr(datasources, "res", calltext)
  
  calltext <- call("limmaDS2", as.symbol(eval(Set)), quote(res), type, contrasts, 
                   coef, annotCols, robust, sort.by)
  # ans <- datashield.aggregate(datasources, calltext)
  ans <- datashield.aggregate(datasources, calltext)
  class(ans) <- c("dsLimma", class(ans))
  return(ans)
}
