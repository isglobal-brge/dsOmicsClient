% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ds.limma.R
\name{ds.limma}
\alias{ds.limma}
\title{Server-side Differential Gene Expression analysis using limma}
\usage{
ds.limma(
  model,
  Set,
  type.data = "microarray",
  contrasts = NULL,
  levels = "design",
  coef = 2,
  sva = FALSE,
  annotCols = NULL,
  method = "ls",
  robust = FALSE,
  normalization = "none",
  voomQualityWeights = FALSE,
  big = FALSE,
  sort.by = "none",
  datasources = NULL
)
}
\arguments{
\item{model}{formula indicating the condition (left side) and other covariates to be adjusted for
(i.e. condition ~ covar1 + ... + covar2). The fitted model is: feature ~ condition + covar1 + ... + covarN}

\item{Set}{name of the DataSHIELD object to which the ExpresionSet or RangedSummarizedExperiment has been assigned}

\item{type.data}{optional parameter that allows the user to specify the number of CPU cores to use during}

\item{sva}{logical value}

\item{annotCols}{the column names of the annotation available in the ExpresionSet or RangedSummarizedExperiment (see fData() function)}

\item{method}{String indicating the method used in the regression: "ls" or
"robust". (Default: "ls")}

\item{robust}{Logical value indicating whether robust method is applied in the eBayes function of limma. Default is FALSE.}

\item{normalization}{String indicating the normalize method used when using voom for RNAseq data
(see normalized.method argument in limma::vomm for possible values)
#' @param voomQualityWeights Logical value indicating whether limma::voomWithQualityWeights should be used instead of
limma::voom. Default is FALSE and hence the pipeline uses limma::voom to transform RNAseq data.}

\item{datasources}{a list of \code{\link{DSConnection-class}} objects obtained after login.
If the \code{datasources} argument is not specified
the default set of connections will be used: see \code{\link{datashield.connections_default}}.}
}
\description{
This function performs a non-disclosive
Differential Gene Expression Analysis using \code{limma} package from Bioconductor.
}
\details{
Implementation of Bioconductor's \code{limma} in DataSHIELD using \code{MEAL} package
}
