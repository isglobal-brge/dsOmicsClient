% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ds.genoDimensions.R
\name{ds.genoDimensions}
\alias{ds.genoDimensions}
\title{Get main dimensions of Genotype data}
\usage{
ds.genoDimensions(x, datasources = NULL)
}
\arguments{
\item{x}{\code{character} Name of the \code{GenotypeData} or \code{GdsGenotypeReader} object on the server}

\item{datasources}{a list of \code{DSConnection-class} objects obtained after login}
}
\value{
\code{list} with results
}
\description{
Get the number of SNPs, number of scans and number of chromosomes on the genotype file
}
