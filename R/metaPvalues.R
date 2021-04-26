#' @title Meta-analysis of p-values
#' 
#' @description Performs meta-analys of pvalues using the sum of logs method (Fisher's method)
#'
#' @param x a \code{dsOmics} object obtained from \code{ds.limma}, \code{ds.GWAS} or \code{ds.PLINK} functions applied o 2 or more studies
#' @param ... other arguments of \code{fitNullModel} function in GENESIS package
#' 
#' @return a matrix with features p-values of each study and its combination 
#' 
#' @author Gonzalez, JR.
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom tibble add_column
#' @import dplyr
#' @import metap 

metaPvalues <- function(x,  ...){
  
  if(length(x)==1)
    stop('Nothing to be meta-analyzed. There is only a single study')


  if (inherits(x, 'dsLimma')){
    ff <- function(x, y){
      inner_join(x, y, 
                 by='id')
    }
  } else if (inherits(x, 'dsGWAS')){
    ff <- function(x, y){
      inner_join(x, y, 
                 by='variant.id')
    }
  } else{
      stop("Object should be of class 'dsLimma', 'dsOmics' or 'dsPLINK' ")
  }
  
  pvals <- Reduce(ff, x)
  ii <- grep("id|P.Value", colnames(pvals))
  pvals <- pvals[,ii]
  colnames(pvals)[-1] <- names(x)
  p.meta <- unlist(apply(pvals[,-1], 1, function(x) metap::sumlog (x)$p))
  ans <- tibble::add_column(pvals, p.meta=p.meta)%>%arrange(p.meta)
  
  return(ans)
}
