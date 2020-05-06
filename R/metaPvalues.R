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
#' @import dplyr
#' @import metap 

metaPvalues <- function(x,  ...){
  
  if(length(x)==1)
    stop("Nothing to be meta-analyzed. There is only a single study")
  
  if (inherits(x, "dsLimma")){
    ff <- function(x, y){
      inner_join(x%>%select("id", "P.Value"), 
                 y%>%select("id", "P.Value"), 
                 by="id")
    }
  } else if (inherits(x, "dsGWAS")){
    ff <- function(x, y){
      inner_join(x%>%select("variant.id", "Score.pval"), 
                 y%>%select("variant.id", "Score.pval"), 
                 by="id")
    }
  } else{
      stop("Object should be of class 'dsLimma', 'dsOmics' or 'dsPLINK'")
  }
  
  pvals <- Reduce(ff, x)
  colnames(pvals)[-1] <- names(x)
  p.meta <- unlist(apply(pvals[,-1], 1, function(x) sumlog (x)$p))
  ans <- tibble(pvals, p.meta=p.meta)%>%arrange(p.meta)
  
  return(ans)
}
