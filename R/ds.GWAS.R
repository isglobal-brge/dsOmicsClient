#' @title Genome-wide association analysis (GWAS)
#' 
#' @description Performs GWAS using GENESIS
#' @param genoData a \code{GenotypeData} object which is a container for storing genotype data
#' from a GWAS toghether with the metadata associated with the subjects (i.e. phenotypes and/or covariates)
#' and SNPs
#' @param model formula indicating the condition (left side) and other covariates to be adjusted for 
##' (i.e. condition ~ covar1 + ... + covar2). The fitted model is: snp ~ condition + covar1 + ... + covarN
#' @param family A description of the generalized linear model used. The defatul is "binomial" that is defined
#' for case/control studies. Quantitative traits can be analyzed by using "gaussian". Other values are accepted.
#' @param snpBlock an integer specifying the number of SNPs in an iteration block. See \code{GenotypeIterator} 
#' function in GWASTools package.
#'  
#' @param ... other arguments of \code{fitNullModel} function in GENESIS package
#' @return a matrix with SNPs ordered by p-values
#' 
#' @author Gonzalez, JR.
#'
#' @export 
#' 

ds.GWAS <- function(genoData, model, family="binomial", snpBlock=10000, datasources=NULL, ...){
  
  family.ini <- family
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  mt <- all.vars(model)
  variable_name <- mt[1] 
  if (length(mt)>1)
    covariable_names <- paste(mt[-1], collapse=",")
  else
    covariable_names <- NULL
  if(length(genoData) > 1){
    vars <- lapply(genoData, function(x){
      ds.varLabels(object = x, datasources = datasources)
    })
    if(!isTRUE(do.call(all.equal, vars))){
      stop('The phenotypes of the supplied GDS files do not match, make sure to use the same phenotypes table on the ds.GenotypeData function')
    } else {vars <- unlist(vars[[1]])}
    
  } else {vars <- unlist(ds.varLabels(object = genoData, datasources = datasources))}
  if(!all(mt %in% vars)){
    stop('[',mt[which(!(mt %in% vars))],
         '] Not a variable(s) name of [', genoData,']')
  }
  
  if(length(genoData) > 1){
    res <- lapply(genoData, function(x){
      cally <- paste0("GWASDS(", x, "," , deparse(variable_name), ",", deparse(covariable_names), 
                      ",", deparse(family.ini), ",", snpBlock, ")")
      DSI::datashield.aggregate(datasources, as.symbol(cally))
    })
    res <- do.call(c, res)
    keys <- unique(names(res))
    ans <- NULL
    for(i in keys){
      ii <- which(names(res) == i)
      ans[[i]] <- do.call(rbind, res[ii])
    }
  } else {
    cally <- paste0("GWASDS(", genoData, "," , deparse(variable_name), ",", deparse(covariable_names), 
                    ",", deparse(family.ini), ",", snpBlock, ")")
    ans <- DSI::datashield.aggregate(datasources, as.symbol(cally))
  }
  
  class(ans) <- c("dsGWAS", class(ans))
  
  return(ans)
  
}
