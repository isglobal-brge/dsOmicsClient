#' @title P-Z plot
#' 
#' @details Plot to reveal issues with beta estimates, standard errors and P values.
#' The plots compare P values reported in the association GWAS with P values calculated from 
#' Z-statistics (P.ztest) derived from the reported beta and standard error. A study with no issues 
#' should show perfect concordance. It will generate a plot per study, where concordance of all the study SNPs 
#' will be displayed.
#'
#' @param x \code{data.frame / list of data.frames (output of ds.GWAS)} Single table of results or list of tables 
#' (default output of \code{ds.GWAS}).
#' @param beta_column \code{character} (default \code{"Est"}) Name of the column containing the betas.
#' @param se_column \code{character} (default \code{"Est.SE"}) Name of the column containing the SE.
#' @param pval_column \code{character} (default \code{"p.value"}) Name of the column containing the p-values.
#'
#' @return A ggplot object
#' @export

pzPlot <- function(x, beta_column = "Est", se_column = "Est.SE", pval_column = "p.value"){
  
  # Single cohort
  if(inherits(x, "data.frame")){
    # Check that all column names are present, if present: plot, if not: throw error
    if(all(c(beta_column, se_column, pval_column) %in% colnames(x))){
      return(.pzPlotgg(x, beta_column, se_column, pval_column))
    } else {
      `%!in%` <- Negate(`%in%`)
      bad_colnames <- c(beta_column, se_column, pval_column)[which(c(beta_column, se_column, pval_column) %!in% colnames(x))]
      stop('[', paste0(bad_colnames, collapse = ", "), '] column(s) not found on the input object.')
    }
  }
  # Multi-cohort
  if(inherits(x, "list")){
    # Check that all column names are present on each table, if present: plot, if not: throw error
    if(all(unlist(lapply(x, function(x){all(c(beta_column, se_column, pval_column) %in% colnames(x))})))){
      pzPlots <- lapply(names(x), function(y){
        .pzPlotgg(x[[y]], beta_column, se_column, pval_column, y)
      })
      return(gridExtra::grid.arrange(grobs = pzPlots, ncol = 2))
    } else {
      `%!in%` <- Negate(`%in%`)
      bad_colnames <- unique(unlist(lapply(x, function(x){
        c(beta_column, se_column, pval_column)[which(c(beta_column, se_column, pval_column) %!in% colnames(x))]
      })))
      stop('[', paste0(bad_colnames, collapse = ", "), '] column(s) not found on the input object.')
    }
    
  }
}

#' @title Auxiliary function for \code{pzPlot}
#'
#' @export

.pzPlotgg <- function(x_single, beta_column, se_column, pval_column, title = NULL){
  ggplot(data = x_single) +
    geom_point(aes_string(y = paste0("-log10(", pval_column, ")"),
                          x = paste0("-log10(2*pnorm(abs(", beta_column, "/", se_column, "),lower.tail=FALSE))"))) +
    xlab(expression(-log[10](P.ztest))) + ylab(expression(-log[10](P))) + ggtitle(title)
}
