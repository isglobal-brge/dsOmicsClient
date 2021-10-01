#' @title Lambda-N plot
#' 
#' @details Plot to reveal issues with population stratification. Population stratification can 
#' either inflate or deflate association P values and can be grasped by the genomic control (GC) 
#' inflation factor (lambdaGC). Studies with lambdaGC > 1.2 should be revised.
#' 
#' @param x \code{list of data.frames (output of ds.GWAS)} List of tables (default output of \code{ds.GWAS}).
#' @param pval_column \code{character} (default \code{"p.value"}) Name of the column containing the p-values.
#' @param n_column \code{character} (default \code{"n.obs"}) Name of the column containing the number of observations.
#' @param threshold \code{numeric} (default \code{1.2}) Threshold to be plotted. Cohorts with a lambda above 
#' this threshold will be labeled to identify them easily.
#'
#' @return
#' @export

lambdaNPlot <- function(x, pval_column = "p.value", n_column = "n.obs", threshold = 1.2){

  # Check that all column names are present on each table, if present: continue, if not: throw error
  if(!all(unlist(lapply(x, function(x){all(c(pval_column, n_column) %in% colnames(x))})))){
    `%!in%` <- Negate(`%in%`)
    bad_colnames <- unique(unlist(lapply(x, function(x){
      c(pval_column, n_column)[which(c(pval_column, n_column) %!in% colnames(x))]
    })))
    stop('[', paste0(bad_colnames, collapse = ", "), '] column(s) not found on the input object.')
  }
    
  # Get lambdaGC from each cohort
  lambdaGC <-lapply(x, function(x){
    chisq <- qchisq(1-x[[pval_column]],1)
    median(chisq)/qchisq(0.5,1)
  })
  Nmax <- lapply(x, function(x){
    max(x[[n_column]])
  })
  lambdaGC <- unlist(lambdaGC)
  Nmax <- unlist(Nmax)
  data <- data.frame(lambdaGC, Nmax) %>% tibble::rownames_to_column("cohort")
  # Plot lambdaGC
  ggplot(data) +
    geom_point(aes(x = sqrt(Nmax), y = lambdaGC)) +
    geom_hline(yintercept = 1.0) +
    geom_hline(yintercept = 1.2, colour = "red") + 
    geom_label_repel(data=data[data$lambdaGC >= {{threshold}},], 
                     aes(label=cohort, x=sqrt(Nmax), y=lambdaGC)) +
    ylab(expression(lambda[GC])) +
    xlab(expression(sqrt(N[max])))
  
}
