#' @title SE-N Plot
#' 
#' @details Plot to reveal issues with trait transformations. A collection of studies with no 
#' issues will fall on a line.
#'
#' @param x \code{list of data.frames (output of ds.GWAS)} List of tables (default output of \code{ds.GWAS}).
#' @param se_column \code{character} (default \code{"Est.SE"}) Name of the column containing the SE.
#' @param n_column \code{character} (default \code{"n.obs"}) Name of the column containing the number of observations.
#'
#' @return A ggplot object
#' @export

seNPlot <- function(x, se_column = "Est.SE", n_column = "n.obs"){
  
  # Check that all column names are present on each table, if present: continue, if not: throw error
  if(!all(unlist(lapply(x, function(x){all(c(se_column, n_column) %in% colnames(x))})))){
    `%!in%` <- Negate(`%in%`)
    bad_colnames <- unique(unlist(lapply(x, function(x){
      c(se_column, n_column)[which(c(se_column, n_column) %!in% colnames(x))]
    })))
    stop('[', paste0(bad_colnames, collapse = ", "), '] column(s) not found on the input object.')
  }
  
  medSEinv <- lapply(x, function(x){
    1/median(x[[se_column]])
  })
  Nmax <- lapply(x, function(x){
    max(x[[n_column]])
  })
  medSEinv <- unlist(medSEinv)
  Nmax <- unlist(Nmax)
  data <- data.frame(medSEinv, Nmax) %>% tibble::rownames_to_column("cohort")
  
  ggplot(data, aes(x = sqrt(Nmax), y = medSEinv)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    ylab("1/median(SE)") +
    xlab(expression(sqrt(N[max])))
}
