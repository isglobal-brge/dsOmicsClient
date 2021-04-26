#' @title Plot the results of a \code{ds.PCASNPS()}
#'
#' @param res \code{PCASNPS} Object returned from \code{ds.PCASNPS()}
#' @param xcomp \code{numeric} Principal component to plot on the X axis
#' @param ycomp \code{numeric} Principal component to plot on the Y axis
#'
#' @return \code{ggplot} object
#' @export

plotPCASNPS <- function(res, xcomp = 1, ycomp = 2){
  plt <- ggplot2::ggplot(res[[1]]) +  # Does plot with results from first study server!
    ggplot2::geom_point(aes_string(x = names(res[[1]])[xcomp], 
                                   y = names(res[[1]])[ycomp]))
  plt
}