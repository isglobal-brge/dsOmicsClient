#' @title Plot the results of a \code{ds.PCASNPS()}
#'
#' @param res \code{PCASNPS} Object returned from \code{ds.PCASNPS()}
#' @param xcomp \code{numeric} Principal component to plot on the X axis
#' @param ycomp \code{numeric} Principal component to plot on the Y axis
#' @param group 
#' @param geno 
#' @param datasources 
#'
#' @return \code{ggplot} object
#' @export

plotPCASNPS <- function(res, group = NULL,
                        xcomp = 1, ycomp = 2, datasources = NULL){
  
  if(!is.null(group)){
    if (is.null(datasources)) {
      datasources <- DSI::datashield.connections_find()
    }
    cally <- paste0('plotPCASNPSDS(', res$set, ',"', group, '")')
    grouping <- as.factor(DSI::datashield.aggregate(datasources, as.symbol(cally))[[1]])
  }
  plt <- ggplot2::ggplot(res[[1]][[1]]) +  # Does plot with results from first study server!
    ggplot2::geom_point(aes_string(x = names(res[[1]][[1]])[xcomp], 
                                   y = names(res[[1]][[1]])[ycomp]))
  if(!is.null(group)){plt <- plt + ggplot2::aes(color = grouping)}
  plt
}
