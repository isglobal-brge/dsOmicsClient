#' @title Plot the results of a \code{ds.PCA()}
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

plotPCA <- function(object, group = NULL, xcomp = 1, ycomp = 2, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  data <- dsBaseClient::ds.scatterPlot(x = paste0(object$pca_res, "$Dim.", xcomp), 
                               y = paste0(object$pca_res, "$Dim.", ycomp), 
                               type = "combine",  datasources = datasources, 
                               return.coords = TRUE)
  dev.off()
  
  if(!is.null(group)){
    cally <- paste0('plotPCADS(', object$pca_res, ', ', object$geno[1],', "', group, '")') # Assumes all geno objects share same phenotypes file
    res <- Reduce(c, datashield.aggregate(datasources, cally))
    grouping <- as.factor(res)
    plt <- ggplot2::ggplot(data.frame(data$pooled.coordinates)) +  
      ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = grouping))
    # data$pooled.coordinates <- cbind(data.frame(data$pooled.coordinates), grouping)
  } else {
    plt <- ggplot2::ggplot(data.frame(data$pooled.coordinates)) +  
      ggplot2::geom_point(ggplot2::aes(x = x, y = y))
  }
  plt
}
