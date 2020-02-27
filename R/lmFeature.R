#'
#' @title Fits a LM to assess association between the features (outcome) and grouping variable (e.g case/control, condition, ...)
#' @description To be supplied
#' @param feature
#' @param vars
#' @param data
#' @param cellEstim
#' @param sva
#' @param datasources
#' 
#' @return a vector with effect estimates, standard error and associated p-value
#' @author Gonzalez, JR.
#'

lmFeature <- function(feature, vars, Set, cellCountsAdjust,
                        datasources){
  
  cally <- paste0("selFeatureDS(", Set, ",", deparse(feature), ",", 
                  deparse(vars), ")")  
  datashield.assign(datasources, 'dat', as.symbol(cally))
                      
  if (isTRUE(cellCountsAdjust)){
    ds.cbind(c('dat', 'cell.counts'), newobj='dat')
  }
  
  
  mm <- stats::as.formula(paste(feature, "~ ", 
                         paste(ds.colnames('dat')[[1]][-1], collapse="+")))
  mod <- ds.glm(mm, family='gaussian', data='dat', viewIter = FALSE)
  metrics <- base::as.data.frame(mod$coefficients[2, c(1,2,4)])
  names(metrics) <- feature
  return(metrics)
}
