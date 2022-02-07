#' @title Subset ExpressionSet
#' 
#' @description Subset ExpressionSet using a categorical variable of the covariates as filter. Can also subset by 
#' complete cases
#'
#' @param eSet \code{character} ExpressionSet to subset
#' @param objective_variable \code{character} (default \code{NULL}) Name of the covariate on the ExpressionSet to use as filter
#' @param objective_value \code{character} (default \code{NULL}) Name of the value from the \code{objective_variable} to filter. 
#' The resulting subset will be the individuals that match this value. 
#' To put in in code, it can be represented as: \code{subset <- expressionSet[expressionSet$objective_variable == objective_value,]}
#' @param complete_cases \code{bool} (default \code{FALSE}) If \code{TRUE} only the complete cases will be included on the subset. 
#' This option can be used with \code{objective_variable} and \code{objective_value} or without them, if those arguments are not 
#' present, the subset will be the complete cases of the whole ExpressionSet.
#' @param newobj.name \code{character} (default \code{NULL}) Name of the subseted ExpressionSet. If \code{NULL}, the variable 
#' \code{"subsetted_ExpressionSet"} will be used.
#' @param datasources a list of \code{\link{DSConnection-class}} (default \code{NULL}) objects obtained after login
#'
#' @return This function does not have an output. It creates (or overwrites) a data frame on the study server.
#' @export

ds.subsetExpressionSet <- function(eSet, objective_variable = NULL, objective_value = NULL, complete_cases = FALSE,
                                   newobj.name = NULL, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if(is.null(newobj.name)){
    newobj.name <- "subsetted_ExpressionSet"
  }
  
  checkForExpressionSet(eSet, datasources)
  
  cally <- paste0("subsetExpressionSetDS(", eSet, ",'", 
                  if(!is.null(objective_variable)){objective_variable}else{"NULL"}, 
                  "','", 
                  if(!is.null(objective_value)){objective_value}else{"NULL"}, 
                  "',", complete_cases, ")")
  DSI::datashield.assign.expr(datasources, newobj.name, cally)
  
}