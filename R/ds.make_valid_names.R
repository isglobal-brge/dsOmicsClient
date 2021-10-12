#' @title Rename column names and character variables
#'
#' @description Passes the column names and the character variables through \code{make.names}, this creates 
#' names without any special characters so there are no problems for accessing any columns or any category 
#' when using the DataSHIELD parser.
#'
#' @param x \code{character} Name of the table on the server to be modified
#' @param newobj \code{character} (default \code{NULL}) Name of the new object on the server. If \code{NULL} 
#' the input object \code{x} will be overwritten
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. 
#' If the \code{datasources} argument is not specified
#'
#' @export
#'

ds.make_valid_names <- function(x, newobj=NULL, datasources=NULL){
  
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }
  
  if(is.null(newobj)){
    newobj <- x
  }
  
  if(!inherits(x, "character")){
    stop('[', deparse(substitute(x)), '] is not of class "character"')
  }
  
  if (!all(unlist(dsBaseClient::ds.exists(x, datasources)))){
    stop('[', x, '] is not defined on all the study servers')
  }
  
  calltext <- call("make_valid_namesDS", x)
  DSI::datashield.assign.expr(datasources, newobj, calltext)
  
}