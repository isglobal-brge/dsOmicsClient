#' @title Convert a Methylation ExpressionSet from EPIC to 450k Array (DataSHIELD Client Function)
#'
#' @description This function subsets a methylation ExpressionSet from an EPIC array (850k probes) to a 450k array or vice versa in a DataSHIELD environment.
#'
#' @param eSet_symbol A character string representing the server-side symbol of the ExpressionSet object containing either EPIC or 450k array data.
#' @param objective_array A character string specifying the target array type. Can be either "450k" or "epic".
#' @param new.obj A character string specifying the name of the new server-side object created. If NULL (default), the same name as `eSet_symbol` is used.
#' @param datasources A DataSHIELD connections object containing the information for one or multiple data sources.
#' 
#' @return A server-side symbol representing the ExpressionSet object containing the subset of probes based on the target array type.
#' @export
#'
#' @examples
#' # Assuming you have a DataSHIELD login and 'eSet' is available on the server
#' # Convert the eSet from EPIC array to 450k array
#' eSet_450k_symbol <- ds.methylation_array_convert(eSet_symbol = "eSet", objective_array = "450k")
#' # Convert the eSet from 450k array to EPIC array
#' eSet_epic_symbol <- ds.methylation_array_convert(eSet_symbol = "eSet", objective_array = "epic")

ds.methylation_array_convert <- function(eSet_symbol, objective_array, new.obj = NULL,
                                         datasources = NULL) {
  if (is.null(datasources)) {
    datasources <-  DSI::datashield.connections_find()
  }
  
  if (!is.character(eSet_symbol) || length(eSet_symbol) != 1) {
    stop("eSet_symbol must be a single character string")
  }
  
  if (!objective_array %in% c("450k", "epic")) {
    stop("objective_array must be either '450k' or 'epic'")
  }
  
  if (!is.null(new.obj) && (!is.character(new.obj) || length(new.obj) != 1)) {
    stop("new.obj must be either NULL or a single character string")
  }
  
  # Use eSet_symbol as the default name if new.obj is NULL
  if (is.null(new.obj)) {
    new.obj <- eSet_symbol
  }
  
  # Prepare the server function call
  func_call <- call("methylation_array_convertDS", as.name(eSet_symbol), objective_array)
  
  # Assign the result to a new server-side symbol
  DSI::datashield.assign.expr(
    conns = datasources, 
    symbol = new.obj, 
    expr = func_call
  )
}
