#' @title PLINK arguments in the DataSHIELD required format
#' 
#' @description Writes PLINK commands in the DataSHIELD required format (e.g. remove --)
#' 
#' 
#' @param args PLINK arguments
#' 
#' @return ...
#' 
#' @author Yannick Marcon
#'
#'
tolist <- function(args) {
  tokens <- strsplit(args, split = " ")[[1]]
  currentKey <- NULL
  params <- list()
  for (i in 1:length(tokens)) {
    token <- tokens[[i]]
    if (startsWith(token, "--")) {
      if (!is.null(currentKey)) {
        params <- append(params, currentKey)
        currentKey <- gsub("--", "", token)
      } else {
        currentKey <- gsub("--", "", token)
      }
    } else if (!is.null(currentKey)) {
      params[currentKey] <- token
      currentKey <- NULL
    }
  }
  if (!is.null(currentKey)) {
    params <- append(params, currentKey)
  }
  params
}
