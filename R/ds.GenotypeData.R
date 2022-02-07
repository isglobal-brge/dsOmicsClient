#' @title Creates a GenotypeData object at each server
#'
#' @description The GenotypeData class is a container for storing genotype data from a GWAS
#' together with the metadata associated with the subjects (i.e. phenotypes and covariates)
#' and SNPs involved in the study.
#'
#' @details This function contains the functionality of re-encoding a column of the covars
#' table from a binomial (plus optional strings to be set to NA) to 0/1 encoding. This is
#' due because the subsequent analysis functionalities such as \code{ds.GWaS()} and \code{ds.glmSNP()}
#' require this encoding to study binomial outcomes. To use this funcionality, input the arguments
#' \code{case_control_column}, \code{case} and \code{control} (optionally \code{na_string} can be
#' used to specify strings that will be set as NA, this is useful for phenotypes that are encoded as:
#' Yes/No/Don't know/Don't answer/...). Input \code{NULL} to this arguments if this functionality is
#' not to be used.
#' 
#' @param x a \code{GdsGenotypeReader} object (see GWASTools). It is the object that is
#' obtained in the opal servers
#' from a resource of type VCF2GDS
#' @param covars a data.frame or a tibble having the metadata of samples (i.e. phenotypes and/or covariates)
#' @param columnId \code{character} Column of the covars that contains the IDs
#' @param sexId \code{character} (default \code{NULL}) Column of the covars that contains the sex phenotype
#' @param male_encoding \code{character} (default \code{"male"}) String used to encode the male sex
#' phenotype on the covars table
#' @param female_encoding \code{character} (default \code{"female"}) String used to encode the female sex
#' phenotype on the covars table
#' @param case_control_column \code{character} (default \code{NULL}) Column that holds the
#' case/control to relevel to 0/1
#' @param case \code{character} (default \code{NULL}) Encoding of the case of the \code{case_control_column}
#' @param control \code{character} (default \code{NULL}) Encoding of the control of the \code{case_control_column}
#' @param newobj.name \code{character} (default \code{NULL})
#' @param datasources a list of \code{DSConnection-class} objects obtained after login. If the <datasources> the default set of connections will be used: see \code{datashield.connections_default}.
#'
#' @export
#' @examples
#'

ds.GenotypeData <- function(x, covars, columnId, sexId = NULL, male_encoding = "male", female_encoding = "female",
                            case_control_column = NULL, case = NULL, control = NULL,
                            newobj.name = NULL, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if (is.null(newobj.name)) {
    newobj.name <- paste0(x, ".Data")
  }
  
  if(!is.null(case_control_column) & (is.null(case)|is.null(control))){
    stop('If a case/control column is to be releveled, please input both arguments
         [case] and [control]')
  }

  cally <- paste0("GenotypeDataDS(",
                  paste(x, covars,
                    if(is.null(columnId)){"NULL"}else{paste0("'", paste0(charToRaw(columnId), collapse = ""), "'")},
                    if(is.null(sexId)){"NULL"}else{paste0("'", paste0(charToRaw(sexId), collapse = ""), "'")},
                    paste0("'", paste0(charToRaw(male_encoding), collapse = ""), "'"),
                    paste0("'", paste0(charToRaw(female_encoding), collapse = ""), "'"),
                    if(is.null(case_control_column)){"NULL"}else{paste0("'", paste0(charToRaw(case_control_column), collapse = ""), "'")},
                    if(is.null(case)){"NULL"}else{paste0("'", paste0(charToRaw(case), collapse = ""), "'")},
                    if(is.null(control)){"NULL"}else{paste0("'", paste0(charToRaw(control), collapse = ""), "'")}, 
                    sep=","),
                  ")")

  DSI::datashield.assign(datasources,  symbol = newobj.name, as.symbol(cally))
  
  test.obj.name <- newobj.name
  calltext <- call("testObjExistsDS", test.obj.name)
  object.info <- DSI::datashield.aggregate(datasources, calltext)
  num.datasources <- length(object.info)
  obj.name.exists.in.all.sources <- TRUE
  obj.non.null.in.all.sources <- TRUE
  for (j in 1:num.datasources) {
    if (!object.info[[j]]$test.obj.exists) {
      obj.name.exists.in.all.sources <- FALSE
    }
    if (object.info[[j]]$test.obj.class == "ABSENT") {
      obj.non.null.in.all.sources <- FALSE
    }
  }
  if (obj.name.exists.in.all.sources && obj.non.null.in.all.sources) {
    return.message <- paste0("Data object <", test.obj.name, 
                             "> correctly created in all specified data sources")
  }
  else {
    return.message.1 <- paste0("Error: A valid data object <", 
                               test.obj.name, "> does NOT exist in ALL specified data sources")
    return.message.2 <- paste0("It is either ABSENT and/or has no valid content/class,see return.info above")
    return.message <- list(return.message.1, return.message.2)
    return.info <- object.info
    return(list(return.info = return.info, 
                return.message = return.message))
  }
}
