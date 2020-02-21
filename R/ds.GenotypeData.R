##' @title Creates a GenotypeData object at each server
##' 
##' @description The GenotypeData class is a container for storing genotype data from a GWAS 
##' together with the metadata associated with the subjects (i.e. phenotypes and covariates) and SNPs involved in the study.
##' 
##' @param x a \code{GdsGenotypeReader} object (see GWASTools). It is the object that is obtained in the opal servers
##' from a resource of type VCF2GDS
##' @param covars a data.frame or a tibble having the metadata of samples (i.e. phenotypes and/or covariates)
##' @param data name of the DataSHIELD object to which the genotype (snpMatrix) and phenotypic data (data.frame) has been assigned
##' @param datasources a list of \code{DSConnection-class} objects obtained after login. If the <datasources> the default set of connections will be used: see \code{datashield.connections_default}.
##' @export
##' @examples
##' 

ds.GenotypeData <- function(x, covars, columnId, newobj.name = NULL, datasources=NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  if (is.null(newobj.name)) {
    newobj.name <- paste0(x, ".Data")
  }
  
  cally <- paste0("GenotypeDataDS(", x, "," , covars, ",", columnId, ")")
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
