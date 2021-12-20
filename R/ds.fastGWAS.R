#' @title Fast genome-wide association analysis (GWAS)
#'
#' @description Performs a distributed GWAS using the methodology proposed on the
#' paper https://doi.org/10.1186/1471-2105-14-166. This implementation has been adapted
#' to DataSHIELD as it can rely on a model 0 provided by the \code{ds.glm} function and
#' the other steps require mostly performing colSums, which are non-disclosive. However,
#' the results are not 100% precise, but they capture all heterogenity and individual
#' differences between cohorts, as the method does not rely on a meta-analysis,
#' the obtained results of this methodology are the same when used on a distributed
#' DataSHIELD environment than they would be having all the data gathered at the same computer.
#'
#' @param genoData \code{character vector} of objects on the server side object which is a container for storing genotype data
#' from a GWAS toghether with the metadata associated with the subjects (i.e. phenotypes and/or covariates)
#' and SNPs
#' @param formula \code{character or fomula} formula indicating the condition (left side) and other covariates to be adjusted for
##' (i.e. condition ~ covar1 + ... + covar2). The fitted model is: snp ~ condition + covar1 + ... + covarN
#' @param family \code{character} A description of the generalized linear model used. "binomial" is defined
#' for case/control studies. Quantitative traits can be analyzed by using "gaussian"
#' @param do.par \code{bool} (default \code{FALSE}) Whether to use parallelization on the servers, to do so the servers
#' have to have the package \code{doParallel} installed and run on a POSIX OS (Mac, Linux, Unix, BSD); Windows
#' is not supported. This parallelization computes in parallel each \code{genoData} object, therefore it is only useful
#' when the genoData is divided by chromosome.
#' @param n.cores \code{numeric} (default \code{NULL}) Numbers of cores to use when \code{do.par} is \code{TRUE}. If
#' \code{NULL} the number of cores used will be the maximum available minus one.
#' @param snpBlock \code{numeric} (default \code{20000L}) Block size for dividing the genotype data, it equals to the
#' number of SNPs used on each iteration, depending on the servers RAM it may perform faster using lower or greater
#' block sizes, do some testing to assess it.
#' @param datasources a list of \code{DSConnection-class} objects obtained after login. If the <datasources>
#' the default set of connections will be used: see \code{datashield.connections_default}.
#'
#' @return \code{data.frame} With the results of the GWAS
#' @export
#'

ds.fastGWAS <- function(genoData, formula, family, do.par = FALSE, n.cores = NULL, snpBlock = 20000L, datasources = NULL){

  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  # Get objective and covars of formula
  objective_variable <- all.vars(as.formula(formula))[1]
  covars <- all.vars(as.formula(formula))[-1]

  # Get amount of covars (k)
  k <- length(covars) + 1 # The +1 is to take into account the intercept

  # TODO Add check that covars and objective_variables are in the genoData objects
  # TODO Add dsBaseClient:::isDefined(datasources, object)
  # and dsBaseClient:::checkClass(datasources, object) for the genoData objects

  # Convert the geno gds files into geno block iterators. This is to operate in SNP blocks so there is no need
  # of loading the whole geno info to the RAM at the same time, this improves performance and allows low performance
  # machines to be able to run this analysis
  DSI::datashield.assign.expr(datasources, "S_fastGWAS",
                              paste0("fastGWAS_S(c(",paste(genoData, collapse = ","),"), 'geno', ", snpBlock, ")"))
  # Extract phenotype information from the geno objects
  DSI::datashield.assign.expr(datasources, "PHENO_fastGWAS",
                              paste0("fastGWAS_S(c(",paste(genoData, collapse = ","),"), 'pheno')"))

  # Remove individuals from the phenotype that have NAs on any of the objective variable or covars. This has to
  # be done because this method does not allow for NA values (on the pheno and geno)
  DSI::datashield.assign.expr(datasources, "PHENO_fastGWAS",
                              paste0("fastGWAS_PHENO_removeNAindiv(PHENO_fastGWAS, '", objective_variable, "', c('",
                                     paste(covars, collapse = "','"), "'))"))

  # Get the geno mean of each individual to impute the missing SNPs
  DSI::datashield.assign.expr(datasources, "Smeans_fastGWAS", paste0("fastGWAS_S_means(S_fastGWAS, PHENO_fastGWAS)"))

  # Get the model0
  if(family == "binomial"){
    # Binomial case
      mod0 <- ds.glm(formula, family = "binomial", data = "PHENO_fastGWAS", datasources = datasources)$coefficients[,1]
  } else if (family == "gaussian") {
      mod0 <- ds.glm(formula, family = "gaussian", data = "PHENO_fastGWAS", datasources = datasources)$coefficients[,1]
  } else {
    stop("Incorrect family selected, options are ['binomial', 'gaussian']")
  }

  # Get fitted values
  DSI::datashield.assign.expr(datasources, "fitted.values_fastGWAS",
                              paste0("fastGWAS_getFitted.values(PHENO_fastGWAS, c('",
                                     paste(covars, collapse = "','"), "'), c('",
                                     paste(names(mod0)[-1], collapse = "','"), "'), '",
                                     family, "', c(", paste(mod0, collapse = ","), "))"))

  # Get residual values
  DSI::datashield.assign.expr(datasources, "YY_fastGWAS",
                              paste0("fastGWAS_getResiduals(PHENO_fastGWAS, fitted.values_fastGWAS, '",
                                     objective_variable, "', '", family,"')"))

  # Get mean of residuals
  mean_yy <- ds.mean("YY_fastGWAS",
          type = "combine", datasources = datasources)$Global.Mean[1]

  # Substract mean_yy to the model (new.glm.obj) residuals
  ds.make(toAssign = paste0("YY_fastGWAS - ", mean_yy),
          newobj = "YC_fastGWAS",
          datasources = datasources)

  if(length(datasources) > 1){
    # Get the colsums and crossprod all at once
    fastGWAS_S1_S2_BetaDEN1 <- Reduce(c, DSI::datashield.aggregate(datasources,
                                                                   paste0("fastGWAS_ColSums(table1 = YC_fastGWAS, geno = S_fastGWAS, do.par = ",
                                                                          if(do.par){"TRUE"}else{"FALSE"}, ", n.cores = ",
                                                                          if(is.null(n.cores)){"NULL"}else{n.cores}, ", ",
                                                                          "type = 'GWASsums', pheno = PHENO_fastGWAS, means = Smeans_fastGWAS)")))
    # Extract the results into three different variables
    fastGWAS_S1_S2_BetaDEN1a <- data.frame(fastGWAS_S1_S2_BetaDEN1)
    cmd <- parse(text = paste0("c('",paste(colnames(fastGWAS_S1_S2_BetaDEN1a)
                                           [grepl("^crossprod", colnames(fastGWAS_S1_S2_BetaDEN1a))],
                                           collapse = "','"),"')"))
    b <- rowSums(fastGWAS_S1_S2_BetaDEN1a[,eval(cmd)])
    cmd <- parse(text = paste0("c('",paste(colnames(fastGWAS_S1_S2_BetaDEN1a)
                                           [grepl("^sums", colnames(fastGWAS_S1_S2_BetaDEN1a))],
                                           collapse = "','"),"')"))
    S1_fastGWAS <- rowSums(fastGWAS_S1_S2_BetaDEN1a[,eval(cmd)])
    cmd <- parse(text = paste0("c('",paste(colnames(fastGWAS_S1_S2_BetaDEN1a)
                                           [grepl("^squared_sums", colnames(fastGWAS_S1_S2_BetaDEN1a))],
                                           collapse = "','"),"')"))
    S2_fastGWAS <- rowSums(fastGWAS_S1_S2_BetaDEN1a[,eval(cmd)])
  } else {
    # Get the colsums and crossprod all at once
    fastGWAS_S1_S2_BetaDEN1 <- DSI::datashield.aggregate(datasources,
                                                         paste0("fastGWAS_ColSums(table1 = YC_fastGWAS, geno = S_fastGWAS, do.par = ",
                                                                if(do.par){"TRUE"}else{"FALSE"}, ", n.cores = ",
                                                                if(is.null(n.cores)){"NULL"}else{n.cores}, ", ",
                                                                "type = 'GWASsums', pheno = PHENO_fastGWAS, means = Smeans_fastGWAS)"))[[1]]
    # Extract the results into three different variables
    b <- fastGWAS_S1_S2_BetaDEN1$crossprod
    S1_fastGWAS <- fastGWAS_S1_S2_BetaDEN1$sums
    S2_fastGWAS <- fastGWAS_S1_S2_BetaDEN1$squared_sums
  }

  # Get number of individuals
  N_IND <- ds.length("YC_fastGWAS", "combine", datasources = datasources)[[1]]

  # Get DEN1
  DEN1 <- S2_fastGWAS - (S1_fastGWAS ^ 2) / N_IND

  # Get betas
  B <- b / DEN1

  # Get sum of YC_fastGWAS squared
  YC_sum_squared <- Reduce("+", DSI::datashield.aggregate(datasources,
                                                          "fastGWAS_ColSums(YC_fastGWAS, 'square_vect')"))

  # Compute sigma squared
  SIG <- (YC_sum_squared - B ^ 2 * DEN1) / (N_IND-k-2)

  # Compute error
  ERR <- sqrt(SIG * (1/DEN1))

  # Compute pvalues
  PVAL <- 2 * pnorm(-abs(B / ERR))

  # Get genotype information to complete GWAS results table to fit ds.GWAS output
  rs <- DSI::datashield.aggregate(datasources[1], paste0("getVariable(c(", paste(genoData, collapse = ","), "), 'snp.rs.id')"))[[1]]
  alleles <- DSI::datashield.aggregate(datasources[1], paste0("getVariable(c(", paste(genoData, collapse = ","), "), 'snp.allele')"))[[1]]
  position <- DSI::datashield.aggregate(datasources[1], paste0("getVariable(c(", paste(genoData, collapse = ","), "), 'snp.position')"))[[1]]
  chromosome <- DSI::datashield.aggregate(datasources[1], paste0("getVariable(c(", paste(genoData, collapse = ","), "), 'snp.chromosome')"))[[1]]
  reference_allele <- substring(alleles, 1, 1)
  alternate_allele <- substring(alleles, 3, 3)

  # Gather information and return ordered by p.value
  return(tibble(rs = rs, chr = chromosome, pos = position, p.value = PVAL, Est = B, Est.SE = ERR,
                reference_allele = reference_allele, alternate_allele = alternate_allele) %>% arrange(p.value))
}
