#' Title
#'
#' @param x 
#' @param genoData 
#' @param formula 
#' @param family 
#' @param datasources 
#'
#' @return
#' @export
#'
#' @examples
ds.fastGWAS <- function(genoData, formula, family, snpBlock, datasources = NULL){
  
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  
  # Get objective and covars of formula
  objective_variable <- all.vars(as.formula(formula))[1]
  covars <- all.vars(as.formula(formula))[-1]
  
  # Get amount of covars (k)
  k <- length(covars) + 1 # TODO el +1 es per tindre en compte el intercept, sha de treure si no hi ha intercep!!!
  
  # Put geno data on a table named S_fastGWAS (extract from GDSs) (also include the phenotypes
  # from the gds objects, to perform the models)
  # TODO fer un check que la funcio funcione be si nomes hi ha un genodata enlloc de varis!, aixo 
  # sha de mirar tant en aquesta funcion com en les altres !
  DSI::datashield.assign.expr(datasources, "S_fastGWAS", 
                              paste0("fastGWAS_S(c(",paste(genoData, collapse = ","),"), 'geno', ", snpBlock, ")"))
  DSI::datashield.assign.expr(datasources, "PHENO_fastGWAS", 
                              paste0("fastGWAS_S(c(",paste(genoData, collapse = ","),"), 'pheno')"))
  
  # Impute missing genotype with the missing SNP values with the sample mean of the observed genotypes. Also, 
  # remove individuals with incomplete phenotype objective variable ~ covariates (method does not allow NAs)
    # First, remove individuals from the PHENO
    DSI::datashield.assign.expr(datasources, "PHENO_fastGWAS", 
                                paste0("fastGWAS_PHENO_removeNAindiv(PHENO_fastGWAS, '", objective_variable, "', c('", 
                                       paste(covars, collapse = "','"), "'))"))
    # # Second, get means of genotype to impute missings
    DSI::datashield.assign.expr(datasources, "Smeans_fastGWAS", paste0("fastGWAS_S_impute(S_fastGWAS, PHENO_fastGWAS)"))
  
  if(family == "binomial"){
    # Binomial case 
      # With intercept
      mod0 <- ds.glm(formula, family = "binomial", data = "PHENO_fastGWAS", datasources = datasources)$coefficients[,1]
      
      # Without intercept
      # TODO
  } else if (family == "gaussian") {
      # With intercept
      mod0 <- ds.glm(formula, family = "gaussian", data = "PHENO_fastGWAS", datasources = datasources)$coefficients[,1]
      
      # Without intercept
      # TODO
  } else {
    stop("Incorrect family selected, options are ['binomial', 'gaussian']")
  }
  
  # Get fitted values
  DSI::datashield.assign.expr(conns, "fitted.values_fastGWAS", 
                              paste0("fastGWAS_getFitted.values(PHENO_fastGWAS, '", 
                                     family, "', ",
                                     paste(mod0, collapse = ","), ", '", 
                                     paste(c(objective_variable, covars), 
                                           collapse = "','"), "')"))
                                     
  # Get residual values
  DSI::datashield.assign.expr(conns, "YY_fastGWAS",
                              paste0("fastGWAS_getResiduals(PHENO_fastGWAS, fitted.values_fastGWAS, '", 
                                     objective_variable, "', '", family,"')"))
  
  # # Get colsums of the geno data
  #   # TODO maybe check before reducing that the colnames (snp rs) do coincide???
  # S1_fastGWAS <- Reduce("+", DSI::datashield.aggregate(datasources, 
  #                                                      "fastGWAS_ColSums(S_fastGWAS, PHENO_fastGWAS, 'std_geno', Smeans_fastGWAS)"))
  # 
  # # Get colsums squared of the geno data
  # S2_fastGWAS <- Reduce("+", DSI::datashield.aggregate(datasources, 
  #                                                      "fastGWAS_ColSums(S_fastGWAS, PHENO_fastGWAS, 'square_geno', Smeans_fastGWAS)"))
  
  # Get mean of residuals
  mean_yy <- ds.mean("YY_fastGWAS", 
          type = "combine", datasources = datasources)$Global.Mean[1]
  
  # Substract mean_yy to the model (new.glm.obj) residuals
  ds.make(toAssign = paste0("YY_fastGWAS - ", mean_yy),
          newobj = "YC_fastGWAS", 
          datasources = datasources)
  
  # Get the colsums and crossprod all at once
  fastGWAS_S1_S2_BetaDEN1 <- Reduce(c, DSI::datashield.aggregate(datasources,
                                       "fastGWAS_ColSums(table1 = YC_fastGWAS, geno = S_fastGWAS, 
                                       type = 'GWASsums', pheno = PHENO_fastGWAS, means = Smeans_fastGWAS)"))
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

  # Get number of individuals
  N_IND <- ds.length("YC_fastGWAS", "combine", datasources = datasources)[[1]]

  # S2 S1 and N_IND are on the client
  DEN1 <- S2_fastGWAS - (S1_fastGWAS ^ 2) / N_IND
  
  # Crossproduct of YC_fastGWAS, S_fastGWAS
  B <- b / DEN1

  # Get sum of YC_fastGWAS squared
  YC_sum_squared <- Reduce("+", DSI::datashield.aggregate(datasources, 
                                                          "fastGWAS_ColSums(YC_fastGWAS, NULL, 'square_vect')"))
  
  # Compute sigma squared
  SIG <- (YC_sum_squared - B ^ 2 * DEN1) / (N_IND-k-2)
  
  # Compute error
  ERR <- sqrt(SIG *(1/DEN1))
  
  # Compute pvalues
  PVAL <- 2 * pnorm(-abs(B / ERR))
  
  # Get genotype information to complete GWAS results table to fit ds.GWAS output
  rs <- do.call(c,lapply(genoData, function(x){
    DSI::datashield.aggregate(datasources, paste0("getVariable(", x, ", 'snp.rs.id')"))[[1]]
  }))
  alleles <- do.call(c,lapply(genoData, function(x){
    DSI::datashield.aggregate(datasources, paste0("getVariable(", x, ", 'snp.allele')"))[[1]]
  }))
  position <- do.call(c,lapply(genoData, function(x){
    DSI::datashield.aggregate(datasources, paste0("getVariable(", x, ", 'snp.position')"))[[1]]
  }))
  chromosome <- do.call(c,lapply(genoData, function(x){
    DSI::datashield.aggregate(datasources, paste0("getVariable(", x, ", 'snp.chromosome')"))[[1]]
  }))
  reference_allele <- substring(alleles, 1, 1)
  alternate_allele <- substring(alleles, 3, 3)
  
  # Gather information and return ordered by p.value
  return(tibble(rs = rs, chr = chromosome, pos = position, p.value = PVAL, Est = B, Est.SE = ERR,
                reference_allele = reference_allele, alternate_allele = alternate_allele) %>% arrange(p.value))
}
