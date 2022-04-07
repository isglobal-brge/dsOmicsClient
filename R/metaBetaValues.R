#' @title Meta-analysis of beta values
#' 
#' @param x a \code{dsOmics} object obtained from \code{ds.limma}, \code{ds.GWAS} or \code{ds.PLINK} functions applied o 2 or more studies
#'
#' @return a matrix with features p-values of each study and its combination
#' @export

metaBetaValues <- function(x){
  
  if(length(x)==1)
    stop('Nothing to be meta-analyzed. There is only a single study')
  
  if (inherits(x, 'dsLimma')){
    stop()
  } else if (inherits(x, 'dsGWAS')){
    ff <- function(x, y){
      inner_join(x, y, 
                 by=c('rs', 'chr', 'pos'))
    }
    grep_ii <- "rs|chr|pos|Est$|Est.[^SE]"
    grep_jj <- "rs|chr|pos|Est.SE"
  } else if (inherits(x, 'dsPLINK')){
    stop()
  } else {
    stop("Object should be of class 'dsLimma', 'dsOmics' or 'dsPLINK' ")
  }
  
  # Join the datasets
  joined_x <- Reduce(ff, x)
  
  # Get the betas
  ii <- grep(grep_ii, colnames(joined_x))
  betas <- joined_x[,ii]
  colnames(betas)[-c(1:3)] <- names(x)
  
  # Get the beta standard errors
  jj <- grep(grep_jj, colnames(joined_x))
  betas.se <- joined_x[,jj]
  colnames(betas.se)[-c(1:3)] <- names(x)
  
  # Mix betas and betas.se into a single data frame
  meta_data_og <- Reduce(ff, list(betas, betas.se))
  meta_data <- meta_data_og
  colnames(meta_data)[-c(1:3)] <- c(paste0("betas.", names(x)), paste0("SE.betas.", names(x)))
  
  # Calculate b.F
  for(i in 1:length(names(x))){
    betas.cohort <- as.symbol(paste0("betas.", names(x)[i]))
    SE.betas.cohort <- as.symbol(paste0("SE.betas.", names(x)[i]))
    colname1 <- as.symbol(paste0("b.F1.", names(x)[i]))
    colname2 <- as.symbol(paste0("b.F2.", names(x)[i]))
    meta_data <- meta_data %>% mutate({{colname1}} := {{betas.cohort}} / {{SE.betas.cohort}}^2)
    meta_data <- meta_data %>% mutate({{colname2}} := 1 / {{SE.betas.cohort}}^2)
  }
  
  meta_data <- meta_data %>% mutate(sum.b.F1 = rowSums(across(starts_with("b.F1."))))
  meta_data <- meta_data %>% mutate(sum.b.F2 = rowSums(across(starts_with("b.F2."))))
  meta_data <- meta_data %>% mutate(b.F = sum.b.F1 / sum.b.F2)
  
  # Calculate se.F
  meta_data <- meta_data %>% mutate(se.F = 1 / sqrt(sum.b.F2))
  
  # Calculate p.F
  meta_data <- meta_data %>% mutate(p.F = pchisq((b.F / se.F)^2, df = 1, lower = F))
  
  # Get OR, confidence interval and p-value
  meta_data <- meta_data %>% mutate(OR = exp(b.F), 
                       low = exp(b.F - 1.96*se.F), 
                       up = exp(b.F + 1.96*se.F), 
                       pval = p.F)
  
  # Return results
  results <- meta_data %>% select(b.F, se.F, p.F, OR, low, up, pval) %>% 
    add_column(rs = joined_x$rs, .before = "b.F") %>%
    add_column(chr = joined_x$chr, .after = "rs") %>%
    add_column(pos = joined_x$pos, .after = "chr") %>%
    arrange(pval)
  # class(results) <- c(class(results), "GWASBetaMetaAnalysis")
  return(results)
}
