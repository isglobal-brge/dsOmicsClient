#' @title EAF plot
#' 
#' @details Plotting the reported EAFs against a reference set, such as from the HapMap or 
#' 1000 Genomes projects, or from one specific study, can help visualize patterns that 
#' pinpoint strand issues, allele miscoding or the inclusion of individuals whose self-reported 
#' ancestry did not match their genetic ancestry.
#'
#' @param x \code{data.frame / list of data.frames (output of ds.GWAS)} Single table of results or list of tables 
#' (default output of \code{ds.GWAS}).
#' @param reference \code{data.frame} Reference table that at least contains: SNP rs IDs and minor (or major) 
#' allele frequency of reference.
#' @param rs_x \code{character} (default \code{"rs"}) Name of the column that contains the rs IDs on \code{x}
#' @param rs_reference \code{character} Name of the column that contains the rs IDs on \code{reference}
#' @param freq_x \code{character} (default \code{"freq"}) Name of the column that contains the minor (or major) 
#' allele frequencies on \code{x}
#' @param freq_reference \code{character} Name of the column that contains the minor (or major) 
#' allele frequencies on \code{reference}
#' @param inverse_freq_x \code{bool} (default \code{TRUE}) Use \code{TRUE} if \code{x} contains major allele 
#' frequencies, use \code{FALSE} if \code{x} contains minor allele frequencies. The default is set to 
#' \code{TRUE} since the default output of \code{ds.GWAS} frequencies is major allele frequencies.
#' @param inverse_freq_reference \code{bool} (default \code{FALSE}) Use \code{TRUE} if \code{x} contains major allele 
#' frequencies, use \code{FALSE} if \code{x} contains minor allele frequencies. The default is set to 
#' \code{TRUE} since usually reference panels report minor allele frequencies.
#'
#' @return A ggplot object
#' @export

eafPlot <- function(x, reference, rs_x = "rs", rs_reference, freq_x = "freq", freq_reference, inverse_freq_x = T,
                    inverse_freq_reference = F){
  
  # Check columns do exist
    # Columns on x
      # Single cohort
      if(inherits(x, "data.frame")){
        if(!all(c(rs_x, freq_x) %in% colnames(x))){
          `%!in%` <- Negate(`%in%`)
          bad_colnames <- c(rs_x, freq_x)[which(c(rs_x, freq_x) %!in% colnames(x))]
          stop('[', paste0(bad_colnames, collapse = ", "), '] column(s) not found on the input object.')
        }
      }
  
      # Multi cohort
      if(inherits(x, "list")){
        if(!all(unlist(lapply(x, function(x){all(c(rs_x, freq_x) %in% colnames(x))})))){
          `%!in%` <- Negate(`%in%`)
          bad_colnames <- unique(unlist(lapply(x, function(x){
            c(rs_x, freq_x)[which(c(rs_x, freq_x) %!in% colnames(x))]
          })))
          stop('[', paste0(bad_colnames, collapse = ", "), '] column(s) not found on the input object.')
        }
      }
  
    # Columns on reference
    if(!all(c(rs_reference, freq_reference) %in% colnames(reference))){
      `%!in%` <- Negate(`%in%`)
      bad_colnames <- c(rs_reference, freq_reference)[which(c(rs_reference, freq_reference) %!in% colnames(reference))]
      stop('[', paste0(bad_colnames, collapse = ", "), '] column(s) not found on the reference object.')
    }
  
  # if freq_x and freq_reference are the same, modify one to avoid collisions on the merge!
  if(freq_x == freq_reference){
    colnames(reference)[which(colnames(reference) == freq_reference)] <- "freq_aux_alt"
    freq_reference <- "freq_aux_alt"
  }
  
  # Single cohort
  if(inherits(x, "data.frame")){
    merged_x <- merge(x, reference, by.x = rs_x, by.y = rs_reference)
    return(.eafPlot(merged_x, freq_x, freq_reference, inverse_freq_x, inverse_freq_reference))
  }
  # Multi-cohort
  if(inherits(x, "list")){
    merged_x <- lapply(x, function(x){merge(x, reference, by.x = rs_x, by.y = rs_reference)})
    eafPlots <- lapply(names(merged_x), function(y){
      .eafPlot(merged_x[[y]], freq_x, freq_reference, inverse_freq_x, inverse_freq_reference, y)
    })
    return(gridExtra::grid.arrange(grobs = eafPlots, ncol = 2))
  }
}

#' @title Auxiliary function for \code{eafPlot}
#'
#' @export

.eafPlot <- function(x, freq_x, freq_reference, inverse_freq_x, inverse_freq_reference, title = NULL){
  ggplot(x) + 
    geom_point(aes_string(y = if(inverse_freq_x){paste0("1-",freq_x)}else{freq_x}, 
                          x = if(inverse_freq_reference){paste0("1-",freq_reference)}else{freq_reference})) +
    xlab("EAF .reference") +
    ylab("EAF .study") +
    ggtitle(title)
}
