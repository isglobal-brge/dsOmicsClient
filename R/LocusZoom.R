#' @title LocusZoom
#' 
#' @param gwas_results \code{data.frame} With atleast three columns: Chromosome, position and p-value
#' @param range \code{numeric} (default \code{5e+05}) Range of the plot in genome positions.
#' @param position_zoom \code{numeric} (default \code{NULL}) Center position of the plot if supplied. If \code{NULL}, 
#' the center position is the lowerst p-value of the \code{gwas_results} table.
#' @param GRCh \numeric (default \code{37}) GTCh version to use when querying biomaRt.
#' @param pvalue \character (default \code{"p.value"}) Column name of \code{gwas_results} that contains the p-values.
#' @param chromosome \character (default \code{"chr"}) Column name of \code{gwas_results} that contains 
#' the chromosome names.
#' @param position \character (default \code{"pos"}) Column name of \code{gwas_results} that contains the positions.
#'
#' @return
#' @export

LocusZoom <- function(gwas_results, range = 5e+05, position_zoom = NULL, GRCh = 37, pvalue = "p.value",
                      chromosome = "chr", position = "pos"){
  require("magrittr")
  require("tidyverse")
  require("patchwork")
  # Get range information from the results of the GWAS, this is the zoomed region to be plotted.
  # By default it takes the lower pvalue as point of interest, however a custom position can be supplied (position_zoom)
  if(is.null(position_zoom)){
    dat.bmi.sel <- gwas_results %>% dplyr::slice(which.min(eval(as.symbol(pvalue))))
    sel.chr <- unlist(dat.bmi.sel[, chromosome])
    sel.pos <- unlist(dat.bmi.sel[, position])
    dat.bmi.sel.region <- gwas_results %>% dplyr::filter(eval(as.symbol(chromosome)) == sel.chr, 
                                                  dplyr::between(eval(as.symbol(position)), sel.pos - range, sel.pos + range))
  } else {
    dat.bmi.sel <- gwas_results %>% dplyr::slice(which(eval(as.symbol(position)) == position_zoom))
    sel.chr <- unlist(dat.bmi.sel[, chromosome])
    sel.pos <- position_zoom
    dat.bmi.sel.region <- gwas_results %>% dplyr::filter(eval(as.symbol(chromosome)) == sel.chr, 
                                                         dplyr::between(eval(as.symbol(position)), sel.pos - range, sel.pos + range))
  }
  
  # Generate plot of the GWAS
  
  p1 <- ggplot(data = dat.bmi.sel.region) + 
    geom_point(aes(eval(as.symbol(position)), -log10(eval(as.symbol(pvalue)))), shape = 1) + theme_bw() +
    labs(subtitle = paste("Chromosome", 
                          sel.chr, "from", format((sel.pos - range), big.mark = "'"), 
                          "to", format((sel.pos + range), big.mark = "'"), "bp")) + 
    ylab(expression(-log[10](p-value))) + theme(text = element_text(size = 11))
  
  # Get genes information on the plotted region
  
  # TODO Option to use annotation from bioC package as per https://github.com/isglobal-brge/R-GADA/blob/master/R/plotCNVs.R
  # since biomaRt uses internet connection and can fail sometimes.
  
  gene.ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = GRCh)
  
  out.bm.genes.region <- biomaRt::getBM(
    attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'), 
    filters = c('chromosome_name','start','end'), 
    values = list(sel.chr, sel.pos - range, sel.pos + range), 
    mart = gene.ensembl)
  
  # Clean and reorder results. Extract info forthe plot
  
  out.bm.genes.region <- out.bm.genes.region %>% mutate(external_gene_name = forcats::fct_reorder(external_gene_name, 
                                                                                         start_position, .desc = TRUE))
  out.bm.genes.region <- out.bm.genes.region %>% filter(external_gene_name != "")
  plot.range <- c(min(sel.pos - range, out.bm.genes.region$start_position), 
                  max(sel.pos + range, out.bm.genes.region$end_position))
  out.bm.genes.region <- out.bm.genes.region %>% 
    mutate(gene_biotype_fac = forcats::fct_relevel(as.factor(gene_biotype), "protein_coding"), 
           external_gene_name = forcats::fct_reorder2(external_gene_name, start_position, gene_biotype_fac, .desc = TRUE))
  
  ## Generate the genes plot
  p2 <- ggplot2::ggplot(data = out.bm.genes.region) + 
    geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, 
                       colour = gene_biotype_fac, group = gene_biotype_fac)) +
    coord_flip() + ylab("") + theme_bw() +
    ylim(plot.range) + 
    geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, 
                  colour = gene_biotype_fac), size = 3.2, fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5,
              show.legend = F) + 
    labs(color = "Gene Biotype") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), 
          strip.text.y = element_text(angle = 0),
          legend.position="bottom", 
          panel.grid.major.y = element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 2))) + 
    expand_limits(y=c(-1, 1)) +
    scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1))) +
    theme(text = element_text(size = 11))
  
  # Combine two plots using patchwork
  
  p1 <- p1 + xlab("") + theme(axis.title.x=element_blank(), 
                              axis.text.x=element_blank(), 
                              axis.ticks.x=element_blank()) + xlim(plot.range)
  return(p1 + p2 + plot_layout(ncol = 1, heights = c(6, 6)))
}
