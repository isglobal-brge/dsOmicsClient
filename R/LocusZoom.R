#' @title LocusZoom
#' 
#' @param gwas_results \code{data.frame} With atleast three columns: Chromosome, position and p-value
#' @param range \code{numeric} (default \code{5e+05}) Range of the plot in genome positions.
#' @param position_zoom \code{numeric} (default \code{NULL}) Center position of the plot if supplied. If \code{NULL}, 
#' the center position is the lowerst p-value of the \code{gwas_results} table.
#' @param use_biomaRt \code{bool} (default \code{TRUE}) Set to \code{TRUE} to use biomaRt to retrieve the annotation of the 
#' genes. If set to \code{FALSE} TxDb.Hsapiens.UCSC.hgXX.knownGene will be used. Use the argument \code{TxDb_version} to 
#' select between hg19 or hg18. Note that when using TxDb the correspondent packages needs to be installed, providing 
#' a internet-free annotation option.
#' @param TxDb_version \code{character} (default \code{hg18}) TxDb version to use. This argument is only used when 
#' \code{use_biomaRt} is set to \code{FALSE}.
#' @param GRCh \code{numeric} (default \code{NULL}) GRCh version to connect to if not the current GRCh38, currently this can 
#' only be 37, if is different from 37 the GRCh38 version will be used. (Parameter passed to biomaRt::useEnsembl). 
#' This argument is only used when \code{use_biomaRt} is set to TRUE.
#' @param pvalue \code{character} (default \code{"p.value"}) Column name of \code{gwas_results} that contains the p-values.
#' @param chromosome \code{character} (default \code{"chr"}) Column name of \code{gwas_results} that contains 
#' the chromosome names.
#' @param rs \code{character} (default \code{"rs"}) Column name of \code{gwas_results} that contains 
#' the rs ID of the SNPs. This column will only be used if \code{snp_threshold} is different than \code{NULL}, therefore 
#' there is no need of supplying this argument if there is no interest to plot the labels of the rs IDs.
#' @param position \code{character} (default \code{"pos"}) Column name of \code{gwas_results} that contains the positions.
#' @param range_filter \code{numeric} (default \code{NULL}). If not \code{NULL}, only the annotations with a range larger 
#' than \code{range_filter} will be displayed.
#' @param draw_genes \code{bool} (default \code{TRUE}). Select if the genes of the region plot are to be plotted below the p-values 
#' plot.
#' @param snp_threshold \code{numeric} (default \code{5}) Threshold of -log_{10}p-value to draw labels with the SNP IDs. 
#' If \code{NULL} no labels will be drawn.
#'
#' @return
#' @export

LocusZoom <- function(gwas_results, range = 5e+05, position_zoom = NULL, use_biomaRt = TRUE, TxDb_version = c("hg18", "hg19"),
                      GRCh = NULL, pvalue = "p.value", chromosome = "chr", position = "pos", rs = "rs", range_filter = NULL,
                      draw_genes = TRUE, snp_threshold = 5){
  require("magrittr")
  require("tidyverse")
  require("patchwork")
  TxDb_version <- match.arg(TxDb_version)
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
  # Calcuate -log10(pvalue)
  dat.bmi.sel.region <- dat.bmi.sel.region %>% mutate(pval_calculated = -log10(eval(as.symbol(pvalue)))) 
  
  # Generate plot of the GWAS
  p1 <- ggplot(data = dat.bmi.sel.region) + 
    geom_point(aes(eval(as.symbol(position)), pval_calculated), shape = 1) + theme_bw() +
    labs(subtitle = paste("Chromosome", 
                          sel.chr, "from", format((sel.pos - range), big.mark = "'"), 
                          "to", format((sel.pos + range), big.mark = "'"), "bp")) + 
    ylab(expression(-log[10](p-value))) + theme(text = element_text(size = 11))
  
  # Put labels + threshold line of top SNPs if snp_threshold is different from NULL
  if(!is.null(snp_threshold)){
    require(ggrepel) # package to add labels without overlapping
    p1 <- p1 + geom_hline(yintercept = snp_threshold, color = "red") +
      geom_label_repel(data=dat.bmi.sel.region[dat.bmi.sel.region$pval_calculated >= snp_threshold,], 
                     aes_string(label=rs, x=position, y="pval_calculated"), size=3, force=1.3)
  }
  
  # Get genes information on the plotted region
  if(use_biomaRt & draw_genes){
    gene.ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = GRCh)
    
    out.bm.genes.region <- biomaRt::getBM(
      attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'), 
      filters = c('chromosome_name','start','end'), 
      values = list(sel.chr, sel.pos - range, sel.pos + range), 
      mart = gene.ensembl)
  } else if(!use_biomaRt & draw_genes) {
    if (TxDb_version == "hg19" & require("TxDb.Hsapiens.UCSC.hg19.knownGene")) {
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    else if (TxDb_version == "hg18" & require("TxDb.Hsapiens.UCSC.hg18.knownGene"))
      txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
    else {
      warning("Genes are not shown since TxDb database is not installed in your computer")
      draw_genes <- FALSE
    }
    if(draw_genes){
      # Create GRanges from chromosome, min and max positions of the region
      x1 <- paste0("chr", unlist(dat.bmi.sel.region[chromosome][1,]), ":", min(dat.bmi.sel.region[position]), "-", 
                   max(dat.bmi.sel.region[position]))
      x1 <- as(x1, "GRanges")
      
      # Find overlaps between txdb genes and GRanges
      allg <- genes(txdb)
      out.bm.genes.region <- subsetByOverlaps(allg, x1)
      if(!require(Homo.sapiens)){
        warning("Genes are not shown since Homo.sapiens database is not installed in your computer")
        draw_genes <- FALSE
      }
      out.bm.genes.region$external_gene_name <- mapIds(Homo.sapiens::Homo.sapiens, 
                                                       keys=out.bm.genes.region$gene_id,
                                                       keytype="ENTREZID",
                                                       column="SYMBOL")
      out.bm.genes.region <- GenomicRanges::as.data.frame(out.bm.genes.region)
      names(out.bm.genes.region) <- c("seqnames", "start_position", "end_position", "width", "strand", 
                                      "gene_id", "external_gene_name")
    }
  }
  
  # Get plot range
  plot.range <- if(draw_genes){
    c(min(sel.pos - range, out.bm.genes.region$start_position), 
      max(sel.pos + range, out.bm.genes.region$end_position))
  }else{
    c(sel.pos - range, sel.pos + range)
  }
  
  if(draw_genes){
    # Clean and reorder results. Extract info for the plot
    out.bm.genes.region <- out.bm.genes.region %>% mutate(external_gene_name = forcats::fct_reorder(external_gene_name, 
                                                                                                    start_position, .desc = TRUE))
    out.bm.genes.region <- out.bm.genes.region %>% filter(external_gene_name != "")
    out.bm.genes.region <- tryCatch({
      out.bm.genes.region %>% 
        mutate(gene_biotype_fac = forcats::fct_relevel(as.factor(gene_biotype), "protein_coding"), 
               external_gene_name = forcats::fct_reorder2(external_gene_name, start_position, gene_biotype_fac, .desc = TRUE))
    }, error = function(w){out.bm.genes.region})
    # Apply filter for small ranges, bypass if NULL
    if(!is.null(range_filter)){
      out.bm.genes.region <- out.bm.genes.region %>% filter(range_filter < end_position - start_position)
    }
    ## Generate the genes plot
    if(use_biomaRt) { # biomaRt database contains the gene_biotype_fac grouping, txdb does not
      p2 <- ggplot2::ggplot(data = out.bm.genes.region) + 
        geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, 
                           colour = gene_biotype_fac, group = gene_biotype_fac), size = 3) +
        geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, 
                      colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size = 3,
                  show.legend = F) +
        scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1)))
    } else {
      p2 <- ggplot2::ggplot(data = out.bm.genes.region) + 
        geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position), size = 3) +
        geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name), 
                  fontface = 2, alpha = I(0.7), hjust = "right", size = 3,
                  show.legend = F)
    }
    p2 <- p2 + coord_flip() + ylab("") + theme_bw() + ylim(plot.range) + 
      labs(color = "Gene Biotype") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(), 
            strip.text.y = element_text(angle = 0),
            legend.position="bottom", 
            panel.grid.major.y = element_blank()) + 
      guides(color = guide_legend(override.aes = list(size = 2))) + 
      expand_limits(y=c(-1, 1)) +
      theme(text = element_text(size = 11))
  }

  # Combine two plots using patchwork if draw_genes is TRUE
  p1 <- p1 + xlab("") + theme(axis.title.x=element_blank(), 
                              axis.text.x=element_blank(), 
                              axis.ticks.x=element_blank()) + xlim(plot.range)
  
  return(if(draw_genes){p1 + p2 + plot_layout(ncol = 1)}else{p1})
}
