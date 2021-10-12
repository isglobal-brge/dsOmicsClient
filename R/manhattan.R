#' @title Mannhattan Plot
#' 
#' @description Performs Mannhattan plot for GWAS, DGE or EWAS analysis
#'
#' @param x an object of class dsGWAS 
#' @param featureCol column having feature name (e.g. SNP, CpG, gene, ..)
#' @param chrCol column having chromosome name
#' @param positionCol column having basepair position
#' @param sig genome-wide significant level. Default 5x10-8
#' @param sugg suggestive significant level. Default 1x10-6
#' @param threshold numeric value or NA for no label; label snps with pvalue < threshold. Default is NA
#' @param hlight list of features names to highlight. For no SNP higlighting, use NA (default)
#' @param col vector listing color/colors desired for manhattan plot. Default c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B")
#' @param ylims y limits to be visualized in -log10 scale. Default is c(1.5,9)
#' @param title plot title, with quotations (e.g. "My Manhattan Plot")
#' @param ... additional ggplot arguments (ylim, cex, etc. can be altered within the plot section of the function)
#'
#' @return a ggplot object having the Manhattan plot
#' 
#' @author modified from here: https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function
#'
#' @export
#' @import ggrepel
#' @import RColorBrewer

manhattan <- function(x, featureCol = 1, chrCol = 2, posCol = 3, pvalCol = 6, sig= 5e-8,
                      sugg = 1e-6, threshold=NA, hlight=NA, 
                      col=c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B"), 
                      ylims = c(1.5,9), title="Manhattan plot", ...){
  df <- x[ , c(featureCol, chrCol, posCol, pvalCol)]
  names(df) <- c("SNP", "CHR", "BP", "P")
  df <- dplyr::filter(df, CHR%in%c(1:22))
  df$BP <- as.numeric(df$BP)
  df$CHR <- as.numeric(df$CHR)

  df.tmp <- df %>% 
      # Compute chromosome size
      group_by(CHR) %>% 
      summarise(chr_len=max(BP)) %>% 
      
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len) %>%
      
      # Add this info to the initial dataset
      left_join(df, ., by=c("CHR"="CHR")) %>%
      
      # Add a cumulative position of each SNP
      arrange(CHR, BP) %>%
      mutate( BPcum=BP+tot) %>%
      
      # Add highlight and annotation information
      mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
      mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
    
    # get chromosome center positions for x-axis
    axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    
    ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
      # Show all points
      geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
      scale_color_manual(values = rep(col, 22 )) +
      
      # custom X axis:
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center, guide = guide_axis(check.overlap = T) ) +
      scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
      
      # add plot and axis titles
      ggtitle(paste0(title)) +
      labs(x = "Chromosome") +
      
      # add genome-wide sig and sugg lines
      geom_hline(yintercept = -log10(sig)) +
      geom_hline(yintercept = -log10(sugg), linetype="dashed") +
      
      # Add highlighted points
      #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
      
      # Add label using ggrepel to avoid overlapping
      geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
      
      # Custom the theme:
      theme_bw(base_size = 22) +
      theme( 
        plot.title = element_text(hjust = 0.5),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      ) 
    }
