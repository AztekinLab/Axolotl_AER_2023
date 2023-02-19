###
## USAGE
##
## input(mandatory)
## - df: the seurat FindMarker output
##
##

library(ggplot2)
library(dplyr)

volcano_plot <- function(df, logFC_threshold = 0.25, padj_threshold = 0.05,
                         label_up = "Up", label_down = "Down", label_no = "Not Sig.",
                         col_up = "#BC3C28", col_down = "#0072B5", col_no = "lightgrey",
                         title = "Volcano plot", xtitle = "log(FoldChange)", ytitle = "-log10(Adjust P-value)",
                         ytext_size = 10, xtext_size = 10, ytitle_size = 10, xtitle_size = 10,
                         label_top = 10, label_size = 4, dot_size = 1, box.padding = 0.5,
                         return_p = TRUE) {
  df$gene <- rownames(df)

  data <- data.frame(
    logFC = df$avg_log2FC,
    padj = df$p_val_adj,
    gene = df$gene
  )

  data$sig[(data$padj > padj_threshold) | (abs(data$logFC) < logFC_threshold)] <- label_no
  data$sig[data$padj <= padj_threshold & data$logFC >= logFC_threshold] <- label_up
  data$sig[data$padj <= padj_threshold & data$logFC <= -logFC_threshold] <- label_down
  data$sig <- factor(data$sig, levels = c(label_up, label_down, label_no))

  data$logpadj <- (-1 * log10(data$padj))
  data$logpadj[data$logpadj == Inf] <- 300


  # set xlim limits
  x_lim <- max(abs(data$logFC))


  # plot
  p <- ggplot(data, aes(logFC, logpadj, color = sig)) +
    geom_point(size = dot_size) + # scale_size_manual(values = c(1.5,2.5))+
    # xlim(-1,1) +  ylim(0,310)+
    labs(title = title, x = xtitle, y = ytitle) +
    scale_color_manual(values = c(col_up, col_down, col_no)) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = 4) +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = 4) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(size = 0),
      axis.text.y = element_text(size = ytext_size),
      axis.text.x = element_text(size = xtext_size),
      axis.title.y = element_text(size = ytitle_size),
      axis.title.x = element_text(size = xtitle_size),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    guides(colour = guide_legend(title = NULL))

  if (!is.null(label_top)) {
    if(is.numeric(label_top)){
label_gene <- data %>%
      group_by(sig) %>%
      top_n(label_top, abs(logFC)) %>%
      filter(sig != label_no)
      label_gene <- label_gene$gene
    }else{
      label_gene <- label_top
    }


    data$delabel <- NA
    data$delabel[data$gene %in% label_gene] <- data$gene[data$gene %in% label_gene]

    p <- p + ggrepel::geom_text_repel(
      label = data$delabel, size = label_size,
      box.padding = 0.5,
      max.overlaps = 50,
      # arrow = arrow(length = unit(0.02, "npc")),
      show.legend = FALSE
    )
  }

  return(p)
}
