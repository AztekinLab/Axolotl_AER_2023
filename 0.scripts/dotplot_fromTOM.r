library(reshape2)

dotplot <- function(countn,
                    genes = rownames(countn),
                    cells = 1:length(colnames(countn)),
                    condition, conditionList = "none", norm = "max",
                    aspect.ratio = 1, xAngle = 45, yAngle = 0,
                    xSize = 10, ySize = 10,
                    condition_plot = NULL,
                    xAdj = 1, plot = TRUE, title = "title", return = FALSE,
                    normcells = 1:length(colnames(countn)), dot.scale = 2,
                    collow = "lightgrey", colhigh = "blue") {
  genes <- rev(genes)


    subset <- log10(countn[na.omit(match(genes, rownames(countn))), cells] + 1)
    nonzero <- countn[na.omit(match(genes, rownames(countn))), cells] > 0

  # rm zero expression
  # nonzero <- countn[na.omit(match(genes, rownames(countn))), cells] > 0

  avgd <- (t(aggregate(t(as.matrix(subset)), list(condition[cells]), mean)))
  clusters <- avgd[1, ]
  avgd <- avgd[-1, , drop = FALSE]

  if (is.vector(avgd)) {
    nr <- 1
  } else {
    nr <- nrow(avgd)
  }

  avgd <- matrix(as.numeric(unlist(avgd)), nrow = nr)
  rownames(avgd) <- rownames(subset)
  colnames(avgd) <- clusters

  pct <- (t(aggregate(t(as.matrix(nonzero)), list(condition[cells]), mean)))
  clusters <- pct[1, ]
  pct <- pct[-1, , drop = FALSE]

  if (is.vector(pct)) {
    nr <- 1
  } else {
    nr <- nrow(pct)
  }

  pct <- matrix(as.numeric(unlist(pct)), nrow = nr)
  rownames(pct) <- rownames(nonzero)
  colnames(pct) <- clusters


  if (norm == "max") {
    avgd_n <- (avgd / apply(avgd, 1, function(x) max(x)))
    # modified by Jixing, 2022.11.15
    # to keep genes of zero expression
    avgd_n[which(is.na(rowSums(avgd_n))), ] <- 0
  } else if (norm == "sum") {
    avgd_n <- avgd / rowSums(avgd)
  } else if (norm == "totalMax") {
    subsetMax <- log10(countn[na.omit(match(genes, rownames(countn))), normcells] + 1)
    avgdMax <- (t(aggregate(t(as.matrix(subsetMax)), list(condition[normcells]), mean)))
    clustersMax <- avgdMax[1, ]
    avgdMax <- avgdMax[-1, ]
    avgdMax <- matrix(as.numeric(unlist(avgdMax)), nrow = nrow(avgdMax))
    rownames(avgdMax) <- rownames(subsetMax)
    colnames(avgdMax) <- clustersMax
    avgd_n <- (avgd / apply(avgdMax, 1, function(x) max(x)))
  } else {
    avgd_n <- avgd
  }

  if (conditionList != "none") {
    avgd_n <- avgd_n[, match(conditionList, colnames(avgd_n)), drop = F]
    pct <- pct[, match(conditionList, colnames(pct)), drop = F]
  }

  if (dim(avgd_n)[2] > 1) {
    if (length(which(is.na(rowSums(avgd_n)))) > 0) {
      pct <- pct[-which(is.na(rowSums(avgd_n))), ]
      avgd_n <- avgd_n[-which(is.na(rowSums(avgd_n))), ]
    }
  }


  melted_ <- melt(t(avgd_n))
  pct_ <- melt(t(pct))
  melted_$pct <- 100 * pct_$value

  if (!is.null(condition_plot)) {
    melted_ <- melted_[melted_$Var1 %in% condition_plot, ]
  }

  p <- ggplot(data = melted_) +
    theme_bw() +
    geom_point(mapping = aes(x = Var1, y = Var2, size = pct, col = value)) +
    scale_color_gradient(low = collow, high = colhigh) +
    scale_size(range = c(0, dot.scale), limits = c(0, 100)) +
    theme(
      axis.text.x = element_text(angle = xAngle, hjust = xAdj, size = xSize),
      axis.text.y = element_text(angle = yAngle, size = ySize),
      aspect.ratio = aspect.ratio,
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    theme(panel.background = element_blank()) +
    theme(panel.grid = element_blank())

  p <- p + ggtitle(title)
  if (plot == TRUE) {
    print(p)
  }
  if (return) {
    return(p)
  }
}
