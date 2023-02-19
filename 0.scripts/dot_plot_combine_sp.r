library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(cowplot)

ortholog <- read.table("/Users/jzhong/Desktop/Desktop–SV-96M-002/EPFL/project/axolotl_aer_final/rmd_pl_figs/0.gene.list/orthologs.txt", header = T, sep = " ")

## ========== function definition ================
Human_fun <- function(mat) {
    return(mat)
}

Mouse_fun <- function(mat) {
    mat$features.plot <- as.character(mat$features.plot)
    idx <- which(mat$features.plot %in% ortholog[, 4])
    mat$features.plot[idx] <- ortholog[match(mat$features.plot[idx], ortholog[, 4]), 2]
    no_idx <- which(!(mat$features.plot %in% ortholog[, 4]))

    if (length(no_idx) > 0) {
        mat$features.plot[no_idx] <- toupper(mat$features.plot[no_idx])
    }
    return(mat)
}

Chicken_fun <- function(mat) {
    mat$features.plot <- as.character(mat$features.plot)
    idx <- which(mat$features.plot %in% ortholog[, 6])
    mat$features.plot[idx] <- ortholog[match(mat$features.plot[idx], ortholog[, 6]), 2]
    no_idx <- which(!(mat$features.plot %in% ortholog[, 6]))

    if (length(no_idx) > 0) {
        mat$features.plot[no_idx] <- toupper(mat$features.plot[no_idx])
    }
    return(mat)
}


Frog_fun <- function(mat) {
    mat$features.plot <- as.character(mat$features.plot)
    mat$features.plot <- toupper(gsub("\\.[LS]$", "", mat$features.plot))
    mat <- mat %>%
        group_by(features.plot, id) %>%
        filter(avg.exp == max(avg.exp))
    return(mat)
}


Axolotl_fun <- function(mat) {
    mat$features.plot <- as.character(mat$features.plot)
    mat$features.plot <- toupper(gsub("\\.[0-9]+$", "", mat$features.plot))
    mat <- mat %>%
        group_by(features.plot, id) %>%
        filter(avg.exp == max(avg.exp))
    return(mat)
}


## ============ extract data for plotting =========

fetch_data <- function(seu, marker_sp_list,
                       group.by = group.by,
                       split.by = split.by,
                       idents = NULL) {
    species_ <- gsub("_.*", "", seu$orig.ident[1])
    print(seu$orig.ident[1])

    Idents(seu) <- group.by
    if (is.null(idents) || idents %in% unique(seu@active.ident)) {
        p <- DotPlot(seu,
            features = unique(marker_sp_list[[species_]]), assay = "RNA",
            idents = idents,
            group.by = group.by,
            split.by = split.by,
            cols = "RdPu",
            scale = FALSE
        )
        mat <- p$data
        sp_fun <- eval(parse(text = paste0(species_, "_fun")))
        mat <- sp_fun(mat)

        print(dim(mat))
        return(mat[, 1:5])
    } else {
        return(NULL)
    }
}


# scale by dataset
scaling_exp <- function(total_mat, scale.by) {
    for (i in unique(total_mat[[scale.by]])) {
        idx <- which(total_mat[[scale.by]] == i)

        x <- total_mat[idx, "avg.exp"]
        # total_mat[which(total_mat$id==i), 'avg.exp'] <- (x-min(x))/(max(x)-min(x))
        total_mat[idx, "avg.exp"] <- x / max(x)
    }
    return(total_mat)
}



mat_formating <- function(total_mat, le_celltype, le_gene){
  total_mat$CellType <- gsub("_.*", "", total_mat$id)
  total_mat <- total_mat[total_mat$CellType != "Delected", ]
  
  total_mat <- total_mat %>%
    group_by(CellType, features.plot) %>%
    summarise(
      avg.exp = max(avg.exp),
      pct.exp = max(pct.exp),
      sample = gsub("(.*_)(.*_)", "\\2", id)
    )
  
  
  total_mat <- scaling_exp(total_mat, scale.by = "sample")
  total_mat <- scaling_exp(total_mat, scale.by = "features.plot")
  
  
  
  total_mat$CellType <- factor(total_mat$CellType, levels = le_celltype)
  total_mat$features.plot <- factor(total_mat$features.plot,
                                    levels = le_gene
  )
  
  return(total_mat)
}







dot_plot <- function(total_mat, x, y) {
    p <- ggplot(total_mat, aes(!!enquo(x), !!enquo(y))) +
        geom_point(aes(size = pct.exp, color = avg.exp)) +
        geom_point(aes(size = ifelse(pct.exp == 0, NA, pct.exp)), shape = 21, colour = "black", stroke = 0.5) +
        guides(
            size = guide_legend(title = "Percent Expressed"),
            color = guide_colorbar(title = "Average Expression")
        ) +
        scale_color_distiller(
            palette = "RdPu", direction = 0, na.value = "grey",
            limits = range(0, 1)
        ) +
      scale_size_continuous(limits = range(0, 100)) + 
        theme_cowplot() +
        theme(
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(colour = "lightgrey", size = 0.2)
        )
    p
}
