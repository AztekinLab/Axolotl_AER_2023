library(dplyr)
library(MetaNeighbor)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(RColorBrewer)
library(stringi)

plot_AUROC_heatmap <- function(AUROC.scores) {
    r_keep <- rowSums(is.na(AUROC.scores)) != nrow(AUROC.scores)
    c_keep <- colSums(is.na(AUROC.scores)) != ncol(AUROC.scores)
    AUROC.scores <- AUROC.scores[r_keep, c_keep]


    ## row/col labels
    r_lab <- stri_replace_all_regex(rownames(AUROC.scores),
        pattern = c(".*_", "\\|.*"), replacement = c("", ""), vectorize = FALSE
    )
    c_lab <- r_lab


    # row/col annotation
    cols <- brewer.pal(5, 'Set2')
    names(cols) <- c('Human',"Mouse", "Chicken", "Frog", "Axolotl")
    sp <- gsub("_.*", "", rownames(AUROC.scores))

    c_anno <- HeatmapAnnotation(
        Species = anno_points(rep(1, times = nrow(AUROC.scores)),
            ylim = c(0, 1),
            size = unit(3, "mm"),
            height = unit(1, "mm"),
            border = FALSE,
            axis = FALSE,
            gp = gpar(col = cols[sp])
        )
    )


    r_anno <- rowAnnotation(
        Species = anno_points(rep(1, times = nrow(AUROC.scores)),
            ylim = c(0, 1),
            size = unit(3, "mm"),
            width = unit(1, "mm"),
            border = FALSE,
            axis = FALSE,
            gp = gpar(col = cols[sp])
        )
    )



    p <- Heatmap(AUROC.scores,
        col = c("steelblue", "white", "red3"), border = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (AUROC.scores[i, j] > 0.9) {
                grid.text("*", x, y, gp = gpar(fontsize = 5))
            }
        },

        # for clustering
        row_split = sub(".*\\|", "", rownames(AUROC.scores)),
        column_split = sub(".*\\|", "", rownames(AUROC.scores)),
        clustering_method_rows = "ward.D2",
        clustering_distance_rows = "spearman",
        clustering_method_columns = "ward.D2",
        clustering_distance_columns = "spearman",
        show_row_dend = TRUE, show_column_dend = TRUE,
        row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
        # row_names_side = 'left', column_names_side = 'top',
        rect_gp = gpar(col = "white", lwd = 1),
        width = ncol(AUROC.scores) * unit(3, "mm"),
        height = nrow(AUROC.scores) * unit(3, "mm"),

        # annotation
        row_labels = r_lab, column_labels = c_lab,
        bottom_annotation = c_anno,
        right_annotation = r_anno
    )

    return(p)
}


run_MN_US <- function(seu_obj, study_id, cell_type, dir,
                      methods = c("HVG_MN", "HVG_Seurat", "PCA")) {
    if (methods == "HVG_MN") {
        se <- SummarizedExperiment(seu_obj@assays$RNA@counts,
            colData = seu_obj@meta.data
        )

        # get var_genes
        min_rec <- length(unique(colData(se)[[study_id]]))
        repeat{
            var_genes <- variableGenes(
                dat = se,
                exp_labels = colData(se)[[study_id]],
                min_recurrence = min_rec
            )
            min_rec <- min_rec - 1
            if (length(var_genes) > 200) {
                break
            }
        }
    } else if (methods == "HVG_Seurat") {
        se <- SummarizedExperiment(seu_obj@assays$RNA@counts,
            colData = seu_obj@meta.data
        )

        # get var_genes
        tryCatch(
            {
                var_genes <- VariableFeatures(seu_obj, assay = "integrated")
            },
            error = function(e) {
                var_genes <- VariableFeatures(seu_obj, assay = "RNA")
            }
        )
    } else if (methods == "PCA") {
        embed <- Embeddings(seu_obj)
        se <- SummarizedExperiment(t(embed),
            colData = seu_obj@meta.data
        )

        var_genes <- colnames(embed)
    } else {
        print("Methods should be one of the HVG_MN, HVG_Seurat, PCA! ")
    }

    print(paste0("len of hvg: ", length(var_genes)))
    quant <- quantile(1:length(var_genes))
    quant <- quant[-1]
    quant <- ceiling(quant)


    for (i in quant) {
        print(paste0("testing: ", i))
        AUROC.scores <- MetaNeighborUS(
            var_genes = var_genes[1:i],
            dat = se,
            study_id = colData(se)[[study_id]],
            cell_type = colData(se)[[cell_type]],
            fast_version = TRUE
        )

        fn <- paste0("AUROC_", methods, "_", i, ".csv")
        write.csv(AUROC.scores, file = file.path(dir, fn))

        p <- plot_AUROC_heatmap(AUROC.scores)

        fn <- paste0("AUROC_", methods, "_", i, ".pdf")
        pdf(file = file.path(dir, fn), width = 10, height = 10)
        print(p)
        dev.off()
    }
}

rm_cells <- function(seu_obj, cond1, cond2, min_cells = 10) {
    cond1 <- enquo(cond1)
    cond2 <- enquo(cond2)

    cond1_str <- quo_name(cond1)
    cond2_str <- quo_name(cond2)

    print(cond1)
    print(cond1_str)

    df <- seu_obj@meta.data %>%
        dplyr::count(!!cond1, !!cond2) %>%
        dplyr::filter(n < min_cells)

    rm_cell_df <- sapply(1:nrow(df), function(i) {
        seu_obj@meta.data[[cond1_str]] == df[i, cond1_str] &
            seu_obj@meta.data[[cond2_str]] == df[i, cond2_str]
    })

    rm_cells <- apply(rm_cell_df, 1, any)
    return(seu_obj[, !rm_cells])
}
