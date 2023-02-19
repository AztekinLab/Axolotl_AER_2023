library(Seurat)
source("/work/gr-aztekin/3.project/axolotl_AER_final/3.MN/pl_MN_utils.r")


# MetaNeighbor Celltype-wise
print("loading data")
multi <- readRDS("../2.clustering/integration/downsampled_seurat_ordered_cc_woC6.rds")
Idents(multi) <- "CellType2"

need <- c("AER", "CT", "BasalEctoderm", "Muscle")
multi_sub <- subset(multi, idents = need)

# rm celltype of less than 10 cells
# min_cells <- 10
# multi_sub <- rm_cells(multi_sub,
#     cond1 = orig.ident,
#     cond2 = CellType2,
#     min_cells = min_cells
# )

for (i in c("HVG_MN", "HVG_Seurat", "PCA")) {
    print(paste("Running", i))

    tryCatch(dir.create(i))
    run_MN_US(multi_sub, "orig.ident", "CellType2",
        dir = i,
        methods = i
    )
}



need <- c("AER", "BasalEctoderm")
multi_sub <- subset(multi, idents = need)


for (i in c("HVG_MN", "HVG_Seurat", "PCA")) {
    print(paste("Running", i))

    dir <- paste0("AER_vs_BasalEctoderm_", i)
    tryCatch(dir.create(dir))
    run_MN_US(multi_sub, "orig.ident", "CellType2",
        dir = dir,
        methods = i
    )
}
