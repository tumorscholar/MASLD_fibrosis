###############################################################################
# 20 - Extended Data Figure 1 - cell-type annotation dotplots (RNA + ADT)
###############################################################################

library(Seurat)
library(SeuratExtend)
library(dplyr)

source("00_nature_style.R")

seu_rds <- "/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds"

fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/ExtendedData"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

SeuObj <- readRDS(seu_rds)

Idents(SeuObj) <- "cell_type_with_cluster"

###############################################################################
# RNA markers
###############################################################################

rna_markers <- SeuObj@misc$markers

top5_rna <- rna_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

DefaultAssay(SeuObj) <- "RNA"

rna_features <- unique(top5_rna$gene)

p_rna <- DotPlot2(
  SeuObj,
  features = rna_features,
  color_scheme = "BuRd"
) +
  theme_nature(7) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

save_nature(
  p_rna,
  file.path(fig_dir, "ExtendedDataFig1_RNA_markers_v7"),
  15,
  30
)

###############################################################################
# ADT markers
###############################################################################

adt_markers <- SeuObj@misc$markersADT

top3_adt <- adt_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3) %>%
  ungroup()

DefaultAssay(SeuObj) <- "ADTonly"

adt_features <- unique(top3_adt$gene)

p_adt <- DotPlot2(
  SeuObj,
  features = adt_features,
  color_scheme = "BuRd"
) +
  theme_nature(7) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

save_nature(
  p_adt,
  file.path(fig_dir, "ExtendedDataFig1_ADT_markers_v7"),
  15,
  30
)

cat("\nDone Extended Data Figure 1.\n")
