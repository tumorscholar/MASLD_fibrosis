###############################################################################
# 16 - Build the CellChat input object (shared T cells + non-T clusters)
#
# Combines the tissue-shared expanded T cells (one unified "T cell" group) with
# all non-T-cell clusters, labelled by cell-type name rather than cluster index.
# The individual T / NK clusters are dropped, as the shared expanded T cells
# replace them as the single "T cell" group for the CellChat analysis.
#
# Input : SeuObjV7.rds (full annotated object; $cell_type, $Tissue, $Stage)
#         expandedTcellFinalV7.rds (shared-class annotation; c_Shared = shared)
# Output: SharedTcells_NonTcell_V7.rds  (Idents-ready $cluster_cellchat by name)
###############################################################################
library(Seurat)
library(dplyr)

seu_rds  <- "/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds"
exp_rds  <- "/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinalV7.rds"
out_rds  <- "/data/Blizard-AlazawiLab/rk/seurat/SharedTcells_NonTcell_V7.rds"

SeuObj       <- readRDS(seu_rds)
expanded_obj <- readRDS(exp_rds)

# barcodes of shared-between-solid-tissue expanded clones (unified "T cell")
shared_barcodes <- colnames(expanded_obj)[
 expanded_obj$expanded_shared_class == "c_Shared_between_solid_tissues"]
cat("shared expanded T cells:", length(shared_barcodes), "\n")

# T / NK lineage cell types to drop (replaced by the unified shared "T cell")
Tlineage <- c("Naive/TCM T cells", "Cytotoxic CD8 T cells", "Activated CD8 T cells",
              "Cytotoxic NK cells", "Memory T cells", "Activated T cells")

# label by cell-type NAME; T/NK -> drop; shared clones -> "T cell"
SeuObj$cluster_cellchat <- as.character(SeuObj$cell_type)
SeuObj$cluster_cellchat[SeuObj$cell_type %in% Tlineage]      <- NA
SeuObj$cluster_cellchat[colnames(SeuObj) %in% shared_barcodes] <- "T cell"

SeuObj <- subset(SeuObj, cells = colnames(SeuObj)[!is.na(SeuObj$cluster_cellchat)])

# "T cell" first, then the non-T cell types alphabetically
other <- sort(setdiff(unique(SeuObj$cluster_cellchat), "T cell"))
SeuObj$cluster_cellchat <- factor(SeuObj$cluster_cellchat, levels = c("T cell", other))

cat("\n== CellChat groups (by name) ==\n"); print(table(SeuObj$cluster_cellchat))
saveRDS(SeuObj, out_rds)
cat("\nSaved:", out_rds, "(", ncol(SeuObj), "cells )\n")
# End of the script