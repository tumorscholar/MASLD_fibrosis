###############################################################################
# 28 - Extended Data Figure 10 - trajectory (Slingshot) of tissue-shared expanded CD8 T cells
#
# Slingshot pseudotime is fitted across the tissue-shared expanded T-cell
# population (expanded_shared_class == "c_Shared_between_solid_tissues"), using
# tissue of origin as the cluster label, to describe the transcriptional ordering
# of shared clones across liver, SAT and VAT (supports the migration model).
#
# Input : expandedTcellFinalV7.rds (expanded T cells with clonal-sharing class)
# Output: PCA embedding with the fitted Slingshot curve, and the pseudotime
#         distribution by tissue, both in the shared Nature figure style.
###############################################################################
suppressPackageStartupMessages({
  library(Seurat); library(slingshot); library(SingleCellExperiment); library(ggplot2)
})
source("00_nature_style.R")                       # theme_nature(), tissue_cols, save_nature()

fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure3"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- subset the tissue-shared expanded T cells and rebuild a clean object ----
seurat_obj_subset <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinalV7.rds")
keep <- seurat_obj_subset$expanded_shared_class == "c_Shared_between_solid_tissues"
counts <- GetAssayData(subset(seurat_obj_subset, expanded_shared_class == "c_Shared_between_solid_tissues"),
                       slot = "counts")
expanded_tcells <- CreateSeuratObject(counts = counts)
expanded_tcells@meta.data <- seurat_obj_subset@meta.data[keep, ]
for (a in c("ADT", "HTO")) if (a %in% names(expanded_tcells@assays)) expanded_tcells[[a]] <- NULL

# ---- standard processing prior to trajectory fitting ------------------------
expanded_tcells <- NormalizeData(expanded_tcells, verbose = FALSE)
expanded_tcells <- FindVariableFeatures(expanded_tcells, verbose = FALSE)
expanded_tcells <- ScaleData(expanded_tcells, verbose = FALSE)
expanded_tcells <- RunPCA(expanded_tcells, npcs = 30, verbose = FALSE)

# ---- Slingshot: tissue of origin as the cluster label, PCA as the embedding --
sce <- as.SingleCellExperiment(expanded_tcells)
sce <- slingshot(sce, clusterLabels = expanded_tcells$Tissue, reducedDim = "PCA")
expanded_tcells$slingshot_pseudotime <- slingPseudotime(sce)[, 1]

tis_lev <- intersect(c("LIVER", "SAT", "VAT"), unique(as.character(expanded_tcells$Tissue)))

# ---- PCA embedding with the fitted trajectory curve (Nature style) ----------
pca <- as.data.frame(SingleCellExperiment::reducedDim(sce, "PCA")[, 1:2])
colnames(pca) <- c("PC1", "PC2")
pca$Tissue <- factor(as.character(expanded_tcells$Tissue), levels = tis_lev)
curves <- slingCurves(sce)
crv_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
  s <- curves[[i]]$s[curves[[i]]$ord, 1:2]
  data.frame(PC1 = s[, 1], PC2 = s[, 2], curve = i)
}))
p_pca <- ggplot(pca, aes(PC1, PC2, colour = Tissue)) +
  geom_point(size = 0.4, alpha = 0.7) +
  geom_path(data = crv_df, aes(PC1, PC2, group = curve), colour = "black", linewidth = 0.5) +
  scale_colour_manual(values = tissue_cols, name = "Tissue") +
  labs(title = "Trajectory of tissue-shared expanded CD8 T cells") +
  theme_nature(7)
save_nature(p_pca, file.path(fig_dir, "Fig3_trajectory_PCA_sharedTcells_v7"), 8.5, 7)

# ---- pseudotime distribution by tissue (Nature style) -----------------------
vdf <- data.frame(Tissue = factor(as.character(expanded_tcells$Tissue), levels = tis_lev),
                  pseudotime = expanded_tcells$slingshot_pseudotime)
p_vln <- ggplot(vdf, aes(Tissue, pseudotime, fill = Tissue)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.2, fill = "white") +
  scale_fill_manual(values = tissue_cols, guide = "none") +
  labs(x = NULL, y = "Slingshot pseudotime", title = "Pseudotime across tissues") +
  theme_nature(7)
save_nature(p_vln, file.path(fig_dir, "Fig3_trajectory_pseudotime_violin_sharedTcells_v7"), 6, 7)

cat("Done trajectory (tissue-shared expanded CD8 T cells): PCA + curve and pseudotime-by-tissue.\n")
# End of the script
