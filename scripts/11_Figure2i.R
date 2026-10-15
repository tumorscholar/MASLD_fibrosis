###############################################################################
# Figure 2i - clonal occupancy per cell cluster
#
# For each cluster, the number of cells belonging to Small / Medium / Large /
# Hyperexpanded clones (stacked counts). Clone size is read from the annotated
# object produced by the shared TCR builder.
#
# Input : annotated Seurat object with clone-size annotation
# Output: occupancy table and Figure 2i (PDF + TIFF)
###############################################################################

.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))   # scRepertoire library location
library(scRepertoire)
library(Seurat)
library(tidyverse)
source("00_nature_style.R")            # theme_nature(), save_nature()

# ---- CONFIG (v7) ------------------------------------------------------------
clonal_rds <- "/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7_clonal.rds"
out_dir    <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure2"   # figure panels
tab_dir    <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure2"    # data tables
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(clonal_rds)

# ---- pick the cluster/annotation column to put on the x-axis ---------------
# EDIT if your Fig 1 annotation lives under a different column name.
cand <- c("cell_type", "CellType", "annotation", "seurat_clusters")
clus_col <- cand[cand %in% colnames(obj@meta.data)][1]
if (is.na(clus_col)) { obj$ident <- as.character(Idents(obj)); clus_col <- "ident" }
cat("x-axis (cluster) column:", clus_col, "\n")

# cloneSize must exist (added by the builder via combineExpression)
if (!"cloneSize" %in% colnames(obj@meta.data))
  stop("cloneSize not found - run the builder (combineExpression) first.")
cat("\ncloneSize levels:", paste(levels(obj$cloneSize), collapse = ", "), "\n")

# ---- export the occupancy table BEFORE plotting -----------------------------
occ_tab <- clonalOccupy(obj, x.axis = clus_col, proportion = FALSE,
                        na.include = FALSE, exportTable = TRUE)
cat("\n=== Clonal occupancy (cell counts per cluster x cloneSize) ===\n")
print(as.data.frame(occ_tab))
write.csv(occ_tab, file.path(tab_dir, "Fig2i_clonalOccupy_counts_V7.csv"), row.names = FALSE)

# ---- build the plot exactly like your original Fig 2i -----------------------
# Your figure: x = cluster IDs (C0..Cn), y = actual No of cells (not scaled),
# and ONLY the expanded bins are shown (Single and None are dropped). Stacking:
# Small at the bottom -> Hyperexpanded on top.
# Cluster number straight from seurat_clusters (0..n). Some annotation labels
# (e.g. "pDC") have no "(C##)" suffix, so DON'T parse the label - use the real
# cluster column; only fall back to label-parsing if seurat_clusters is absent.
if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  cl_num <- as.integer(as.character(obj$seurat_clusters))
} else {
  cl_num <- as.integer(sub(".*\\(C([0-9]+)\\).*$", "\\1", obj@meta.data[[clus_col]]))
}
if (all(is.na(cl_num)))
  stop("Could not derive cluster IDs. Columns available: ",
       paste(colnames(obj@meta.data), collapse = ", "))
cl_lev  <- paste0("C", 0:max(cl_num, na.rm = TRUE))      # keep empty clusters as ticks too

plot_df <- obj@meta.data %>%
  transmute(cluster = factor(paste0("C", cl_num), levels = cl_lev),
            cloneSize = as.character(cloneSize),
            bin = sub(" .*$", "", cloneSize)) %>%
  filter(bin %in% c("Small", "Medium", "Large", "Hyperexpanded")) %>%
  count(cluster, cloneSize, bin, name = "n")

# keep the FULL level labels (with the ranges, as in your legend), ordered
# Hyperexpanded -> Large -> Medium -> Small (first level = top of the stack).
ord_bins  <- c("Hyperexpanded", "Large", "Medium", "Small")
lab_order  <- unique(plot_df$cloneSize)[order(match(sub(" .*$", "", unique(plot_df$cloneSize)), ord_bins))]
plot_df$cloneSize <- factor(plot_df$cloneSize, levels = lab_order)

# your palette (sampled from the original figure)
bin_cols   <- c(Small = "#F08400", Medium = "#D83C18", Large = "#600000", Hyperexpanded = "#7830B4")
clone_cols <- setNames(unname(bin_cols[sub(" .*$", "", lab_order)]), lab_order)

write.csv(plot_df, file.path(tab_dir, "Fig2i_clonalOccupy_plotdata_V7.csv"), row.names = FALSE)

p <- ggplot(plot_df, aes(x = cluster, y = n, fill = cloneSize)) +
  geom_col(width = 0.8) +                                    # stacked counts, no outline
  scale_fill_manual(name = NULL, values = clone_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(drop = FALSE) +                           # show every cluster tick
  labs(x = "clusters", y = "No of cells") +
  theme_nature(7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key.size = unit(0.7, "lines"))

save_nature(p, file.path(out_dir, "Fig2i_clonalOccupy_V7"), w_cm = 18, h_cm = 7)

# End of the script
