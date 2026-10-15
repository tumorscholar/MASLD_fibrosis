###############################################################################
# Script name: 04_Figure1.R
#
# Description:
# This script performs cluster labeling, cell‑type assignment, and visualisation
# of integrated single‑cell RNA‑seq data across all tissues. The input to this
# script is an integrated and clustered Seurat object generated from previous
# preprocessing, integration, and clustering steps.
#
# The primary purpose of this script is to map Seurat cluster IDs to biologically
# meaningful cell‑type labels based on previously identified RNA and ADT marker
# genes, and to generate plots summarising cellular composition across tissues 
# and disease stages.
#
# Key steps performed in this script include:
#  - Loading an integrated Seurat object containing Harmony‑corrected clusters
#  - Creating ordered cluster labels for consistent visualisation
#  - Mapping cluster IDs to cell‑type names using predefined marker‑based
#    annotations
#  - Generating combined cell‑type and cluster labels for plotting
#  - Calculating cell‑type proportions across tissues and fibrosis stages
#
# Visualisations generated in this script include:
#  - UMAP plots coloured by annotated cell types
#  - Stacked bar plots showing tissue contributions per cluster
#  - Proportion plots of cell types across tissues and disease stages
#
# Outputs:
#  - Annotated Seurat object with cell‑type metadata
#  - Figures (UMAPs, bar plots, dot plots)
#  - CSV files summarising cell‑type proportions
#
# Purpose:
# To assign biologically interpretable cell‑type labels to clusters and to
# visualise cellular composition and marker expression across tissues in MASLD.
#

###############################################################################

# Load packages
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
library(tibble)
library(tidyr)
library(ggpubfigs)
library(SeuratExtend)
library(colorspace)
library(ggrepel)

source("00_nature_style.R")

# PUBLICATION THEME — define FIRST before any plots
my_theme <- theme_nature(7) +
  theme(
    axis.text.x     = element_text(size = 7, color = "black"),
    axis.text.y     = element_text(size = 7, color = "black"),
    axis.title.x    = element_text(size = 8, color = "black", face = "bold"),
    axis.title.y    = element_text(size = 8, color = "black", face = "bold"),
    axis.line       = element_line(color = "black", linewidth = 0.4),
    axis.ticks      = element_line(color = "black", linewidth = 0.4),
    panel.background= element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.text      = element_text(size = 7, face = "bold", color = "black"),
    strip.background= element_rect(fill = "white", color = "black", linewidth = 0.4),
    plot.title      = element_text(size = 8, face = "plain", color = "black", hjust = 0),
    legend.position    = "right",
    legend.direction   = "vertical",
    legend.background  = element_rect(fill = "white", color = NA),
    legend.title       = element_text(size = 7, face = "bold", color = "black"),
    legend.text        = element_text(size = 7, color = "black"),
    legend.key.size    = unit(0.7, "lines"),
    legend.key         = element_rect(fill = "white", color = NA),
    legend.spacing.y   = unit(0.15, "cm"),
    legend.box.margin  = margin(2, 2, 2, 2),
    plot.margin        = margin(4, 4, 4, 4, "pt")
  )

# Save figure
save_figure <- function(plot,
                        filename,
                        width_cm  = 18.3,
                        height_cm = 12,
                        path      = "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure1/") {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  cairo_pdf(file.path(path, paste0(filename, ".pdf")),
            width = width_cm/2.54, height = height_cm/2.54)
  print(plot); dev.off()
  
  svglite::svglite(file.path(path, paste0(filename, ".svg")),
                   width = width_cm/2.54, height = height_cm/2.54)
  print(plot); dev.off()
  
  # type="cairo": on a headless node the default (Xlib) backend cannot render the
  # Arial/showtext text and errors with "invalid font type"; cairo matches the PDF.
  tiff(file.path(path, paste0(filename, ".tiff")),
       width=width_cm, height=height_cm,
       units="cm", res=600, compression="lzw", type="cairo")
  print(plot); dev.off()
  
  message("Saved: ", filename)
}

# ---- Load integrated + clustered object (from script 03) --------------------
SeuObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds')

# ---- Cluster labels (C0 ... Cn), derived from the data (no hard-coded count) -
n_clusters     <- length(levels(as.factor(SeuObj$seurat_clusters)))   # 34 -> C0..C33
cluster_levels <- paste0("C", 0:(n_clusters - 1))
SeuObj$cluster <- paste0("C", as.character(SeuObj$seurat_clusters))
SeuObj$cluster <- factor(SeuObj$cluster, levels = cluster_levels)
Idents(SeuObj) <- "cluster"

# ---- Cluster -> cell-type annotation (30 clusters, V7 seeded run) ------------
cluster_to_type <- c(
  "0"  = "Naive/TCM T cells",
  "1"  = "Cytotoxic CD8 T cells",
  "2"  = "Activated CD8 T cells",
  "3"  = "Classical Monocytes",
  "4"  = "Cytotoxic NK cells",
  "5"  = "Capillary Endothelial cells",
  "6"  = "Memory T cells",
  "7"  = "SPP1+ Macrophages",
  "8"  = "Naive B cells",
  "9"  = "FOLR2+ Resident Macrophages",
  "10" = "Venous Endothelial cells",
  "11" = "Non-classical Monocytes",
  "12" = "Fibroblasts",
  "13" = "Neutrophils",
  "14" = "cDC2",
  "15" = "TAGLN Endothelial cells",
  "16" = "Erythrocytes",
  "17" = "Cholangiocytes",
  "18" = "Mast cells",
  "19" = "Plasma cells",
  "20" = "Activated T cells",
  "21" = "Hepatocytes",
  "22" = "Cycling cells",
  "23" = "Inflammatory Monocytes",
  "24" = "Mesothelial cells",
  "25" = "pDC",
  "26" = "TREM2+ Macrophages",
  "27" = "Fibroblasts",
  "28" = "Adipocytes",
  "29" = "Plasma cells"
)

raw_clusters      <- as.character(SeuObj$seurat_clusters)
SeuObj$cell_type  <- unname(cluster_to_type[raw_clusters])
missing_idx <- is.na(SeuObj$cell_type)
if (any(missing_idx)) {
  warning("Clusters not in the map -> labelled Unknown_*")
  SeuObj$cell_type[missing_idx] <- paste0("Unknown_", raw_clusters[missing_idx])
}
SeuObj$cell_type_with_cluster <- paste0(SeuObj$cell_type, " (", SeuObj$cluster, ")")

# ---- Biologically grouped order + matching palette --------------------------
# Order is by lineage; labels/palette are built from cluster_to_type so they
# can never drift out of sync with the annotation.
lineage_order_ids <- c(
  0, 2, 1, 6, 20,        # T cells
  4,                     # NK
  8, 19, 29,             # B / plasma
  3, 11, 23,             # monocytes
  7, 9, 26,              # macrophages
  14, 25,                # DC
  13, 18,                # granulocytes (neutrophil / mast)
  12, 27, 24, 28,        # stromal
  5, 10, 15,             # endothelial
  17, 21,                # epithelial / parenchymal
  16, 22                 # other (RBC / cycling)
)
stopifnot(setequal(lineage_order_ids, 0:(n_clusters - 1)))   # all clusters, once

type_levels <- paste0(cluster_to_type[as.character(lineage_order_ids)],
                      " (C", lineage_order_ids, ")")

SeuObj$cell_type_with_cluster <- factor(SeuObj$cell_type_with_cluster,
                                        levels = type_levels)

# 30 colours, grouped by lineage (same order as lineage_order_ids), named so
# DimPlot maps colour -> cell type by NAME (order-independent).
my_colors <- c(
  "#7f0000","#b30000","#d7301f","#ef6548","#fc8d59",   # T (reds)
  "#6a51a3",                                           # NK (purple)
  "#08519c","#3182bd","#6baed6",                       # B / plasma (blues)
  "#00441b","#238b45","#66c2a4",                       # monocytes (dark green)
  "#41ab5d","#74c476","#a1d99b",                       # macrophages (green)
  "#addd8e","#78c679",                                 # DC (yellow-green)
  "#fec44f","#fe9929",                                 # granulocytes (amber)
  "#8c510a","#bf812d","#dfc27d","#f6e8c3",             # stromal (brown)
  "#01665e","#35978f","#80cdc1",                       # endothelial (teal)
  "#d94801","#fd8d3c",                                 # epithelial/parenchymal (orange)
  "#737373","#bdbdbd"                                  # other (grey)
)
stopifnot(length(my_colors) == n_clusters)
names(my_colors) <- type_levels

# ---- Write the annotation back to the object + a lookup table ---------------
saveRDS(SeuObj, "/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds")
dir.create("/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure1", showWarnings = FALSE, recursive = TRUE)
write.csv(
  data.frame(cluster = paste0("C", names(cluster_to_type)),
             cell_type = unname(cluster_to_type)),
  "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure1/cluster_annotation_V7.csv",
  row.names = FALSE
)

# ---- Stromal fraction (names updated to the new annotation) -----------------
stromal_types <- c("Fibroblasts (C12)", "Mesothelial cells (C24)", "Adipocytes (C27)")
cell_counts   <- table(SeuObj$cell_type_with_cluster)
stromal_count <- sum(cell_counts[names(cell_counts) %in% stromal_types])
cat("Stromal cells:", stromal_count, "of", sum(cell_counts),
    "(", round(stromal_count / sum(cell_counts) * 100, 2), "% )\n")

#### FIG 1B - UMAP ####
Idents(SeuObj) <- "cell_type_with_cluster"
ggp_b <- DimPlot(SeuObj, reduction = "umap", raster = FALSE, label = FALSE,
                 group.by = "cell_type_with_cluster", cols = my_colors) +
  ggtitle("All tissue harmony integrated clusters") +
  my_theme +
  theme(legend.position = "none")

label_data <- ggp_b$data %>%
  group_by(cell_type_with_cluster) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
cat("Label count:", nrow(label_data), "\n")

ggp_b <- ggp_b +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(x = umap_1, y = umap_2, label = cell_type_with_cluster),
    inherit.aes = FALSE, size = 2.0, fontface = "plain",
    color = "black", bg.color = "white", bg.r = 0.12,
    box.padding = 0.4, point.padding = 0.3,
    segment.color = "grey40", segment.size = 0.25, segment.alpha = 0.6,
    max.overlaps = Inf, min.segment.length = 0.3,
    force = 3, force_pull = 0.5, seed = 42
  )
save_figure(ggp_b, "Fig1b_UMAP_labelled", width_cm = 12, height_cm = 12)

#### FIG 1C - Tissue proportion per cluster ####
# ordered cell-type labels, by cluster number C0..Cn
cl_ids <- sort(as.integer(as.character(unique(SeuObj$seurat_clusters))))
type_by_cluster <- paste0(cluster_to_type[as.character(cl_ids)], " (C", cl_ids, ")")

# Build an atomic data.frame first: a Seurat meta.data list-column (or a plyr-
# masked count) makes count() fail with "Argument 'x' is not a vector: list".
.mdc <- data.frame(cell_type_with_cluster = as.character(SeuObj@meta.data$cell_type_with_cluster),
                   Tissue = as.character(SeuObj@meta.data$Tissue), stringsAsFactors = FALSE)
summary_c <- .mdc %>%
  dplyr::count(cell_type_with_cluster, Tissue, name = "n") %>%
  group_by(cell_type_with_cluster) %>%
  mutate(proportion = n / sum(n)) %>% ungroup()
summary_c$cell_type_with_cluster <- factor(summary_c$cell_type_with_cluster,
                                           levels = type_by_cluster)

ggp_c <- ggplot(summary_c, aes(x = cell_type_with_cluster, y = proportion, fill = Tissue)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = tissue_cols) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Tissue composition per cluster", x = NULL,
       y = "Proportion of cells", fill = "Tissue") +
  my_theme +
  theme(plot.title  = element_text(size = 8, face = "plain", color = "black", hjust = 0),
        axis.text.x = element_text(size = 6, color = "black", angle = 90, hjust = 1, vjust = 0.5)) +
  guides(fill = guide_legend(ncol = 1))
# (no bare `ggp_c` auto-print: it renders to R's default pdf() device, which cannot
#  handle the Arial/showtext font on the node - "invalid font type". save_figure
#  below renders via cairo_pdf and is the only render we need.)
save_figure(ggp_c, "Fig1c_TissueProportion", width_cm = 18.3, height_cm = 9)

#### FIG 1D - Cell type proportions by tissue and fibrosis stage ####
df <- data.frame(cell_type_with_cluster = as.character(SeuObj@meta.data$cell_type_with_cluster),
                 Stage  = as.character(SeuObj@meta.data$Stage),
                 Tissue = as.character(SeuObj@meta.data$Tissue),
                 stringsAsFactors = FALSE) %>%
  filter(!is.na(Stage))
df$Stage <- dplyr::recode(as.character(df$Stage),
                          "Fibrosis" = "Fibrosis", "No_fibrosis" = "No Fibrosis", "Healthy" = "Healthy")
df$Stage <- factor(df$Stage, levels = c("Healthy", "No Fibrosis", "Fibrosis"))

df_summary <- df %>%
  group_by(Tissue, Stage, cell_type_with_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Tissue, Stage) %>% mutate(prop = count / sum(count))
df_summary$cell_type_with_cluster <- factor(df_summary$cell_type_with_cluster,
                                            levels = levels(SeuObj$cell_type_with_cluster))
write.csv(df_summary,
          "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure1/all_cluster_proportion_cluster_tissue_stage.csv")

ggp_d <- ggplot(df_summary, aes(x = Stage, y = prop, fill = cell_type_with_cluster)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = my_colors, guide = guide_legend(reverse = FALSE)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ Tissue, ncol = 4) +
  labs(title = "Cell type proportions by tissue and fibrosis stage",
       x = NULL, y = "Proportion of cells", fill = "Cell type") +
  my_theme +
  theme(plot.title       = element_text(size = 8, face = "plain", color = "black", hjust = 0),
        axis.text.x      = element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
        strip.text       = element_text(size = 7, face = "bold", color = "black"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4)) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 2)))
save_figure(ggp_d, "Fig1d_CellTypeProportion", width_cm = 14, height_cm = 14)

# End of the script