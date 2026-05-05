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

#Load seurat object
SeuObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

# Make cluster label
SeuObj$cluster <- paste0("C", as.character(SeuObj$seurat_clusters))

Idents(SeuObj) <- "cluster"

# Arrange the clusters in increasing order
Idents(SeuObj) <- "cluster"
SeuObj@active.ident <- factor (SeuObj@active.ident,
                                 levels = c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22','C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30'))


# Map from raw cluster ID -> cell type name
cluster_to_type <- c(
 "1"  = "Naïve-Activated T cells",
 "0"  = "CD8 Effector Memory",
 "4"  = "CD8 Cytotoxic T cells",
 "7"  = "CD4 T cells",
 "21" = "CD4 Proliferating T cells",
 "8"  = "Naïve B cells",
 "27" = "B cells",
 "19" = "Plasma cells_1",
 "29" = "Plasma cells_2",
 "30" = "Plasma cells_3",
 "3"  = "CD14 Monocytes",
 "11" = "CD16 Monocytes",
 "5"  = "Macrophage",
 "10" = "Myeloid Dendritic cells",
 "14" = "Neutrophil",
 "22" = "Plasmacytoid dendritic cells",
 "23" = "TREM2 Dendritic cells",
 "2"  = "CD56 dim NK cells",
 "12" = "CD56 high NK cells",
 "13" = "Fibroblasts",
 "24" = "Mesothelial cells",
 "25" = "Adipocytes_1",
 "26" = "Adipocytes_2",
 "9"  = "ACKR1 Endothelial cells",
 "16" = "TAGLN Endothelial cells",
 "6"  = "Blood vessels",
 "20" = "Hepatocytes",
 "18" = "Liver Ductal cells",
 "17" = "Erythroid cells",
 "15" = "Mast cells",
 "28" = "Platelets"
)

# Create the cell_type column from the mapping
raw_clusters <- as.character(SeuObj$seurat_clusters)
SeuObj$cell_type <- unname(cluster_to_type[raw_clusters])

# (Optional but recommended) flag any clusters not in the map
missing_idx <- is.na(SeuObj$cell_type)
if (any(missing_idx)) {
 warning("Some clusters are not in the mapping. They will be labeled as 'Unknown'.")
 SeuObj$cell_type[missing_idx] <- paste0("Unknown_", raw_clusters[missing_idx])
}

# Create cell_type_with_cluster
SeuObj$cell_type_with_cluster <- paste0(SeuObj$cell_type, " (", SeuObj$cluster, ")")


# Plot UMAP
Idents(SeuObj) <- "cell_type_with_cluster"

# define a biologically meaningful order
SeuObj$cell_type_with_cluster <- factor(
 SeuObj$cell_type_with_cluster,
 levels = c(
  # T cells
  "Naïve-Activated T cells (C1)",
  "CD8 Effector Memory (C0)",
  "CD8 Cytotoxic T cells (C4)",
  "CD4 T cells (C7)",
  "CD4 Proliferating T cells (C21)",
  # B cells / Plasma
  "Naïve B cells (C8)",
  "B cells (C27)",
  "Plasma cells_1 (C19)",
  "Plasma cells_2 (C29)",
  "Plasma cells_3 (C30)",
  # Myeloid
  "CD14 Monocytes (C3)",
  "CD16 Monocytes (C11)",
  "Macrophage (C5)",
  "Myeloid Dendritic cells (C10)",
  "Plasmacytoid dendritic cells (C22)",
  "TREM2 Dendritic cells (C23)",
  "Neutrophil (C14)",
  # NK
  "CD56 dim NK cells (C2)",
  "CD56 high NK cells (C12)",
  # Stromal
  "Fibroblasts (C13)",
  "Mesothelial cells (C24)",
  "Adipocytes_1 (C25)",
  "Adipocytes_2 (C26)",
  # Endothelial
  "ACKR1 Endothelial cells (C9)",
  "TAGLN Endothelial cells (C16)",
  "Blood vessels (C6)",
  # Liver parenchyma
  "Hepatocytes (C20)",
  "Liver Ductal cells (C18)",
  # Other
  "Erythroid cells (C17)",
  "Mast cells (C15)",
  "Platelets (C28)"
 )
)

# Define stromal cell types
stromal_types <- c("Fibroblasts (C13)", "Mesothelial cells (C24)", "Adipocytes_1 (C25)", "Adipocytes_2 (C26)")

# Count total cells and stromal cells
cell_counts <- table(SeuObj$cell_type_with_cluster)
stromal_count <- sum(cell_counts[names(cell_counts) %in% stromal_types])

# Calculate percentage
stromal_percentage <- stromal_count / sum(cell_counts) * 100

# Print results
cat("Stromal cells:", stromal_count, "of", sum(cell_counts), 
    "(", round(stromal_percentage, 2), "% )\n")


# Define biologically grouped color palette
my_colors <- c(
 # T cells (blue/red tones)
 "#a50f15", "#de2d26", "#fb6a4a", "#fc9272", "#fcbba1",
 
 # B / Plasma (blue tones)
 "#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef",
 
 # Myeloid (greens)
 "#006d2c", "#31a354", "#74c476", "#a1d99b", "#c7e9c0",
 "#e5f5e0", "#edf8e9",
 
 # NK (purple)
 "#6a51a3", "#9e9ac8",
 
 # Stromal (brown/tan)
 "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5",
 
 # Endothelial (teal)
 "#01665e", "#5ab4ac", "#35978f",
 
 # Parenchymal (orange)
 "#d95f0e", "#fdd0a2",
 
 # Other (grey)
 "#636363", "#969696", "#cccccc"
)

#### Fig1.b ####
# Ploting UMAP for cell_type_with_cluster

ggp <- DimPlot(SeuObj, reduction = "umap", raster = FALSE, label = F, group.by = "cell_type_with_cluster",
               cols = my_colors) +
 ggtitle("All tissue harmony integrated clusters") +
 theme_classic(base_size = 12) +
 theme(
  axis.text.x = element_text(angle = 90, hjust = 1),
  axis.line = element_line(color = "black"),
  plot.title = element_blank(),
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.key.size = unit(0.8, "lines"),  # reduce point size in legend
  legend.spacing.y = unit(0.8, "cm")     # reduce spacing between items
 ) +
 guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))  # one column legend, small dots

ggp

# Final TIFF for journal
setwd("/data/home/hdx044/plots/seurat/allTissue")
tiff("allTissue_cell_type.tiff", width=22, height=15, units="cm", res=600, compression="lzw")
print(ggp)
dev.off()

setwd("/data/home/hdx044/plots/seurat/allTissue")
svg("allTissue_cell_type1.svg", width = 22, height = 15)   # width/height in cm by default for svg()
print(ggp)
dev.off()

#### Fig1.c ####
#count plots Tissue wise

Idents(SeuObj) <- "cluster"

# Arrange the cluster according to no
SeuObj@active.ident <- factor (SeuObj@active.ident,
                               levels = c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30'))

# Create tissue summary data frame
summary <- SeuObj@meta.data %>%
  group_by(cluster) %>%
  count(Tissue) %>%
  mutate(proportion = n / sum(n))

# Reorder the cluster factor levels in the summary table to match the Seurat object
summary$cluster <- factor(summary$cluster, levels = levels(Idents(SeuObj)))

# Plot summary stacked bar chart
ggp <- ggplot(summary, aes(x = cluster, y = proportion, fill = Tissue)) +
  geom_bar(stat = "identity", position = "stack") +
  guides(fill=guide_legend(title="Tissue")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_manual(values = friendly_pal("zesty_four"))+
  scale_y_continuous(expand = c(0,0))

ggp <- ggp+ 
  labs(title = "Proportions of cells in each cluster by tissue")+
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.8, "lines"),  # reduce point size in legend
    legend.spacing.y = unit(1, "cm")     # reduce spacing between items
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))  # one column legend, small dots

ggp

# Save high-resolution SVG
ggsave(
  filename = "allTissue_Proportion_by_tissue.svg",
  plot = ggp,
  device = "svg",
  path = "/data/home/hdx044/plots/seurat/allTissue",
  width = 7,     
  height = 3,
  units = "in",
  dpi = 600
)

#### Fig1.d ####

# Extract metadata including Tissue
df <- SeuObj@meta.data %>%
  dplyr::select(cell_type_with_cluster, Stage, Tissue)

# Summarise proportions per Stage and Tissue
df_summary <- df %>%
  group_by(Tissue, Stage, cell_type_with_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Tissue, Stage) %>%
  mutate(prop = count / sum(count))

# Ensure order of factors
df_summary$Stage <- factor(df_summary$Stage, 
                           levels = c("Healthy", "No_fibrosis", "Fibrosis"))
df_summary$cell_type_with_cluster <- factor(
  df_summary$cell_type_with_cluster,
  levels = levels(SeuObj$cell_type_with_cluster)
)

write.csv(df_summary, "/data/home/hdx044/files/seurat/allTissue/all_cluster_proportion_cluster_tissue_stage.csv")

# Palette must match levels
names(my_colors) <- levels(df_summary$cell_type_with_cluster)

# Plot, faceting by Tissue
ggp <- ggplot(df_summary, aes(x = Stage, y = prop, fill = cell_type_with_cluster)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = my_colors, guide = guide_legend(reverse = FALSE)) +
  facet_wrap(~ Tissue, ncol = 4) +   # split plots by tissue
  labs(
    title = "Cell type proportions across tissues and stages",
    x = NULL, y = "Proportion of cells", fill = "Cell types"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size = 12, color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12)
  )+
  guides(fill = guide_legend(ncol = 1))


ggp

# Final TIFF for journal
setwd("/data/home/hdx044/plots/seurat/allTissue")
tiff("allTissue_cell_type_Proportion_by_tissue_stage.tiff", width=25, height=25, units="cm", res=600, compression="lzw")
print(ggp)
dev.off()

# End of the script
