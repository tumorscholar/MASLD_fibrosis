##################################################################################
# Script name: 18_MakingcellchatObj.R
#
# Generation of Seurat object with shared expanded T cells
# for CellChat and downstream cell–cell communication analyses
#
# This script:
#  - Loads a full Seurat object containing all cell populations
#  - Identifies expanded T-cell clonotypes shared across multiple
#    solid tissues using a curated expanded T-cell reference object
#  - Annotates shared expanded T cells in the full dataset
#  - Collapses shared expanded T cells into a unified "Tcell" group
#  - Removes selected non-T-cell clusters not required for analysis
#  - Reorders cluster identities to prioritise T cells
#  - Saves a final Seurat object for CellChat input
#
# Input files:
#  - SeuObjx.rds: Full Seurat object with all cells
#  - expandedTcellFinal.rds: Expanded T-cell object with clone sharing annotations
#
# Output file:
#  - SharedTcells&NonTcell.rds: Annotated Seurat object containing
#    shared expanded T cells and remaining non-T-cell populations
#
# Notes:
#  - Shared expanded T cells are defined as clonotypes detected
#    across more than one solid tissue compartment
#  - Metadata fields 'cellchat' and 'cluster_cellchat' are created
#    to support CellChat-based communication inference
##################################################################################
library(Seurat)
library(dplyr)

# load Full Seurat object with all cells
SeuObj <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds")

# Expanded T cell object with shared clone annotation
expanded_obj <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinal.rds")

# Identify barcodes of shared expanded clonotypes
shared_barcodes <- colnames(expanded_obj)[
  expanded_obj$expanded_shared_class == "c_Shared_between_solid_tissues"
]

# Initialise new metadata column
SeuObj$cellchat <- NA

# Label shared expanded T cells
SeuObj$cellchat[colnames(SeuObj) %in% shared_barcodes] <- "T cell"

# Create cluster label for downstream analyses
# Ensure cluster is character
SeuObj$cluster <- as.character(SeuObj$cluster)

SeuObj$cluster_cellchat <- SeuObj$cluster
SeuObj$cluster_cellchat[SeuObj$cellchat == "T cell"] <- "Tcell"

# Remove unwanted non-T clusters
clusters_to_remove <- c("C0", "C1", "C4", "C7", "C21")
SeuObj$cluster_cellchat[SeuObj$cluster_cellchat %in% clusters_to_remove] <- NA

# Subset to keep shared T cells and remaining clusters
SeuObj <- subset(
  SeuObj,
  subset = !is.na(cluster_cellchat) | cluster_cellchat == "Tcell"
)

# Set factor order (T cells first)
clusters <- as.character(SeuObj$cluster_cellchat)
other_clusters <- setdiff(clusters, "Tcell")
cluster_numbers <- as.numeric(sub("C", "", other_clusters))
ordered_clusters <- paste0("C", sort(cluster_numbers))

SeuObj$cluster_cellchat <- factor(
  SeuObj$cluster_cellchat,
  levels = c("Tcell", ordered_clusters)
)

# Save final object
saveRDS(
  SeuObj,
  file = "/data/Blizard-AlazawiLab/rk/seurat/SharedTcells&NonTcell.rds"
)

# End of the script