
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