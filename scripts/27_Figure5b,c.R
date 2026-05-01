###############################################################################
# Script name: 27_Figure5b,c.R
#
# Purpose:
#   - Perform downstream single-cell analysis of liver FNA samples
#   - Add curated liver-specific metadata (patient, visit, fibrosis stage)
#   - Normalise RNA and ADT expression
#   - Perform dimensionality reduction, Harmony integration, clustering,
#     and cell type annotation
#   - Identify RNA and ADT marker genes for liver cell populations
#   - Generate PCA, UMAP, and cluster visualisations for figures
#
# Input:
#   - Integrated Seurat object after QC and merging:
#       SeuObjFNAx.rds
#
# Main steps:
#   1. Subset Seurat object to liver samples only
#   2. Construct per-sample liver metadata (baseline / follow-up pairing)
#   3. Normalise RNA and ADT data and perform scaling
#   4. Run PCA and Harmony integration
#   5. Cluster cells and compute UMAP embeddings
#   6. Identify cluster-specific RNA and ADT markers
#   7. Annotate clusters with biologically meaningful liver cell identities
#   8. Export publication-quality figures and result tables
#
# Output:
#   - PCA plot (PowerPoint):
#       FNA_LIVER_PCAplot.pptx
#   - Cluster marker tables:
#       FNA_LIVER_cluster_markers_RNA_top20.csv
#       FNA_LIVER_cluster_markers_ADT_top10.csv
#   - Final annotated Seurat object:
#       SeuObjFNA_LIVER.rds
#   - High-resolution UMAP figure:
#       LiverFNAclusters.tiff
#
#
###############################################################################

# Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(officer)
library(rvg)

# Set working directory
setwd("~/seurat")

# Load data
SeuObj <- readRDS("SeuObjFNAx.rds")

# Subset for LIVER only
SeuObj <- subset(SeuObj, subset = Tissue %in% "LIVER")

## Build the per-tissue metadata (8 LIVER samples total)
per_tissue_meta <- data.frame(
 SampleKey          = character(),
 SampleCore         = character(),
 Patient_ID         = character(),
 Tissue             = character(),
 Visit              = character(),
 Stage              = character(),
 Pair_ID            = character(),
 Counterpart_Sample = character(),
 VisitLabel         = character(),
 stringsAsFactors   = FALSE
)

add_row <- function(sample_core, tissue, visit, patient, pair, counterpart, stage) {
 data.frame(
  SampleKey          = paste0(sample_core, "-", tissue),
  SampleCore         = sample_core,
  Patient_ID         = patient,
  Tissue             = tissue,
  Visit              = visit,
  Stage              = stage,
  Pair_ID            = pair,
  Counterpart_Sample = counterpart,
  VisitLabel         = paste0(tolower(visit), " ", pair),
  stringsAsFactors   = FALSE
 )
}

# BASELINE LIVER (4 samples)
per_tissue_meta <- rbind(
 per_tissue_meta,
 add_row("10291-1", "LIVER", "Baseline", "10291-1", "10291-1", "11303", "F1"),
 add_row("10738",   "LIVER", "Baseline", "10738",   "10738",   "11570", "F1"),
 add_row("11040",   "LIVER", "Baseline", "11040",   "11040",   "11816", "F3"),
 add_row("11183",   "LIVER", "Baseline", "11183",   "11183",   "11937", "F3")
)

# FOLLOW-UP LIVER (4 samples)
followup_map <- c("10291-1"="11303", "10738"="11570", "11040"="11816", "11183"="11937")

for (baseline in names(followup_map)) {
 fu <- followup_map[baseline]
 per_tissue_meta <- rbind(
  per_tissue_meta,
  add_row(fu, "LIVER", "Followup", baseline, baseline, baseline, NA)
 )
}

rownames(per_tissue_meta) <- NULL
stopifnot(nrow(per_tissue_meta) == 8)  # 8 LIVER samples total

# Parse SampleCore from Seurat object
SeuObj$SampleCore <- gsub("^GC-WL-", "", SeuObj$Sample)
SeuObj$SampleCore <- sub("-LIVER$", "", SeuObj$SampleCore)

# All are LIVER
SeuObj$Tissue <- "LIVER"

# Build SampleKey
SeuObj$SampleKey <- paste(SeuObj$SampleCore, SeuObj$Tissue, sep = "-")

# Add metadata to Seurat object
for (col in c("Patient_ID","Visit","Stage","Pair_ID","Counterpart_Sample","VisitLabel")) {
 SeuObj[[col]] <- NA
}

m <- match(SeuObj$SampleKey, per_tissue_meta$SampleKey)

SeuObj$Patient_ID         <- per_tissue_meta$Patient_ID[m]
SeuObj$Visit              <- per_tissue_meta$Visit[m]
SeuObj$Stage              <- per_tissue_meta$Stage[m]
SeuObj$Pair_ID            <- per_tissue_meta$Pair_ID[m]
SeuObj$Counterpart_Sample <- per_tissue_meta$Counterpart_Sample[m]
SeuObj$VisitLabel         <- per_tissue_meta$VisitLabel[m]

# CHECK
cat("\nCheck Tissue counts:\n")
print(table(SeuObj$Tissue, useNA="ifany"))

cat("\nCheck NA in VisitLabel:\n")
print(table(is.na(SeuObj$VisitLabel)))

stopifnot(all(!is.na(SeuObj$VisitLabel)))

# Print final metadata (8 LIVER samples)
print(per_tissue_meta[, c("SampleKey","SampleCore","Patient_ID","Visit","Stage","VisitLabel")], row.names = FALSE)

# ADT expression normalisation
DefaultAssay(SeuObj) <- "ADTonly"
SeuObj <- NormalizeData(SeuObj, normalization.method = "CLR", margin = 2)

# Gene expression normalisation
DefaultAssay(SeuObj) <- "RNA"
SeuObj <- NormalizeData(SeuObj)
SeuObj <- FindVariableFeatures(SeuObj)

# Scaling
all.genes <- rownames(SeuObj)
SeuObj <- ScaleData(SeuObj, features = all.genes)

# PCA
SeuObj <- RunPCA(SeuObj)

#### Fig 5.b ####
pca_df <- as.data.frame(Embeddings(SeuObj, reduction = "pca"))
pca_df$VisitLabel <- SeuObj$VisitLabel

pca_mean_df <- pca_df %>%
 group_by(VisitLabel) %>%
 summarise(
  PC_1 = mean(PC_1),
  PC_2 = mean(PC_2),
  .groups = "drop"
 )

legend_order <- c(
 "baseline 10291-1","followup 10291-1",
 "baseline 10738","followup 10738",
 "baseline 11040","followup 11040",
 "baseline 11183","followup 11183"
)

present_levels <- intersect(legend_order, unique(pca_mean_df$VisitLabel))
pca_mean_df$VisitLabel <- factor(pca_mean_df$VisitLabel, levels = present_levels)

visitlabel_colors <- c(
 "baseline 10291-1" = "#1f77b4",  "followup 10291-1" = "#8fbce5",
 "baseline 10738"   = "#2ca02c",  "followup 10738"   = "#97d397",
 "baseline 11040"   = "#ff7f0e",  "followup 11040"   = "#ffbe7a",
 "baseline 11183"   = "#9467bd",  "followup 11183"   = "#c1a3de"
)

ggp <- ggplot(
 pca_mean_df,
 aes(x = PC_1, y = PC_2, color = VisitLabel)
) +
 geom_point(size = 6) +
 geom_text_repel(
  aes(label = VisitLabel),
  size = 4,
  max.overlaps = Inf
 ) +
 scale_color_manual(
  values = visitlabel_colors,
  breaks = present_levels,
  name = "Patient / Visit"
 ) +
 theme_classic(base_size = 14) +
 labs(
  title = "PCA of LIVER samples by patient and visit",
  x = "PC1",
  y = "PC2"
 )

print(ggp)

# Save to PowerPoint

ppt <- read_pptx()
ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme")
ppt <- ph_with(ppt, dml(ggobj = ggp), location = ph_location_fullsize())
print(ppt, target = "/data/home/hdx044/plots/screpertoire/FNA_LIVER_PCAplot.pptx")

# Elbow plot
ElbowPlot(SeuObj, ndims = 50)

# Harmony integration
library(Seurat)
SeuObj <- IntegrateLayers(
 object = SeuObj, method = HarmonyIntegration,
 orig.reduction = "pca", new.reduction = "harmony",
 verbose = FALSE
)

SeuObj[["RNA"]] <- JoinLayers(SeuObj[["RNA"]])

# Clustering
SeuObj <- FindNeighbors(SeuObj, reduction = "harmony", dims = 1:30)
SeuObj <- FindClusters(SeuObj, resolution = 0.5)
SeuObj <- RunUMAP(SeuObj, dims = 1:30, reduction = "harmony")

# Find markers - RNA
DefaultAssay(SeuObj) <- "RNA"
SeuObj@misc$markers <- FindAllMarkers(SeuObj, min.pct = 0.25, logfc.threshold = 0.25)

cat("\nRNA markers per cluster:\n")
print(table(SeuObj@misc$markers$cluster))

top20_markers <- as.data.frame(SeuObj@misc$markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
write.csv(top20_markers, file="/data/home/hdx044/files/seurat/allTissue/FNA/FNA_LIVER_cluster_markers_RNA_top20.csv")

# Find markers - ADT
DefaultAssay(SeuObj) <- "ADTonly"
SeuObj@misc$markersADT <- FindAllMarkers(SeuObj, min.pct = 0.25, logfc.threshold = 0.25)

cat("\nADT markers per cluster:\n")
print(table(SeuObj@misc$markersADT$cluster))

top10_markersADT <- as.data.frame(SeuObj@misc$markersADT %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
write.csv(top10_markersADT, file="/data/home/hdx044/files/seurat/allTissue/FNA/FNA_LIVER_cluster_markers_ADT_top10.csv")

# Save updated object
saveRDS(SeuObj, "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_LIVER.rds")


#### Fig5.c ####

# Ensure the seurat_clusters column is present and numeric/integer
stopifnot(!is.null(SeuObj$seurat_clusters))

# Make a new metadata column 'cluster' with values like C0, C1, C2...
SeuObj$cluster <- paste0("C", as.integer(as.character(SeuObj$seurat_clusters)))

# Quick check
table(SeuObj$seurat_clusters, SeuObj$cluster)

#### Add cluster ID ####
new.cluster.ids <- c(
 "Naïve-Activated T cells",     # 0
 "CD56 high NK cells",          # 1
 "CD4 T cells",                 # 2
 "CD56 dim NK cells",           # 3
 "CD8 Effector Memory",         # 4
 "CD8 Cytotoxic T cells",       # 5
 "CD14 Monocytes",              # 6
 "MAIT cells",                  # 7
 "Erythroid cells",             # 8
 "Naïve B cells",               # 9
 "CD4 memory-like T cells",     #10
 "Neutrophils",                 #11
 "Non‑classical Monocytes",     #12
 "Liver ductal cells",          #13
 "Hepatocytes",                 #14
 "Blood vessels",               #15
 "Mast cells",                  #16
 "Myeloid Dendritic cells",     #17
 "cDC1",                        #18
 "Macrophage",                  #19
 "TREM2 Dendritic cells",       #20
 "Plasma cells",                #21
 "Plasmacytoid Dendritic cells",#22
 "Mast/eosinophil-like hybrid", #23
 "Liver ductal cells",          #24
 "Cycling/Proliferating cells", #25
 "TAGLN Endothelial cells"      #26
)

# Assign new cluster labels
# Create labels like "CD4 T cells (C2)"
new.cluster.ids.suffixed <- paste0(new.cluster.ids, " (C", 0:26, ")")


# Map seurat_clusters 0..26 → "Name (C#)"
SeuObj$ClusterName <- plyr::mapvalues(
 x   = SeuObj$seurat_clusters,
 from= 0:26,
 to  = new.cluster.ids.suffixed
)

sorted_cluster_ids <- c(
 # T cells
 0, 2, 10, 4, 5, 7,
 # NK
 1, 3,
 # B lineage
 9, 21,
 # Myeloid
 6, 12, 17, 18, 22, 19, 20, 11,
 # Stromal
 16, 23, 13, 24,
 # Endothelial
 15, 26,
 # Parenchymal
 14,
 # Other
 8, 25
)

# Generate the ordered labels in "Name (C#)" form:
sorted_clusters_suffixed <- new.cluster.ids.suffixed[sorted_cluster_ids + 1]

names(sorted_colors) <- sorted_clusters_suffixed

sorted_colors <- c(
 # T cells (reds)
 "#a50f15",
 "#de2d26",
 "#fb6a4a",
 "#fc9272",
 "#fcbba1",
 "#fee0d2",
 
 # NK cells (purples)
 "#6a51a3",
 "#9e9ac8",
 
 # B lineage (blues)
 "#08519c",
 "#3182bd",
 
 # Myeloid lineage (greens)
 "#006d2c",
 "#31a354",
 "#74c476",
 "#a1d99b",
 "#c7e9c0",
 "#e5f5e0",
 "#edf8e9",
 "#c6dbef",
 
 # Stromal (brown/tan)
 "#8c510a",
 "#d8b365",
 "#f6e8c3",
 "#c7eae5",
 
 # Endothelial (teals)
 "#01665e",
 "#5ab4ac",
 
 # Parenchymal (hepatocytes)
 "#d95f0e",
 
 # Other (erythroid + cycling)
 "#636363",
 "#cccccc"
)

# Make sure the factor follows your sorted order (with suffixes)
SeuObj$ClusterName <- factor(SeuObj$ClusterName, levels = sorted_clusters_suffixed)

ggp <- DimPlot(SeuObj, reduction = "umap", raster = FALSE, label = F, group.by = "ClusterName",
               cols = sorted_colors) +
 ggtitle("Liver FNA harmony integrated clusters") +
 theme_classic(base_size = 12) +
 theme(
  axis.text.x = element_text(angle = 90, hjust = 1),
  axis.line = element_line(color = "black"),
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
setwd("/data/home/hdx044/plots/seurat/liver/FNA")
tiff("LiverFNAclusters.tiff", width=22, height=15, units="cm", res=600, compression="lzw")
print(ggp)
dev.off()

# End of the script
