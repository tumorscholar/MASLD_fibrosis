###############################################################################
# Script name: 03_integration_clustering_and_fibrosis_de_analysis.R
#
# Description:
# This script performs downstream analysis of quality‑controlled single‑cell
# RNA‑seq data, including metadata augmentation, data integration, clustering,
# and differential expression analysis across fibrosis stages.
#
# The input to this script is a QC‑filtered and ambient RNA‑corrected Seurat
# object generated from previous preprocessing steps. This script focuses on
# integrating samples, identifying cellular clusters, and identifying cluster‑
# specific RNA and ADT markers that are subsequently used for cluster
# annotation.
#
# Key steps performed in this script include:
#  - Loading a QC‑filtered Seurat object
#  - Adding clinical and phenotypic metadata, including diabetic status and
#    fibrosis stage annotations at multiple levels (binary, grouped, and
#    detailed staging)
#  - Normalization and scaling of RNA and ADT (protein) data
#  - Dimensionality reduction using PCA
#  - Batch correction and integration using Harmony
#  - Graph‑based clustering
#
# Marker identification:
#  - Identification of cluster‑specific marker genes (RNA) using FindAllMarkers
#  - Identification of cluster‑specific protein markers (ADT)
#  - Export of RNA and ADT marker tables for downstream use
#
# Note:
# This script does NOT perform cell‑type or cluster annotation directly.
# Identified RNA and ADT markers are intended to be used in subsequent manual
# or supervised annotation steps.
#
# Differential expression analysis:
#  - Comparison of Fibrosis versus No‑fibrosis conditions within each cluster
#  - Separate differential expression analyses for RNA and ADT data
#
# Outputs:
#  - Integrated and clustered Seurat object
#  - Cluster‑level RNA and ADT marker tables
#  - Fibrosis versus No‑fibrosis differential expression results per cluster
#
# Purpose:
# To integrate multi‑tissue scRNA‑seq data, identify cluster‑specific molecular
# signatures, and generate marker sets for downstream cluster annotation and
# fibrosis‑focused analyses.
#

###############################################################################

# Load packages
library(Seurat)
library(tidyverse)

# Set working directory
setwd("~/seurat")

# Load data
SeuObj <- readRDS("SeuObjx.rds")

# Add stage information
#Create lookup table with stage and sample
lookup <- data.frame(Patient_ID = c('11183','11051','10113-1','10202','10380','10634','10738','9680','10205','9999','10113-2','9961','10742',  '10203','10291-2','11327','9932', '9991','11040','11471'),
                     Diabetic = c('No','No','No','No','No','No','No','No','No','No','No','No','No','Yes', 'Yes', 'Yes','Yes','Yes', 'Yes','Yes'))

# Extract meta data
meta <- SeuObj@meta.data

# Create new column in meta data for Stage
meta$Diabetic <- NA

# Run for all samples in a for loop to populate the whole column
for (x in 1:length(lookup$Patient_ID)) {
 meta$Diabetic[grep(lookup$Patient_ID[x], meta$Patient_ID)] <- lookup$Diabetic[x]
}

# Add meta data back to seurat object
SeuObj <- AddMetaData(SeuObj, metadata = meta)

# Add stage details
# Create lookup table with stage and sample
lookup <- data.frame(Patient_ID = c('11471', '11183', '11040','11051','10113-1','10202','10203','10291-2','10380','10634','10738','9680','9932', '9991','10205','9999','10113-2','9961','10742', '11327'),
                     StageSep = c('F3','F3','F3','F1','F1','F1','F1','F1','F1','F1','F1','F1','F2','F2', 'F2','F0','F0','F0','H', 'F1'))

# Extract meta data
meta <- SeuObj@meta.data

# Create new column in meta data for Stage
meta$StageSep <- NA

# Run for all samples in a for loop to populate the whole column
for (x in 1:length(lookup$Patient_ID)) {
 meta$StageSep[grep(lookup$Patient_ID[x], meta$Patient_ID)] <- lookup$StageSep[x]
}

# Add meta data back to seurat object
SeuObj <- AddMetaData(SeuObj, metadata = meta)

# Add stage information
#Create lookup table with stage and sample
lookup <- data.frame(Patient_ID = c('11471', '11327', '11183', '11040','11051','10113-1','10202','10203','10291-2','10380','10634','10738','9680','9932', '9991','10205','9999','10113-2','9961','10742'),
                     Stage = c('Fibrosis', 'Fibrosis','Fibrosis', 'Fibrosis','Fibrosis','Fibrosis','Fibrosis','Fibrosis','Fibrosis','Fibrosis','Fibrosis','Fibrosis','Fibrosis','Fibrosis','Fibrosis', 'Fibrosis', 'No_fibrosis','No_fibrosis','No_fibrosis', 'Healthy'))

# Extract meta data
meta <- SeuObj@meta.data

# Create new column in meta data for Stage
meta$Stage <- NA

# Run for all samples in a for loop to populate the whole column
for (x in 1:length(lookup$Patient_ID)) {
 meta$Stage[grep(lookup$Patient_ID[x], meta$Patient_ID)] <- lookup$Stage[x]
}

# Add meta data back to seurat object
SeuObj <- AddMetaData(SeuObj, metadata = meta)

#### ADT expression normalisation ####
DefaultAssay(SeuObj) <- "ADTonly"
SeuObj <- NormalizeData(SeuObj, normalization.method = "CLR", margin = 2)

#### Gene expression normalisation ####
DefaultAssay(SeuObj) <- "RNA"

# Normalisation
SeuObj <- NormalizeData(SeuObj)

# FindVariableFeatures
SeuObj <- FindVariableFeatures(SeuObj)

# Scaling
all.genes <- rownames(SeuObj)
SeuObj <- ScaleData(SeuObj, features = all.genes)

#PCA
SeuObj <- RunPCA(SeuObj)

# Determine dimensionality of the data
ElbowPlot(SeuObj, ndims = 50)

# Harmony integration
SeuObj <- IntegrateLayers(
 object = SeuObj, method = HarmonyIntegration,
 orig.reduction = "pca", new.reduction = "harmony",
 verbose = FALSE
)

# Save obj in case R studio crashes
saveRDS(SeuObj, '/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

# Re-join layers after integration
SeuObj[["RNA"]] <- JoinLayers(SeuObj[["RNA"]])

# Checking clusters of harmony integration
SeuObj <- FindNeighbors(SeuObj, reduction = "harmony", dims = 1:30)
SeuObj <- FindClusters(SeuObj, resolution = 0.5)
set.seed(1234)
SeuObj <- RunUMAP(SeuObj, dims = 1:30, reduction = "harmony")

# Identify markers of each clusters
SeuObj@misc$markers <- FindAllMarkers(SeuObj, min.pct = 0.25, logfc.threshold = 0.25)
table(SeuObj@misc$markers$cluster)
top20_markers <- as.data.frame(SeuObj@misc$markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
top20_markers
write.csv(top20_markers, file="/data/home/hdx044/files/seurat/allTissue/cluster_markers_RNA_top20_allTissue.csv")

#### sig genes in excel ####

library(openxlsx)

# Ensure identities are seurat_clusters
Idents(SeuObj) <- "seurat_clusters"

# Filter significant genes
sig_markers <- SeuObj@misc$markers %>%
 filter(p_val_adj < 0.05)

# Sort within each seurat cluster by adjusted p-value and log2FC
sig_markers <- sig_markers %>%
 group_by(cluster) %>%
 arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE)

## Make sure cluster is numeric
sig_markers$cluster <- as.numeric(as.character(sig_markers$cluster))

# Create workbook
wb <- createWorkbook()
for (cl in sort(unique(sig_markers$cluster))) {
 addWorksheet(wb, paste0("Cluster_", cl))
 cluster_data <- sig_markers %>% filter(cluster == cl)
 writeData(wb, paste0("Cluster_", cl), cluster_data)
}
saveWorkbook(
 wb,
 "/data/home/hdx044/files/seurat/allTissue/Significant_Gene_Markers_by_seurat_clusters.xlsx",
 overwrite = TRUE
)

DefaultAssay(SeuObj) <- 'ADTonly'
SeuObj@misc$markersADT <- FindAllMarkers(SeuObj, min.pct = 0.25, logfc.threshold = 0.25)
table(SeuObj@misc$markersADT$cluster)
top10_markersADT <- as.data.frame(SeuObj@misc$markersADT %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
top10_markersADT
write.csv(top10_markersADT, file="/data/home/hdx044/files/seurat/allTissue/cluster_markers_ADT_top10_allTissue.csv")

# Ensure identities are seurat_clusters
Idents(SeuObj) <- "seurat_clusters"

# Filter significant genes
sig_markers <- SeuObj@misc$markersADT %>%
 filter(p_val_adj < 0.05)

# Sort within each seurat cluster by adjusted p-value and log2FC
sig_markers <- sig_markers %>%
 group_by(cluster) %>%
 arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE)

## Make sure cluster is numeric
sig_markers$cluster <- as.numeric(as.character(sig_markers$cluster))

# Create workbook
wb <- createWorkbook()
for (cl in sort(unique(sig_markers$cluster))) {
 addWorksheet(wb, paste0("Cluster_", cl))
 cluster_data <- sig_markers %>% filter(cluster == cl)
 writeData(wb, paste0("Cluster_", cl), cluster_data)
}
saveWorkbook(
 wb,
 "/data/home/hdx044/files/seurat/allTissue/Significant_ADT_Markers_by_seurat_clusters.xlsx",
 overwrite = TRUE
)

DefaultAssay(SeuObj) <- 'RNA'

# Add cluster.stage in metadata
SeuObj$Cluster.StageSep <- paste(SeuObj$seurat_clusters , SeuObj$StageSep, sep = "_")

# Save obj in case R studio crashes

saveRDS(SeuObj, '/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

#### Fibrosis gene markers in each cluster ####

Idents(SeuObj) <- "seurat_clusters"
DefaultAssay(SeuObj) <- 'RNA'

SeuObj$Stage <- factor(
 SeuObj$Stage,
 levels = c("No_fibrosis", "Fibrosis")
)

library(dplyr)
library(openxlsx)

clusters <- levels(Idents(SeuObj))

de_list <- list()

for (cl in clusters) {
 
 message("Processing cluster: ", cl)
 
 cells_cl <- WhichCells(
  SeuObj,
  idents = cl,
  expression = Stage %in% c("Fibrosis", "No_fibrosis")
 )
 
 if (length(cells_cl) < 30) next  # safeguard
 
 sub_obj <- subset(SeuObj, cells = cells_cl)
 
 Idents(sub_obj) <- "Stage"
 
 de <- FindMarkers(
  sub_obj,
  ident.1 = "Fibrosis",
  ident.2 = "No_fibrosis",
  logfc.threshold = 0.25,
  min.pct = 0.25,
  test.use = "wilcox"
 )
 
 if (nrow(de) == 0) next
 
 de <- de %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj)
 
 de_list[[paste0("Cluster_", cl)]] <- de
}

wb <- createWorkbook()

for (nm in names(de_list)) {
 addWorksheet(wb, nm)
 writeData(wb, nm, de_list[[nm]])
}

saveWorkbook(
 wb,
 file = "/data/home/hdx044/files/seurat/allTissue/Fibrosis_vs_NoFibrosis_DE_by_cluster.xlsx",
 overwrite = TRUE
)

#### Fibrosis ADT markers in each cluster ####
DefaultAssay(SeuObj) <- 'ADTonly'

SeuObj$Stage <- factor(
 SeuObj$Stage,
 levels = c("No_fibrosis", "Fibrosis")
)

library(dplyr)
library(openxlsx)

clusters <- levels(Idents(SeuObj))

de_list <- list()

for (cl in clusters) {
 
 message("Processing cluster: ", cl)
 
 cells_cl <- WhichCells(
  SeuObj,
  idents = cl,
  expression = Stage %in% c("Fibrosis", "No_fibrosis")
 )
 
 if (length(cells_cl) < 30) next  # safeguard
 
 sub_obj <- subset(SeuObj, cells = cells_cl)
 
 Idents(sub_obj) <- "Stage"
 
 de <- FindMarkers(
  sub_obj,
  ident.1 = "Fibrosis",
  ident.2 = "No_fibrosis",
  logfc.threshold = 0.25,
  min.pct = 0.25,
  test.use = "wilcox"
 )
 
 if (nrow(de) == 0) next
 
 de <- de %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj)
 
 de_list[[paste0("Cluster_", cl)]] <- de
}

wb <- createWorkbook()

for (nm in names(de_list)) {
 addWorksheet(wb, nm)
 writeData(wb, nm, de_list[[nm]])
}

saveWorkbook(
 wb,
 file = "/data/home/hdx044/files/seurat/allTissue/Fibrosis_vs_NoFibrosis_DE_ADT_by_cluster.xlsx",
 overwrite = TRUE
)


# End of the script
