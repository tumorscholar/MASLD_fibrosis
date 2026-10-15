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
setwd("/data/Blizard-AlazawiLab/rk/seurat")

# Load data  (V7 main object from script 02)
SeuObj <- readRDS("SeuObjV7.rds")

# =============================================================================
# Clinical metadata -- single source of truth (from PatientsDetails.xlsx).
# StageSep = fine fibrosis stage (H/F0/F1/F2/F3); Stage = coarse grouping used
# in the paper (Healthy / No_fibrosis [F0] / Fibrosis [F1-F3]).
# 10742 = healthy control (not in the patient sheet). 10291-1 = new liver (F1, diabetic).
# Matching is by EXACT Patient_ID (no partial grep).
# =============================================================================
clinical <- data.frame(
  Patient_ID = c("9680","9932","9961","9991","9999","10113-1","10113-2","10202",
                 "10203","10205","10291-1","10291-2","10380","10634","10738","10742",
                 "11040","11051","11183","11327","11471"),
  StageSep   = c("F1","F2","F0","F2","F0","F1","F0","F1",
                 "F1","F2","F1","F1","F1","F1","F1","H",
                 "F3","F1","F3","F1","F3"),
  Diabetic   = c("No","Yes","No","Yes","No","No","No","Yes",
                 "No","No","Yes","Yes","No","No","No","No",
                 "Yes","No","No","Yes","Yes"),
  stringsAsFactors = FALSE
)
clinical$Stage <- ifelse(clinical$StageSep == "H", "Healthy",
                  ifelse(clinical$StageSep == "F0", "No_fibrosis", "Fibrosis"))

meta <- SeuObj@meta.data
idx  <- match(meta$Patient_ID, clinical$Patient_ID)
meta$Diabetic <- clinical$Diabetic[idx]
meta$StageSep <- clinical$StageSep[idx]
meta$Stage    <- clinical$Stage[idx]
if (any(is.na(idx))) warning("Patient_ID(s) with no clinical match: ",
                             paste(unique(meta$Patient_ID[is.na(idx)]), collapse=", "))
SeuObj <- AddMetaData(SeuObj, metadata = meta[, c("Diabetic","StageSep","Stage")])

#### ADT expression normalisation ####
# Drop isotype controls (Isotype-*) + HTO hashtag oligos (C0251/C0252/C0253)
# from the ADT panel BEFORE CLR, so the per-cell geometric mean is computed on
# the real antibody panel only (130 markers). HuMs.* cross-reactive markers are
# kept; the full 'ADT' assay is left untouched.
adt_drop <- c(grep("^Isotype", rownames(SeuObj[["ADTonly"]]), value = TRUE, ignore.case = TRUE),
              intersect(c("C0251", "C0252", "C0253"), rownames(SeuObj[["ADTonly"]])))
keep_adt <- setdiff(rownames(SeuObj[["ADTonly"]]), adt_drop)
SeuObj[["ADTonly"]] <- CreateAssayObject(
  counts = GetAssayData(SeuObj, assay = "ADTonly", layer = "counts")[keep_adt, , drop = FALSE])
cat("ADTonly: dropped", length(adt_drop), "control/HTO features ->", length(keep_adt), "antibodies\n")

DefaultAssay(SeuObj) <- "ADTonly"
SeuObj <- NormalizeData(SeuObj, normalization.method = "CLR", margin = 2)

#### Gene expression normalisation ####
DefaultAssay(SeuObj) <- "RNA"

# Normalisation
SeuObj <- NormalizeData(SeuObj)

# FindVariableFeatures
SeuObj <- FindVariableFeatures(SeuObj)

# Scaling
SeuObj <- ScaleData(SeuObj, features = VariableFeatures(SeuObj))

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
saveRDS(SeuObj, '/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds')

# Re-join layers after integration
SeuObj[["RNA"]] <- JoinLayers(SeuObj[["RNA"]])

# Checking clusters of harmony integration
SeuObj <- FindNeighbors(SeuObj, reduction = "harmony", dims = 1:30)
SeuObj <- FindClusters(SeuObj, resolution = 0.5)
set.seed(1234)
SeuObj <- RunUMAP(SeuObj, dims = 1:30, reduction = "harmony", seed.use = 1234)

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
saveRDS(SeuObj, '/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds')

#### Fibrosis gene markers in each cluster ####

Idents(SeuObj) <- "seurat_clusters"
DefaultAssay(SeuObj) <- 'RNA'

library(dplyr)
library(openxlsx)

Idents(SeuObj) <- "seurat_clusters"
DefaultAssay(SeuObj) <- "RNA"

clusters <- levels(Idents(SeuObj))
de_list <- list()

for (cl in clusters) {
  
  message("Processing cluster: ", cl)
  
  cells_cl <- WhichCells(
    SeuObj,
    idents = cl,
    expression = Stage %in% c("Fibrosis", "No_fibrosis")
  )
  
  # Skip empty/small clusters
  if (length(cells_cl) < 30) {
    message("Skipping cluster ", cl, ": <30 cells")
    next
  }
  
  sub_obj <- subset(SeuObj, cells = cells_cl)
  
  # Check group sizes
  stage_counts <- table(sub_obj$Stage)
  
  print(stage_counts)
  
  if (length(stage_counts) < 2) {
    message("Skipping cluster ", cl,
            ": only one stage present")
    next
  }
  
  if (any(stage_counts < 3)) {
    message("Skipping cluster ", cl,
            ": fewer than 3 cells in one stage")
    next
  }
  
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
    tibble::rownames_to_column("gene") %>%
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

library(dplyr)
library(openxlsx)

Idents(SeuObj) <- "seurat_clusters"
DefaultAssay(SeuObj) <- "RNA"

clusters <- levels(Idents(SeuObj))
de_list <- list()

for (cl in clusters) {
  
  message("Processing cluster: ", cl)
  
  cells_cl <- WhichCells(
    SeuObj,
    idents = cl,
    expression = Stage %in% c("Fibrosis", "No_fibrosis")
  )
  
  # Skip empty/small clusters
  if (length(cells_cl) < 30) {
    message("Skipping cluster ", cl, ": <30 cells")
    next
  }
  
  sub_obj <- subset(SeuObj, cells = cells_cl)
  
  # Check group sizes
  stage_counts <- table(sub_obj$Stage)
  
  print(stage_counts)
  
  if (length(stage_counts) < 2) {
    message("Skipping cluster ", cl,
            ": only one stage present")
    next
  }
  
  if (any(stage_counts < 3)) {
    message("Skipping cluster ", cl,
            ": fewer than 3 cells in one stage")
    next
  }
  
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
    tibble::rownames_to_column("gene") %>%
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
