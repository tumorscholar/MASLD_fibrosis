## ==========================================================
## CellChat analysis of shared expanded T cells in VAT
## SLURM-compatible, memory-safe script
## ==========================================================

library(Seurat)
library(CellChat)
library(dplyr)
library(future)

options(stringsAsFactors = FALSE)

## ----------------------------------------------------------
## IMPORTANT: disable parallelism to avoid OOM + RNG issues
## ----------------------------------------------------------
plan(sequential)
options(
  future.globals.maxSize = 20 * 1024^3,
  future.rng.onMisuse = "ignore"
)

setwd("/data/Blizard-AlazawiLab/rk/cellchat/VAT/sharedTcells")

## ----------------------------------------------------------
## Load Seurat object
## ----------------------------------------------------------
SeuObj <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/SharedTcells&NonTcell.rds")

## Rename Sample -> samples if present
if ("Sample" %in% colnames(SeuObj@meta.data)) {
  colnames(SeuObj@meta.data)[
    colnames(SeuObj@meta.data) == "Sample"
  ] <- "samples"
}

DefaultAssay(SeuObj) <- "RNA"
Idents(SeuObj) <- "Tissue"

## ----------------------------------------------------------
## Subset VAT tissue
## ----------------------------------------------------------
VATObj <- subset(SeuObj, Tissue == "VAT")
Idents(VATObj) <- "cluster_cellchat"

## ----------------------------------------------------------
## Remove rare clusters
## ----------------------------------------------------------
VATObjcellchat <- subset(
  VATObj,
  idents = c("C17","C25","C26","C27","C28","C30"),
  invert = TRUE
)

## ----------------------------------------------------------
## Split by disease stage
## ----------------------------------------------------------
objList <- SplitObject(VATObjcellchat, split.by = "Stage")

CellChatDB <- CellChatDB.human

## ----------------------------------------------------------
## Helper function
## ----------------------------------------------------------
run_cellchat <- function(seu, outfile) {
  cellchat <- createCellChat(seu, group.by = "ident", assay = "RNA")
  cellchat@DB <- CellChatDB
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  saveRDS(cellchat, file = outfile)
}

## ----------------------------------------------------------
## Run per stage (only if present)
## ----------------------------------------------------------
if ("Healthy" %in% names(objList)) {
  run_cellchat(objList$Healthy, "cellchatVATTcHealthy.rds")
}

if ("No_fibrosis" %in% names(objList)) {
  run_cellchat(objList$No_fibrosis, "cellchatVATTcNoFibrosis.rds")
}

if ("Fibrosis" %in% names(objList)) {
  run_cellchat(objList$Fibrosis, "cellchatVATTcFibrosis.rds")
}

message("✅ VAT CellChat analysis completed successfully.")
``