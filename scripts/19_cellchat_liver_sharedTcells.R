## ==========================================================
## CellChat analysis of shared expanded T cells in LIVER
## SLURM-compatible, memory-safe (sequential) script
## ==========================================================

library(Seurat)
library(CellChat)
library(dplyr)
library(future)

options(stringsAsFactors = FALSE)

## ----------------------------------------------------------
## IMPORTANT: run sequentially to avoid OOM on LIVER
## ----------------------------------------------------------
plan(sequential)
options(future.globals.maxSize = 20 * 1024^3)  # 20 GB

setwd("/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells")

## ----------------------------------------------------------
## Load Seurat object
## ----------------------------------------------------------
SeuObj <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/SharedTcells&NonTcell.rds")

## Rename Sample -> samples if present (CellChat requirement)
if ("Sample" %in% colnames(SeuObj@meta.data)) {
  colnames(SeuObj@meta.data)[
    colnames(SeuObj@meta.data) == "Sample"
  ] <- "samples"
}

DefaultAssay(SeuObj) <- "RNA"
Idents(SeuObj) <- "Tissue"

## ----------------------------------------------------------
## Subset LIVER tissue only
## ----------------------------------------------------------
LIVERObj <- subset(SeuObj, Tissue == "LIVER")
Idents(LIVERObj) <- "cluster_cellchat"

## ----------------------------------------------------------
## Remove rare clusters
## ----------------------------------------------------------
liverObjcellchat <- subset(
  LIVERObj,
  idents = c("C13","C23","C24","C26","C27","C28","C29"),
  invert = TRUE
)

## ----------------------------------------------------------
## Split by disease stage
## ----------------------------------------------------------
objList <- SplitObject(liverObjcellchat, split.by = "Stage")

## ----------------------------------------------------------
## Load CellChat database
## ----------------------------------------------------------
CellChatDB <- CellChatDB.human

## ----------------------------------------------------------
## Helper function to run CellChat
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
## Run CellChat per stage (only if present)
## ----------------------------------------------------------
if ("Healthy" %in% names(objList)) {
  run_cellchat(objList$Healthy, "cellchatLiverTcHealthy.rds")
}

if ("No_fibrosis" %in% names(objList)) {
  run_cellchat(objList$No_fibrosis, "cellchatLiverTcNoFibrosis.rds")
}

if ("Fibrosis" %in% names(objList)) {
  run_cellchat(objList$Fibrosis, "cellchatLiverTcFibrosis.rds")
}

message("✅ LIVER CellChat analysis completed successfully.")