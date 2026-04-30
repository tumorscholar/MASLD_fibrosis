# Load packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(scRepertoire)

#### Load seurat object ####
SeuObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

# Load CSV
expanded_clones <- read.csv("/data/home/hdx044/files/screpertoire/expanded_clonotypes_3098.csv", stringsAsFactors = FALSE)

# Check the structure
head(expanded_clones)

# Extract the expanded TCR sequences
expanded_ctaa <- expanded_clones$CTaa

# Create a new metadata column
SeuObj$Expanded <- ifelse(SeuObj$CTaa %in% expanded_ctaa, "Expanded", "Non-expanded")

# Keep only non expanded cells
SeuObj_NonExpanded <- subset(SeuObj, subset = Expanded == "Non-expanded")

# Keep only non T cells
# Check which clusters exist
table(Idents(SeuObj_NonExpanded))

# Define the clusters you want to remove (T-cell clusters)
tcell_clusters <- c(0, 1, 4, 7, 21)

#### Subset to keep only non-T-cell clusters ####
SeuObj_NonTcell <- subset(SeuObj_NonExpanded, 
                          idents = setdiff(unique(Idents(SeuObj_NonExpanded)), tcell_clusters))

# Verify the new object
table(Idents(SeuObj_NonTcell))

# Optionally save the new Seurat object
saveRDS(SeuObj_NonTcell, file = "/data/Blizard-AlazawiLab/rk/seurat/SeuObj_NonTcell.rds")

