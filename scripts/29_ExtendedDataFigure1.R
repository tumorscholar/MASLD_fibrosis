###############################################################################
# Script name: 29_ExtendedDataFigure1.R
#
# Description:
# This script generate Dot plots of top RNA and ADT markers used to support 
# cluster annotation The input to this
# script is an integrated and clustered Seurat object generated from previous
# preprocessing, integration, and clustering steps.
#
# Purpose:
# To assign biologically interpretable cell‑type labels to clusters and to
# visualise cellular composition and marker expression across tissues in MASLD.
#

###############################################################################

# Load packages
library(Seurat)       
library(SeuratExtend) 
library(dplyr)        
library(ggplot2)       
library(grid)          

#Load seurat object
SeuObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

#### Extended data Fig1.a ####

top_markers <- read.csv(
  "/data/home/hdx044/files/seurat/allTissue/cluster_markers_RNA_top20_allTissue.csv",
  stringsAsFactors = FALSE
)

top5_markers <- top_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

DefaultAssay(SeuObj) <- "RNA"
Idents(SeuObj) <- "cluster"

ggp <- DotPlot2(SeuObj, features = unique(top5_markers$gene), color_scheme = "BuRd")

ggp <- ggp + 
 labs(title = "Top5 genes of each cluster")+
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
 filename = "Top5genesofeachcluster.svg",
 plot = ggp,
 device = "svg",
 path = "/data/home/hdx044/plots/seurat/allTissue",
 width = 8,     
 height = 20,
 units = "in",
 dpi = 600
)


#### Extended data Fig1.b ####

top_markersADT <- read.csv(
  "/data/home/hdx044/files/seurat/allTissue/cluster_markers_ADT_top10_allTissue.csv",
  stringsAsFactors = FALSE
)

top3_markersADT <- top_markersADT %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3, with_ties = FALSE) %>%
  ungroup()

DefaultAssay(SeuObj) <- "ADTonly"

ggp <- DotPlot2(SeuObj, features = unique(top3_markersADT$gene), color_scheme = "BuRd")

ggp <- ggp + 
 labs(title = "Top3 ADTs of each cluster")+
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
 filename = "Top3ADTsofeachcluster.svg",
 plot = ggp,
 device = "svg",
 path = "/data/home/hdx044/plots/seurat/allTissue",
 width = 8,     
 height = 10,
 units = "in",
 dpi = 600
)

# End of the script
