################################################################################
# Script name: ExtendedDataFigure8.R
#
# CellChat heatmap visualisation of differential signalling
# across disease states in LIVER, SAT, and VAT
#
# This script performs:
#  - Loading of precomputed CellChat objects for shared
#    expanded T cells in LIVER, SAT, and VAT tissues
#  - Pairwise merging of CellChat objects between
#    No Fibrosis and Fibrosis conditions for each tissue
#  - Visualisation of differential numbers of inferred
#    cell–cell interactions using CellChat heatmaps
#
# Input files:
#  - LIVER:
#      cellchatLiverTcNoFibrosis.rds
#      cellchatLiverTcFibrosis.rds
#  - SAT:
#      cellchatSATTcNoFibrosis.rds
#      cellchatSATTcFibrosis.rds
#  - VAT:
#      cellchatVATTcNoFibrosis.rds
#      cellchatVATTcFibrosis.rds
#
# Output files:
#  - Signalling_Heatmap_Tc_liver.svg
#  - Signalling_Heatmap_Tc_SAT.svg
#  - Signalling_Heatmap_Tc_VAT.svg
#
# Notes:
#  - All CellChat objects were generated upstream on HPC
#    using memory-safe, sequential execution
#  - This script performs no CellChat inference; it is
#    restricted to comparative visualisation
################################################################################

library(CellChat)
library(ggplot2)
options(stringsAsFactors = FALSE)

# Load liver cellchat object of NoFibrosis and Fibrosis
cellchatF0liver <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells/cellchatLiverTcNoFibrosis.rds")
cellchatF123liver <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells/cellchatLiverTcFibrosis.rds")

# Merge datasets for comparison
object.list <- list(No_fibrosis = cellchatF0liver, Fibrosis = cellchatF123liver)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Heatmap showing differential number of interactions
ht1 <- netVisual_heatmap(cellchat)
ht1
# Open SVG device to save the output
svg(filename = "/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells/Signalling_Heatmap_Tc_liver.svg", width = 6, height = 6)
print(ht1)
dev.off()

# Load SAT cellchat object of Healthy NoFibrosis and Fibrosis
cellchatF0sat <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/SAT/sharedTcells/cellchatSATTcNoFibrosis.rds")
cellchatF123sat <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/SAT/sharedTcells/cellchatSATTcFibrosis.rds")

# Merge datasets for comparison
object.list <- list(NoFibrosis = cellchatF0sat, Fibrosis = cellchatF123sat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Heatmap showing differential number of interactions
ht2 <- netVisual_heatmap(cellchat)
ht2
# Open SVG device to save the output
svg(filename = "/data/Blizard-AlazawiLab/rk/cellchat/SAT/sharedTcells/Signalling_Heatmap_Tc_SAT.svg", width = 6, height = 6)
print(ht2)
dev.off()

# Load VAT cellchat object of NoFibrosis and Fibrosis
cellchatF0 <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/VAT/sharedTcells/cellchatVATTcNoFibrosis.rds")
cellchatF123 <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/VAT/sharedTcells/cellchatVATTcFibrosis.rds")

# Merge datasets for comparison
object.list <- list(No_fibrosis = cellchatF0, Fibrosis = cellchatF123)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Heatmap showing differential number of interactions
ht3 <- netVisual_heatmap(cellchat)
ht3
# Open SVG device to save the output
svg(filename = "//data/Blizard-AlazawiLab/rk/cellchat/VAT/sharedTcells/Signalling_Heatmap_Tc_VAT.svg", width = 6, height = 6)
print(ht3)
dev.off()

# End of the script
