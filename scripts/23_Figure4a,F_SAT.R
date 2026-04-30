library(Seurat)
library(ggplot2)
library(CellChat)
library(future)
library(NMF)
library(ggalluvial)
library(ggsci)
library(ComplexHeatmap)
library(patchwork)
options(stringsAsFactors = FALSE)
future::plan("multisession", workers = 1)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB per worker

# Set working directory
setwd("/data/Blizard-AlazawiLab/rk/cellchat/SAT/sharedTcells")

#### Data input & processing and initialization of CellChat object ####
# Prepare required input data for CellChat analysis

SeuObj <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/SharedTcells&NonTcell.rds")

# Rename the column "Sample" to "samples"
colnames(SeuObj@meta.data)[colnames(SeuObj@meta.data) == "Sample"] <- "samples"

DefaultAssay(SeuObj) <- "RNA"

Idents(SeuObj) <- "Tissue"

# Isolate Liver tissue
SATObj <- subset(SeuObj, Tissue == "SAT")

# Arrange the clusters in increasing order
Idents(SATObj) <- "cluster_cellchat"

# Extract cluster and stage metadata
meta <- data.frame(
 cluster = Idents(SATObj),
 Stage = SATObj$Stage
)
# Filter to only No_Fibrosis and Fibrosis
meta_filtered <- meta %>%
 filter(Stage %in% c("No_fibrosis", "Fibrosis"))

# Count cells per cluster and stage
cell_counts <- meta_filtered %>%
 group_by(cluster, Stage) %>%
 summarise(count = n(), .groups = "drop")

print(cell_counts, n = Inf)

# Keep common clusters of all Stages for comparison
SATObjcellchat <- subset(x = SATObj, idents = c("C13", "C17", "C19", "C22","C23", "C24", "C25", "C27", "C28"), invert = TRUE)

data.input <- SATObjcellchat[["RNA"]]$data # normalized data matrix
labels <- Idents(SATObjcellchat)

# create a dataframe of the cell labels
meta <- data.frame(labels = labels, row.names = names(labels))

# Split obj stage wise for comparison
objList <- SplitObject(SATObjcellchat, split.by = 'Stage')

#### Create a CellChat object for Healthy ####
cellchat <- objList$Healthy
cellchat <- createCellChat(object = cellchat, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

saveRDS(cellchat, file = "cellchatSATTcH.rds")

#### Create a CellChat object for NoFibrosis ####
cellchat <- objList$No_fibrosis
cellchat <- createCellChat(object = cellchat, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

saveRDS(cellchat, file = "cellchatSATTcF0.rds")

#### Create a CellChat object for Fibrosis ####
cellchat<- objList$Fibrosis
cellchat <- createCellChat(object = cellchat, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchatHsat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

saveRDS(cellchat, file = "cellchatSATcF123.rds")

#### Load CellChat object of Healthy NoFibrosis and Fibrosis ####
cellchatHsat <- readRDS("cellchatSATTcH.rds")
cellchatF0sat <- readRDS("cellchatSATTcF0.rds")
cellchatF123sat <- readRDS("cellchatSATcF123.rds")

# Merge F0 and F1_2_3 datasets for comparison
object.list <- list(Healthy = cellchatHsat, NoFibrosis = cellchatF0sat, Fibrosis = cellchatF123sat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Identify cell populations with significant changes in sending or receiving signals
num.link <- unlist(sapply(object.list, function(x) {
 rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
}))

weight.MinMax <- c(min(num.link, na.rm = TRUE), max(num.link, na.rm = TRUE)) 


gg <- list()
for (i in 1:length(object.list)) {
 gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                              title = names(object.list)[i], 
                                              weight.MinMax = weight.MinMax, 
                                              label.size = 6, font.size = 14, font.size.title = 14 ) +  
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) + 
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
  
  theme_classic(base_size = 12) +
  theme(
   axis.text.x = element_text(angle = 45, hjust = 1),
   axis.line = element_line(color = "black"),
   legend.position = "right",
   legend.direction = "vertical",
   legend.title = element_text(size = 10),
   legend.text = element_text(size = 8),
   legend.key.size = unit(1.5, "lines")
  )
}
ggp <- patchwork::wrap_plots(plots = gg)
ggp

ggsave("/data/home/hdx044/plots/cellchat/SAT/expandedTcells/Differntial_interaction_strength_2D_Tc_SAT.png", dpi = 300, width = 18, height = 6)

# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2

svg(filename = "/data/home/hdx044/plots/cellchat/SAT/expandedTcells/Infered_Interaction_SAT.svg", width = 3, height = 6)
print(gg1)
dev.off()

svg(filename = "/data/home/hdx044/plots/cellchat/SAT/expandedTcells/Infered_InteractionStrength_SAT.svg", width = 3, height = 6)
print(gg2)
dev.off()

# Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)

# combining all the identified signaling pathways from different datasets
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i+1]]@netP$pathways, object.list[[i+2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 40, color.heatmap = "Blues")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 40, color.heatmap = "Blues")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 10, height = 40, color.heatmap = "Blues")

combined_heatmap <- draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm", ))

# Save as SVG
svg("/data/home/hdx044/plots/cellchat/SAT/expandedTcells/TcincomingSignallingSAT.svg",
    width = 18, height = 18)
print(combined_heatmap)
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 40, color.heatmap = "Oranges")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 40, color.heatmap = "Oranges")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 40, color.heatmap = "Oranges")

combined_heatmap <- draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm", ))

# Save as SVG
svg("/data/home/hdx044/plots/cellchat/SAT/expandedTcells/TcOutgoingSignallingSAT.svg",
    width = 18, height = 18)
print(combined_heatmap)
dev.off()

# Merge F0 and F1_2_3 datasets for comparison
object.list <- list(NoFibrosis = cellchatF0sat, Fibrosis = cellchatF123sat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Circle plot showing differential number of interactions or interaction strength 
gg1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
gg2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

svg(filename = "/data/home/hdx044/plots/cellchat/SAT/expandedTcells/Differntial_Interaction_SAT.svg", width = 5, height = 5)
print(gg1)
dev.off()

svg(filename = "/data/home/hdx044/plots/cellchat/SAT/expandedTcells/Differntial_InteractionStrength_SAT.svg", width = 5, height = 5)
print(gg2)
dev.off()


# plot
# Generate the two heatmaps
ht1 <- netVisual_heatmap(cellchat)  # Default is "count"
ht2 <- netVisual_heatmap(cellchat, measure = "weight")  # For interaction strength

# Combine the heatmaps side by side with a gap
combined_heatmap <- ht1 + ht2

# Open SVG device to save the output
svg(filename = "/data/home/hdx044/plots/cellchat/SAT/expandedTcells/Signalling_Heatmap1_TcTr_SAT.svg", width = 6, height = 6)
print(ht1)
dev.off()

# Open SVG device to save the output
svg(filename = "/data/home/hdx044/plots/cellchat/SAT/expandedTcells/Signalling_Heatmap2_TcTr_SAT.svg", width = 6, height = 6)
print(ht2)
dev.off()

#### Identify the up-regulated and down-regulated signaling ligand-receptor pairs ####
ptm = Sys.time()
netVisual_bubble(cellchat, sources.use = c(7,13), targets.use = c(1),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = c(7,13), targets.use = c(1),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signalling in SAT", angle.x = 45, remove.isolate = T, dot.size.min = 4,
                        dot.size.max = 5,
                        font.size = 16,
                        font.size.title = 14,
                        color.text = c("blue", "red"))

# Increase legend font size
gg1 <- gg1 + theme(
 legend.title = element_text(size = 16),   # title font
 legend.text  = element_text(size = 14),   # text font
 legend.key.size = unit(1.5, "lines")      # enlarge legend symbols
)

gg1 <- netVisual_bubble(cellchat, sources.use = c(7,13), targets.use = c(1),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signalling in SAT", angle.x = 45, remove.isolate = T, dot.size.min = 5,
                        dot.size.max = 8,
                        font.size = 16,
                        font.size.title = 14,
                        color.text = c("black", "black"))
gg1 <- gg1 +
 coord_flip() +
 scale_color_gradient(
  low = "mistyrose",  # light
  high = "maroon",     # dark maroon
  name = "Communication probability") +
 guides(
  color = guide_colorbar(order = 1),  # heatmap bar
  size  = guide_legend(order = 2)     # p-value dots
 ) +
 theme(
  legend.position = "right",
  legend.direction = "vertical",
  legend.box = "vertical",
  legend.title = element_text(size = 16),
  legend.text  = element_text(size = 14),
  legend.key.size = unit(1.5, "lines")
 )

gg1

# Save as SVG
svg("/data/home/hdx044/plots/cellchat/SAT/expandedTcells/EndothelialTcSignallingSAT.svg",
    width = 16, height = 5)
print(gg1)
dev.off()


ptm = Sys.time()
netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(7,13),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(7,13),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signalling in SAT", angle.x = 45, remove.isolate = T, dot.size.min = 4,
                        dot.size.max = 5,
                        font.size = 16,
                        font.size.title = 14,
                        color.text = c("blue", "red"))

# Increase legend font size
gg1 <- gg1 + theme(
 legend.title = element_text(size = 16),   # title font
 legend.text  = element_text(size = 14),   # text font
 legend.key.size = unit(1.5, "lines")      # enlarge legend symbols
)

gg1 <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(7,13),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signalling in SAT", angle.x = 45, remove.isolate = T, dot.size.min = 5,
                        dot.size.max = 8,
                        font.size = 16,
                        font.size.title = 14,
                        color.text = c("black", "black"))
gg1 <- gg1 +
 coord_flip() +
 scale_color_gradient(
  low = "mistyrose",  # light
  high = "maroon",     # dark maroon
  name = "Communication probability") +
 guides(
  color = guide_colorbar(order = 1),  # heatmap bar
  size  = guide_legend(order = 2)     # p-value dots
 ) +
 theme(
  legend.position = "right",
  legend.direction = "vertical",
  legend.box = "vertical",
  legend.title = element_text(size = 16),
  legend.text  = element_text(size = 14),
  legend.key.size = unit(1.5, "lines")
 )

gg1
# Save as SVG
svg("/data/home/hdx044/plots/cellchat/SAT/expandedTcells/TcEndothelialSignallingSAT.svg",
    width = 16, height = 5)
print(gg1)
dev.off()




