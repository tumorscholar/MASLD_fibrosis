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
setwd("/data/Blizard-AlazawiLab/rk/cellchat/VAT/sharedTcells")

#### Data input & processing and initialization of CellChat object ####
# Prepare required input data for CellChat analysis

SeuObj <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/SharedTcells&NonTcell.rds")

# Rename the column "Sample" to "samples"
colnames(SeuObj@meta.data)[colnames(SeuObj@meta.data) == "Sample"] <- "samples"

DefaultAssay(SeuObj) <- "RNA"

Idents(SeuObj) <- "Tissue"

# Isolate Liver tissue
VATObj <- subset(SeuObj, Tissue == "VAT")

# Arrange the clusters in increasing order
Idents(VATObj) <- "cluster_cellchat"

# Extract cluster and stage metadata
meta <- data.frame(
 cluster = Idents(VATObj),
 Stage = VATObj$Stage
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
VATObjcellchat <- subset(x = VATObj, idents = c("C17", "C25", "C26", "C27", "C28", "C30"), invert = TRUE)

data.input <- VATObjcellchat[["RNA"]]$data # normalized data matrix
labels <- Idents(VATObjcellchat)
# create a dataframe of the cell labels
meta <- data.frame(labels = labels, row.names = names(labels))

# Split obj stage wise for comparison
objList <- SplitObject(VATObjcellchat, split.by = 'Stage')

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

saveRDS(cellchat, file = "cellchatVATcH.rds")

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

saveRDS(cellchat, file = "cellchatVATcF0.rds")

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
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")


saveRDS(cellchat, file = "cellchatVATcF123.rds")


#### Load CellChat object of Healthy NoFibrosis and Fibrosis dataset ####
cellchatH <- readRDS("cellchatVATcH.rds")
cellchatF0 <- readRDS("cellchatVATcF0.rds")
cellchatF123 <- readRDS("cellchatVATcF123.rds")

# Merge F0 and F123 datasets for comparison
object.list <- list(Healthy = cellchatH, No_fibrosis = cellchatF0, Fibrosis = cellchatF123)
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
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) + 
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
  
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

ggsave("/data/home/hdx044/plots/cellchat/VAT/expandedTcells/Differntial_interaction_strength_2D_VATTc.png", dpi = 300, width = 18, height = 6)


# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2

svg(filename = "/data/home/hdx044/plots/cellchat/VAT/expandedTcells/Infered_Interaction_VAT.svg", width = 3, height = 6)
print(gg1)
dev.off()

svg(filename = "/data/home/hdx044/plots/cellchat/VAT/expandedTcells/Infered_InteractionStrength_VAT.svg", width = 3, height = 6)
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
svg("/data/home/hdx044/plots/cellchat/VAT/expandedTcells/TcincomingSignallingVAT.svg",
    width = 18, height = 18)
print(combined_heatmap)
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 40, color.heatmap = "Oranges")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 40, color.heatmap = "Oranges")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 40, color.heatmap = "Oranges")

combined_heatmap <- draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm", ))

# Save as SVG
svg("/data/home/hdx044/plots/cellchat/VAT/expandedTcells/TcOutgoingSignallingVAT.svg",
    width = 18, height = 18)
print(combined_heatmap)
dev.off()

# Merge datasets for comparison
object.list <- list(No_fibrosis = cellchatF0, Fibrosis = cellchatF123)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Circle plot showing differential number of interactions or interaction strength 
gg1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
gg2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

svg(filename = "/data/home/hdx044/plots/cellchat/VAT/expandedTcells/Differntial_Interaction_VAT.svg", width = 5, height = 5)
print(gg1)
dev.off()

svg(filename = "/data/home/hdx044/plots/cellchat/VAT/expandedTcells/Differntial_InteractionStrength_VAT.svg", width = 5, height = 5)
print(gg2)
dev.off()

# Heatmap showing differential number of interactions or interaction strength
# Generate the two heatmaps
ht1 <- netVisual_heatmap(cellchat)  # Default is "count"
ht1
# Open SVG device to save the output
svg(filename = "/data/home/hdx044/plots/cellchat/VAT/expandedTcells/Signalling_Heatmap_Tc_VAT.svg", width = 6, height = 6)
print(ht1)
dev.off()


ht2 <- netVisual_heatmap(cellchat, measure = "weight")  # Default is "count"
ht2
# Open SVG device to save the output
svg(filename = "/data/home/hdx044/plots/cellchat/VAT/expandedTcells/Signalling_Heatmap2_Tc_VAT.svg", width = 6, height = 6)
print(ht2)
dev.off()


#### Identify the up-regulated and down-regulated signaling ligand-receptor pairs ####
ptm = Sys.time()
netVisual_bubble(cellchat, sources.use = c(7,14), targets.use = c(1),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = c(7,14), targets.use = c(1),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signalling in VAT", angle.x = 45, remove.isolate = T, dot.size.min = 5,
                        dot.size.max = 8,
                        font.size = 16,
                        font.size.title = 14,
                        color.text = c("blue", "red"))

# Increase legend font size
gg1 <- gg1 + theme(
 legend.title = element_text(size = 16),   # title font
 legend.text  = element_text(size = 14),   # text font
 legend.position = "right",
 legend.justification = "center",
 legend.key.size = unit(1.5, "lines")      # enlarge legend symbols
) 

gg1 <- netVisual_bubble(cellchat, sources.use = c(7,14), targets.use = c(1),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signalling in VAT", angle.x = 45, remove.isolate = T, dot.size.min = 5,
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
svg("/data/home/hdx044/plots/cellchat/VAT/expandedTcells/EndothelialTcSignallingVAT.svg",
    width = 16, height = 5)
print(gg1)
dev.off()


ptm = Sys.time()
netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(7,14),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(7,14),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signalling in VAT", angle.x = 45, remove.isolate = T, dot.size.min = 5,
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
svg("/data/home/hdx044/plots/cellchat/VAT/expandedTcells/TcEndothelialSignallingVAT.svg",
    width = 16, height = 5)
print(gg1)
dev.off()



