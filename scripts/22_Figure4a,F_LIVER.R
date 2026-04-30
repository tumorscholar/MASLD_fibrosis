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

setwd("/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells")

#### Prepare required input data for CellChat analysis ####

SeuObj <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/SharedTcells&NonTcell.rds")

DefaultAssay(SeuObj) <- "RNA"

Idents(SeuObj) <- "Tissue"

# Isolate LIVER tissue
LIVERObj <- subset(SeuObj, Tissue == "LIVER")

# Rename the column "Sample" to "samples"
colnames(LIVERObj@meta.data)[colnames(LIVERObj@meta.data) == "Sample"] <- "samples"

Idents(LIVERObj) <- "cluster_cellchat"

# Extract cluster and stage metadata
meta <- data.frame(
 cluster = Idents(LIVERObj),
 Stage = LIVERObj$Stage
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
liverObjcellchat <- subset(x = LIVERObj, idents = c("C13", "C23", "C24", "C26", "C27", "C28", "C29"), invert = TRUE)

data.input <- liverObjcellchat[["RNA"]]$data # normalized data matrix
labels <- Idents(liverObjcellchat)

# create a dataframe of the cell labels
meta <- data.frame(labels = labels, row.names = names(labels))

# Split obj stage wise for comparison
objList <- SplitObject(liverObjcellchat, split.by = 'Stage')

#### Create a CellChat object for Healthy ####
cellchat<- objList$Healthy
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

saveRDS(cellchatHliver, file = "cellchatLiverTcH.rds")

#### Create a CellChat object for No_fibrosis ####
cellchat<- objList$No_fibrosis
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

saveRDS(cellchatF0liver, file = "cellchatLiverTcF0.rds")

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

saveRDS(cellchat, file = "cellchatLiverTcF123.rds")


#### Load CellChat object of different dataset ####
cellchatHliver <- readRDS("cellchatLiverTcH.rds")
cellchatF0liver <- readRDS("cellchatLiverTcF0.rds")
cellchatF123liver <- readRDS("cellchatLiverTcF123.rds")

# Merge datasets for comparison
object.list <- list(Healthy = cellchatHliver, No_fibrosis = cellchatF0liver, Fibrosis = cellchatF123liver)
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

ggsave("/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells/Differntial_interaction_strength_2D_liverTc.png", dpi = 300, width = 18, height = 6)


# Merge disease datasets for comparison
object.list <- list(No_fibrosis = cellchatF0liver, Fibrosis = cellchatF123liver)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Circle plot showing differential number of interactions or interaction strength 
gg1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
gg2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

svg(filename = "/data/home/hdx044/plots/cellchat/liver/expandedTcells/Differntial_Interaction_liver.svg", width = 5, height = 5)
print(gg1)
dev.off()

svg(filename = "/data/home/hdx044/plots/cellchat/liver/expandedTcells/Differntial_InteractionStrength_liver.svg", width = 5, height = 5)
print(gg2)
dev.off()

# Heatmap showing differential number of interactions or interaction strength
# Generate the two heatmaps
ht1 <- netVisual_heatmap(cellchat)  # Default is "count"
ht1
# Open SVG device to save the output
svg(filename = "/data/home/hdx044/plots/cellchat/liver/expandedTcells/Signalling_Heatmap_Tc_liver.svg", width = 6, height = 6)
print(ht1)
dev.off()


ht2 <- netVisual_heatmap(cellchat, measure = "weight")  # Default is "count"
ht2
# Open SVG device to save the output
svg(filename = "/data/home/hdx044/plots/cellchat/liver/expandedTcells/Signalling_Heatmap2_TcTr_liver.svg", width = 6, height = 6)
print(ht2)
dev.off()


#### Identify the up-regulated and down-regulated signaling ligand-receptor pairs ####
ptm = Sys.time()
netVisual_bubble(cellchat, sources.use = c(7,13), targets.use = c(1),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = c(7,13), targets.use = c(1),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signalling in LIVER", angle.x = 45, remove.isolate = T, dot.size.min = 5,
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
svg("/data/home/hdx044/plots/cellchat/liver/expandedTcells/EndothelialTcSignallingLiver.svg",
    width = 16, height = 5)
print(gg1)
dev.off()


ptm = Sys.time()
netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(7,13),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(7,13),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signalling in LIVER", angle.x = 45, remove.isolate = T, dot.size.min = 5,
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
svg("/data/home/hdx044/plots/cellchat/liver/expandedTcells/TcEndothelialSignallingLiver.svg",
    width = 16, height = 5)
print(gg1)
dev.off()


# End of the script

