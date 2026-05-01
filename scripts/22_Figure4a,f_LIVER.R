################################################################################
## Script name: 22_Figure4a,f_LIVER.R
##
## CellChat visualisation and comparative analysis – LIVER
##
## This script performs:
## - Loading of precomputed CellChat objects for liver tissue
##   across Healthy, No Fibrosis, and Fibrosis conditions
## - Merging of CellChat objects for cross‑condition comparison
## - Visualisation of global sending/receiving signalling roles
##   using 2D signaling role scatter plots (Fig. 4a)
## - Identification and visualisation of differential
##   ligand–receptor interactions between No Fibrosis and
##   Fibrosis conditions using bubble plots (Fig. 4f)
##
##   Input files:
## - cellchatLiverTcHealthy.rds
## - cellchatLiverTcNoFibrosis.rds
## - cellchatLiverTcFibrosis.rds
##
##   Output files:
## - Differntial_interaction_strength_2D_liverTc.png
## - EndothelialTcSignallingLiver.svg
## - TcEndothelialSignallingLiver.svg
##
##   Notes:
## - All CellChat objects were generated upstream on HPC
##   using sequential execution to ensure memory stability
## - This script performs no CellChat inference, only
##   visualisation and comparative analyses
################################################################################

library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)
library(future)
options(stringsAsFactors = FALSE)

# Load CellChat object of different dataset
cellchatHliver <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells/cellchatLiverTcHealthy.rds")
cellchatF0liver <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells/cellchatLiverTcNoFibrosis.rds")
cellchatF123liver <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells/cellchatLiverTcFibrosis.rds")

####Fig.4.a ####

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

####Fig.4.f ####

# Identify the up-regulated and down-regulated signaling ligand-receptor pairs
ptm = Sys.time()
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
svg("/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells/EndothelialTcSignallingLiver.svg",
    width = 16, height = 5)
print(gg1)
dev.off()

ptm = Sys.time()

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
svg("/data/Blizard-AlazawiLab/rk/cellchat/liver/sharedTcells/TcEndothelialSignallingLiver.svg",
    width = 16, height = 5)
print(gg1)
dev.off()

# End of the script

