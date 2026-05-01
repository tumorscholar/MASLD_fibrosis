###############################################################################
# Script name: 17_Figure3-d-i.R
#  
# Regulon activity profiling of expanded T cells using AUCell
#
# This script performs:
#  - AUCell-based scoring of SCENIC-derived regulons
#  - Comparison of regulon activity between:
#      (i) Expanded T cells shared across solid tissues vs
#          tissue-restricted expanded clones
#     (ii) Shared expanded T cells from Fibrosis vs No Fibrosis
#  - Visualization of differential regulon activity using
#    heatmaps with a common color scale
#
# Input:
#  - Seurat object containing expanded T cells
#  - SCENIC regulons in GMT format
#
# Output:
#  - Heatmaps (Fig. 3i and Fig. 3j)
#  - AUCell score matrix (CSV)
#
# Notes:
#  - Regulon activity is quantified using AUCell AUC scores
#  - Statistical comparisons use Wilcoxon rank-sum tests with
#    FDR correction
###############################################################################

library(Seurat)        
library(AUCell)        
library(data.table)    
library(dplyr)         
library(pheatmap)      
library(Matrix)     

# Load regulons from GMT
regulon_gmt <- "/data/Blizard-AlazawiLab/rk/scenicTcell/results/regulons_from_ctx.gmt"
regulons <- readLines(regulon_gmt) %>%
 strsplit("\t") %>%
 lapply(function(x) list(name = x[1], genes = x[3:length(x)]))

# Convert to a named list
reg_list <- setNames(lapply(regulons, function(x) x$genes),
                     sapply(regulons, function(x) x$name))

# load obj
obj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinal.rds')

obj <- subset(obj, subset = expanded_shared_class != "a_PBMC")

# Combine tissue and clone classification
obj$tissue.type <- paste0(obj$Tissue, obj$expanded_shared_class, "-")

# Ensure it's genes x cells (rows = genes, columns = cells)

ex_mtx <- as.matrix(GetAssayData(obj, slot="counts"))

# Build rankings for AUCell
cells_rankings <- AUCell_buildRankings(ex_mtx, nCores=1, plotStats=TRUE)

# Calculate AUCell scores
cells_AUC_mat <- AUCell_calcAUC(reg_list, cells_rankings, aucMaxRank = nrow(ex_mtx) * 0.05)

# cells_AUC_mat: rows = regulons, columns = cell barcodes

# Save AUCell matrix
write.csv(as.data.frame(getAUC(cells_AUC_mat)),
          "/data/Blizard-AlazawiLab/rk/scenicTcell/results/AUCell_scores.csv")

## Fig3.j ##

# Tissues to compare
tissues <- c("LIVER", "SAT", "VAT")

# Store results
regulon_results <- list()

# AUCell matrix
auc_mat <- getAUC(cells_AUC_mat)

# Compute per-tissue differential AUC
for (tiss in tissues) {
 message("Processing tissue: ", tiss)
 
 grp1_label <- paste0(tiss, "b_One_solid_tissue-")
 grp2_label <- paste0(tiss, "c_Shared_between_solid_tissues-")
 
 grp1_cells <- colnames(obj)[obj$tissue.type == grp1_label]
 grp2_cells <- colnames(obj)[obj$tissue.type == grp2_label]
 
 grp1_cells_match <- intersect(grp1_cells, colnames(auc_mat))
 grp2_cells_match <- intersect(grp2_cells, colnames(auc_mat))
 
 message("  grp1 cells: ", length(grp1_cells_match))
 message("  grp2 cells: ", length(grp2_cells_match))
 
 if (length(grp1_cells_match) < 5 | length(grp2_cells_match) < 5) {
  message("  Skipping ", tiss, " (not enough matched cells)")
  next
 }
 
 auc_grp1 <- auc_mat[, grp1_cells_match, drop=FALSE]
 auc_grp2 <- auc_mat[, grp2_cells_match, drop=FALSE]
 
 mean_grp1 <- rowMeans(auc_grp1)
 mean_grp2 <- rowMeans(auc_grp2)
 
 pvals <- sapply(rownames(auc_mat), function(r) {
  wilcox.test(auc_grp1[r, ], auc_grp2[r, ])$p.value
 })
 
 padj <- p.adjust(pvals, method="fdr")
 diff_vals <- mean_grp2 - mean_grp1
 
 regulon_results[[tiss]] <- data.frame(
  tissue = tiss,
  regulon = rownames(auc_mat),
  mean_grp1 = mean_grp1,
  mean_grp2 = mean_grp2,
  diff = diff_vals,
  p_val = pvals,
  padj = padj
 )
}

# Filter significant regulons
filter_significant_regulons <- function(df, padj_thresh = 0.05) {
 df %>% filter(padj < padj_thresh) %>% arrange(padj)
}
regulon_filtered <- lapply(regulon_results, filter_significant_regulons)

# Collect all significant regulons across tissues
all_regulons <- unique(c(
 regulon_filtered$LIVER$regulon,
 regulon_filtered$SAT$regulon,
 regulon_filtered$VAT$regulon
)) %>% trimws()

# Initialize heatmap matrix
heat_mat_shared_vs_single <- matrix(0, nrow = length(tissues), ncol = length(all_regulons),
                                    dimnames = list(tissues, all_regulons))

# Fill matrix with differential AUC values
for(tiss in tissues){
 reg_df <- regulon_filtered[[tiss]] %>% filter(regulon %in% all_regulons)
 if(nrow(reg_df) > 0){
  for(i in 1:nrow(reg_df)){
   heat_mat_shared_vs_single[tiss, reg_df$regulon[i]] <- reg_df$diff[i]
  }
 }
}

# Replace NA with 0
heat_mat_shared_vs_single[is.na(heat_mat_shared_vs_single)] <- 0

# Plot heatmap
pheatmap(
 heat_mat_shared_vs_single,
 cluster_rows = FALSE,
 cluster_cols = TRUE,
 color = colorRampPalette(c("blue", "white", "red"))(50),
 main = "Differential Regulon AUC (Shared vs Single) T cells",
 fontsize_row = 10,
 fontsize_col = 8,
 breaks = seq(min(heat_mat_shared_vs_single),
              max(heat_mat_shared_vs_single),
              length.out = 51)
)


## Fig3.j ##

# AUCell matrix
auc_mat <- getAUC(cells_AUC_mat)

# Tissues to compare
tissues <- c("LIVER", "SAT", "VAT")

# Store results
regulon_results_shared <- list()

# Loop through tissues
for (tiss in tissues) {
 message("Processing tissue: ", tiss)
 
 # Select shared expanded clones
 shared_cells <- colnames(obj)[obj$tissue.type == paste0(tiss, "c_Shared_between_solid_tissues-")]
 
 # Split by fibrosis stage
 grp1_cells <- intersect(shared_cells, colnames(auc_mat)[obj$Stage[colnames(auc_mat)] == "Fibrosis"])
 grp2_cells <- intersect(shared_cells, colnames(auc_mat)[obj$Stage[colnames(auc_mat)] == "No_fibrosis"])
 
 message("  Fibrosis cells: ", length(grp1_cells))
 message("  No fibrosis cells: ", length(grp2_cells))
 
 # Skip if insufficient cells
 if (length(grp1_cells) < 5 | length(grp2_cells) < 5) {
  message("  Skipping ", tiss, " (not enough cells)")
  next
 }
 
 # Subset AUCell
 auc_grp1 <- auc_mat[, grp1_cells, drop=FALSE]
 auc_grp2 <- auc_mat[, grp2_cells, drop=FALSE]
 
 # Compute mean AUCell
 mean_grp1 <- rowMeans(auc_grp1)
 mean_grp2 <- rowMeans(auc_grp2)
 
 # Wilcoxon test per regulon
 pvals <- sapply(rownames(auc_mat), function(r){
  wilcox.test(auc_grp1[r,], auc_grp2[r,])$p.value
 })
 
 # FDR correction
 padj <- p.adjust(pvals, method = "fdr")
 
 # Difference (fibrosis minus no fibrosis)
 diff_vals <- mean_grp1 - mean_grp2
 
 # Results table
 regulon_results_shared[[tiss]] <- data.frame(
  tissue = tiss,
  regulon = rownames(auc_mat),
  mean_fibrosis = mean_grp1,
  mean_no_fibrosis = mean_grp2,
  diff = diff_vals,
  p_val = pvals,
  padj = padj
 )
}

# Filter significant regulons
filter_sig <- function(df, padj_thresh = 0.05){
 df %>% filter(padj < padj_thresh)
}
regulon_filtered_shared <- lapply(regulon_results_shared, filter_sig)

# Collect all significant regulons across tissues
all_regulons <- unique(unlist(lapply(regulon_filtered_shared, function(x) x$regulon)))

# Initialize heatmap matrix
heat_mat_fib_vs_nofib <- matrix(0, nrow = length(tissues), ncol = length(all_regulons),
                                dimnames = list(tissues, all_regulons))

# Fill matrix with diff values
for(tiss in tissues){
 df <- regulon_filtered_shared[[tiss]]
 if(nrow(df) > 0){
  for(i in 1:nrow(df)){
   heat_mat_fib_vs_nofib[tiss, df$regulon[i]] <- df$diff[i]
  }
 }
}

# Replace NA with 0
heat_mat_fib_vs_nofib[is.na(heat_mat_fib_vs_nofib)] <- 0

# Plot heatmap
pheatmap(
 heat_mat_fib_vs_nofib,
 cluster_rows = FALSE,
 cluster_cols = TRUE,
 color = colorRampPalette(c("blue", "white", "red"))(50),
 main = "Differential Regulon AUC (Fibrosis - No fibrosis) Shared T cells",
 fontsize_row = 10,
 fontsize_col = 8,
 breaks = seq(min(heat_mat_fib_vs_nofib),
              max(heat_mat_fib_vs_nofib),
              length.out = 51)
)

## Plot both fig 3.j,i in same scale ##

# Compute global min/max across both matrices
global_min <- min(c(heat_mat_shared_vs_single, heat_mat_fib_vs_nofib), na.rm = TRUE)
global_max <- max(c(heat_mat_shared_vs_single, heat_mat_fib_vs_nofib), na.rm = TRUE)

# Define common breaks
common_breaks <- seq(global_min, global_max, length.out = 51)

# Plot first heatmap
ggp1 <- pheatmap(
  heat_mat_shared_vs_single,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  breaks = common_breaks,
  main = "Differential Regulon AUC\n(Shared vs Single) T cells",  # \n splits title into 2 lines
  fontsize_row = 10,
  fontsize_col = 8,
  main.fontface = "plain"                               
)

setwd("/data/home/hdx044/plots/scenic")
svg("RegulonSharedVsSingle.svg", width = 6, height = 4)
print(ggp1)
dev.off()


# Plot second heatmap with same color axis
ggp2 <- pheatmap(
 heat_mat_fib_vs_nofib,
 cluster_rows = FALSE,
 cluster_cols = TRUE,
 color = colorRampPalette(c("blue", "white", "red"))(50),
 breaks = common_breaks,
 main = "Differential Regulon AUC\n(Fibrosis - No Fibrosis) Shared T cells",
 fontsize_row = 10,
 fontsize_col = 8,
 main.fontface = "plain"
)

ggp2

setwd("/data/home/hdx044/plots/scenic")
svg("RegulonFibrosisVsNoFibrosis.svg", width = 6, height = 4)
print(ggp2)
dev.off()

# End of the script
