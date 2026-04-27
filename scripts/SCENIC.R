# -----------------------------
# 1. Load libraries
# -----------------------------
library(AUCell)
library(dplyr)
library(data.table)

# -----------------------------
# 2. Load regulons from GMT
# -----------------------------
regulon_gmt <- "/data/Blizard-AlazawiLab/rk/scenic/motifs/regulons_from_ctx.gmt"
regulons <- readLines(regulon_gmt) %>%
 strsplit("\t") %>%
 lapply(function(x) list(name = x[1], genes = x[3:length(x)]))

# Convert to a named list
reg_list <- setNames(lapply(regulons, function(x) x$genes),
                     sapply(regulons, function(x) x$name))

# -----------------------------
# 3. Load expression matrix
# -----------------------------

# load obj
obj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinal.rds')

obj <- subset(obj, subset = expanded_shared_class != "a_PBMC")

# Combine tissue and clone classification
obj$tissue.type <- paste0(obj$Tissue, obj$expanded_shared_class, "-")

# Ensure it's genes x cells (rows = genes, columns = cells)
# Example if from Seurat:

ex_mtx <- as.matrix(GetAssayData(obj, slot="counts"))

# -----------------------------
# 4. Build rankings for AUCell
# -----------------------------
cells_rankings <- AUCell_buildRankings(ex_mtx, nCores=1, plotStats=TRUE)

# -----------------------------
# 5. Calculate AUCell scores
# -----------------------------
cells_AUC_mat <- AUCell_calcAUC(reg_list, cells_rankings, aucMaxRank = nrow(ex_mtx) * 0.05)

# cells_AUC_mat: rows = regulons, columns = cell barcodes

# -----------------------------
# 6. Save AUCell matrix
# -----------------------------
write.csv(as.data.frame(getAUC(cells_AUC_mat)),
          "/data/Blizard-AlazawiLab/rk/scenic/AUCell_scores.csv")


#### Shared vs Single Clone Regulon Analysis ####

library(dplyr)
library(pheatmap)

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


#### Shared Clone Analysis: Fibrosis vs No Fibrosis ####

library(dplyr)
library(pheatmap)

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
 grp1_cells <- intersect(shared_cells, colnames(auc_mat)[obj$Stage[colnames(auc_mat)] == "a_Fibrosis"])
 grp2_cells <- intersect(shared_cells, colnames(auc_mat)[obj$Stage[colnames(auc_mat)] == "b_No_Fibrosis"])
 
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


# Assume you already have:
# heat_mat_shared_vs_single
# heat_mat_fib_vs_nofib

# Step 1: Compute global min/max across both matrices
global_min <- min(c(heat_mat_shared_vs_single, heat_mat_fib_vs_nofib), na.rm = TRUE)
global_max <- max(c(heat_mat_shared_vs_single, heat_mat_fib_vs_nofib), na.rm = TRUE)

# Step 2: Define common breaks
common_breaks <- seq(global_min, global_max, length.out = 51)

library(pheatmap)
library(officer)
library(rvg)
library(grid)

# Step 3: Plot first heatmap
ggp1 <- pheatmap(
 heat_mat_shared_vs_single,
 cluster_rows = FALSE,
 cluster_cols = TRUE,
 color = colorRampPalette(c("blue", "white", "red"))(50),
 breaks = common_breaks,
 main = "Differential Regulon AUC (Shared vs Single) T cells",
 fontsize_row = 10,
 fontsize_col = 8
)

# Step 4: Plot second heatmap with same color axis
ggp2 <- pheatmap(
 heat_mat_fib_vs_nofib,
 cluster_rows = FALSE,
 cluster_cols = TRUE,
 color = colorRampPalette(c("blue", "white", "red"))(50),
 breaks = common_breaks,
 main = "Differential Regulon AUC (Fibrosis - No fibrosis) Shared T cells",
 fontsize_row = 10,
 fontsize_col = 8
)

ppt <- read_pptx()
ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")

ppt <- ph_with(
 ppt,
 dml(code = {
  grid.newpage()
  grid.draw(ggp1$gtable)
 }),
 location = ph_location(
  left   = 0,
  top    = 0,
  width  = 5,
  height = 3
 )
)

ppt <- ph_with(
 ppt,
 dml(code = {
  grid.newpage()
  grid.draw(ggp2$gtable)
 }),
 location = ph_location(
  left   = 5,
  top    = 0,
  width  = 5,
  height = 3
 )
)

print(
 ppt,
 target = "/data/home/hdx044/plots/scenic/Regulon_heatmaps_Tcells.pptx"
)

