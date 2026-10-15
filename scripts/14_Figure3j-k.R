###############################################################################
# Figure 3j,k
#
# Transcription factor activity analysis of tissue-expanded T cells using
# DoRothEA regulons and VIPER.
#
# Figure 3j
# Shared-between-solid-tissues clones vs One-solid-tissue clones
#
# Figure 3k
# Fibrosis vs No Fibrosis within Shared-between-solid-tissues clones
#
# Input:
#   expandedTcellFinalV7.rds
#
# Output:
#   Dorothea_TF_activity.csv
#   Figure3j_significant_TFs.csv
#   Figure3k_significant_TFs.csv
#
# Method:
#   - DoRothEA human regulons (confidence A-C)
#   - VIPER TF activity inference
#   - Wilcoxon rank-sum test
#   - Benjamini-Hochberg FDR correction
#
###############################################################################

source("00_nature_style.R")

out_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure3"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure3"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# Load packages
###############################################################################

library(Seurat)
library(dorothea)
library(viper)
library(dplyr)
library(Matrix)
library(pheatmap)
library(tidyr)

###############################################################################
# Robust heatmap saver: writes PDF + SVG + TIFF for one matrix.
# Each format is opened, drawn and closed independently (tryCatch + guaranteed
# dev.off), so a failure in one device can NEVER stop the others - this is why
# 3j previously produced no TIFF while 3k did. Row clustering is only requested
# when there are >= 2 rows (pheatmap errors on a single-row matrix); columns
# (tissues) are kept in fixed order and never clustered.
###############################################################################
save_heat <- function(mat, file, title) {
  mat <- as.matrix(mat)
  if (nrow(mat) < 1 || ncol(mat) < 1) {
    message("  [save_heat] empty matrix - skipping: ", basename(file)); return(invisible(NULL))
  }
  cr <- nrow(mat) >= 2
  draw <- function()
    pheatmap(mat,
             color        = colorRampPalette(c("blue","white","red"))(100),
             cluster_rows = cr,
             cluster_cols = FALSE,
             main         = title)
  dev_write <- function(open_fun, path) {
    tryCatch({ open_fun(); draw() },
             error = function(e) message("  [save_heat] ", basename(path),
                                         " FAILED: ", conditionMessage(e)))
    try(grDevices::dev.off(), silent = TRUE)
  }
  dev_write(function() grDevices::pdf(paste0(file, ".pdf"),  width = 5, height = 6), paste0(file, ".pdf"))
  dev_write(function() grDevices::svg(paste0(file, ".svg"),  width = 5, height = 6), paste0(file, ".svg"))
  dev_write(function() grDevices::tiff(paste0(file, ".tiff"), width = 5, height = 6,
                                       units = "in", res = 600, compression = "lzw"), paste0(file, ".tiff"))
  message("Saved heatmap (pdf/svg/tiff): ", basename(file))
}

###############################################################################
# Load tissue-expanded T-cell object
#
# PBMC-only clones are excluded to focus on tissue-expanded populations.
###############################################################################

obj <- readRDS(
  "/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinalV7.rds"
)

obj <- subset(
  obj,
  subset = expanded_shared_class != "a_PBMC"
)

DefaultAssay(obj) <- "RNA"

###############################################################################
# Load high-confidence DoRothEA TF-target interactions
#
# Confidence levels:
#   A = highest confidence
#   B = high confidence
#   C = medium confidence
###############################################################################

data(dorothea_hs)

regulons <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

cat(
  "Unique TFs:",
  length(unique(regulons$tf)),
  "\n"
)

###############################################################################
# Infer TF activity using VIPER
###############################################################################

regulon_obj <- df2regulon(regulons)

# viper()/aREA require a dense numeric matrix; a sparse dgCMatrix makes
# is.na(eset) return a Matrix object and aREA fails ("argument to 'which' is
# not logical"). Coerce to a base matrix.
expr_mat <- as.matrix(GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
))


tf_activity <- viper(
  expr_mat,
  regulon_obj,
  minsize = 4,
  eset.filter = FALSE,
  verbose = TRUE
)

write.csv(
  tf_activity,
  file.path(tab_dir, "Dorothea_TF_activity.csv")
)

cat(
  "TF activity matrix:",
  dim(tf_activity)[1],
  "TFs x",
  dim(tf_activity)[2],
  "cells\n"
)

###############################################################################
# Figure 3j
# Shared vs One-solid-tissue clones
###############################################################################

obj$tissue.type <- paste0(
  obj$Tissue,
  obj$expanded_shared_class,
  "-"
)

tissues <- c("LIVER","SAT","VAT")

res_j <- list()

for(tiss in tissues){
  
  grp1 <- colnames(obj)[
    obj$tissue.type ==
      paste0(tiss,"c_Shared_between_solid_tissues-")
  ]
  
  grp2 <- colnames(obj)[
    obj$tissue.type ==
      paste0(tiss,"b_One_solid_tissue-")
  ]
  
  grp1 <- intersect(grp1, colnames(tf_activity))
  grp2 <- intersect(grp2, colnames(tf_activity))
  
  cat(
    tiss,
    " Shared:", length(grp1),
    " One-tissue:", length(grp2),
    "\n"
  )
  
  if(length(grp1) < 5 || length(grp2) < 5)
    next
  
  pvals <- sapply(
    rownames(tf_activity),
    function(tf){
      wilcox.test(
        tf_activity[tf, grp1],
        tf_activity[tf, grp2]
      )$p.value
    }
  )
  
  res_j[[tiss]] <- data.frame(
    TF = rownames(tf_activity),
    diff = as.numeric(
      rowMeans(tf_activity[, grp1, drop = FALSE]) -
        rowMeans(tf_activity[, grp2, drop = FALSE])
    ),
    padj = p.adjust(pvals,"fdr")
  )
}

all_j <- bind_rows(
  lapply(
    names(res_j),
    function(tiss){
      
      res_j[[tiss]] %>%
        mutate(tissue = tiss)
      
    }
  )
)

sig_j <- all_j %>%
  filter(padj < 0.05) %>%
  arrange(padj)

###############################################################################
# Figure 3j
#
# Shared-between-solid-tissues clones
# versus
# One-solid-tissue clones
###############################################################################

write.csv(
  all_j,
  file.path(tab_dir, "Figure3j_all_TFs.csv"),
  row.names = FALSE
)

write.csv(
  sig_j,
  file.path(tab_dir, "Figure3j_significant_TFs.csv"),
  row.names = FALSE
)

print(table(sig_j$tissue))

###############################################################################
# Figure 3k
#
# Fibrosis versus No Fibrosis
# within Shared-between-solid-tissue clones
###############################################################################

stage_named <- obj$Stage
names(stage_named) <- colnames(obj)

shared_of <- function(t){
  
  colnames(obj)[
    obj$tissue.type ==
      paste0(
        t,
        "c_Shared_between_solid_tissues-"
      )
  ]
}

res_k <- list()

for(tiss in tissues){
  
  grp1 <- intersect(
    shared_of(tiss),
    names(stage_named)[
      stage_named == "Fibrosis"
    ]
  )
  
  grp2 <- intersect(
    shared_of(tiss),
    names(stage_named)[
      stage_named == "No_fibrosis"
    ]
  )
  
  grp1 <- intersect(grp1, colnames(tf_activity))
  grp2 <- intersect(grp2, colnames(tf_activity))
  
  cat(
    tiss,
    " Fibrosis:", length(grp1),
    " NoFibrosis:", length(grp2),
    "\n"
  )
  
  if(length(grp1) < 5 || length(grp2) < 5)
    next
  
  pvals <- sapply(
    rownames(tf_activity),
    function(tf){
      wilcox.test(
        tf_activity[tf, grp1],
        tf_activity[tf, grp2]
      )$p.value
    }
  )
  
  res_k[[tiss]] <- data.frame(
    TF = rownames(tf_activity),
    diff = as.numeric(
      rowMeans(tf_activity[, grp1, drop = FALSE]) -
        rowMeans(tf_activity[, grp2, drop = FALSE])
    ),
    padj = p.adjust(pvals,"fdr")
  )
}

all_k <- bind_rows(
  lapply(
    names(res_k),
    function(tiss){
      
      res_k[[tiss]] %>%
        mutate(tissue = tiss)
      
    }
  )
)

sig_k <- all_k %>%
  filter(padj < 0.05) %>%
  arrange(padj)

write.csv(
  all_k,
  file.path(tab_dir, "Figure3k_all_TFs.csv"),
  row.names = FALSE
)

write.csv(
  sig_k,
  file.path(tab_dir, "Figure3k_significant_TFs.csv"),
  row.names = FALSE
)

print(table(sig_k$tissue))

###############################################################################
# Figure 3j heatmap
###############################################################################

top_j <- sig_j %>%
  group_by(tissue) %>%
  slice_head(n = 20)

heat_j <- top_j %>%
  select(TF, tissue, diff) %>%
  distinct() %>%
  pivot_wider(
    names_from = tissue,
    values_from = diff,
    values_fill = 0
  )

heat_j_mat <- as.data.frame(heat_j)
rownames(heat_j_mat) <- heat_j_mat$TF
heat_j_mat$TF <- NULL

save_heat(heat_j_mat, file.path(out_dir, "Figure3j_Dorothea_VIPER"),
          "Shared vs One-solid-tissue")


###############################################################################
# Figure 3k heatmap
###############################################################################

top_k <- sig_k %>%
  group_by(tissue) %>%
  slice_head(n = 20)

heat_k <- top_k %>%
  select(TF, tissue, diff) %>%
  distinct() %>%
  pivot_wider(
    names_from = tissue,
    values_from = diff,
    values_fill = 0
  )

heat_k_mat <- as.data.frame(heat_k)
rownames(heat_k_mat) <- heat_k_mat$TF
heat_k_mat$TF <- NULL

save_heat(heat_k_mat, file.path(out_dir, "Figure3k_Dorothea_VIPER"),
          "Fibrosis vs No Fibrosis")

# End of the script