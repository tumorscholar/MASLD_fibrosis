################################################################################
# Script name: 35_ExtendedDataFigure10.R
# 
# Expanded Clone Gene Signature in Liver T Cells
# Purpose:
#   - Integrate TCR repertoire information into a processed
#     liver FNA Seurat object using scRepertoire
#   - Subset TCR-positive (CTaa⁺) T cells
#   - Quantify absolute log-normalized RNA expression of
#     expanded clone–associated genes
#   - Visualize expression per patient visit (baseline vs follow-up)
#
# Key Features:
#   -TCR contigs integrated via scRepertoire::combineExpression
#   -Barcode harmonization performed on TCR data
#     (Seurat barcodes left unchanged)
#   -Expression values derived from Seurat RNA assay,
#     slot = "data" (log-normalized, no scaling)
#   -Mean expression calculated per VisitLabel
#   -Heatmap displays absolute log-normalized values
#
# Input:
#   - CellRanger TCR contig CSV files (LIVER FNAs)
#   - Processed Seurat object:
#       /data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_LIVER.rds
#
# Output:
#   - Heatmap (PDF):
#       ExpandedCloneGeneSignatureLogValuesinallLiverTcells.pdf
#
################################################################################

library(Seurat)
library(scRepertoire)
library(dplyr)
library(stringr)

setwd("/data/home/hdx044/files/screpertoire/demux_contig/TCR")

# Read TCR contig files (LIVER FNAs only)

contig_list <- list(
  read.csv("GC-WL-10738-LIVER_LIVER_TCR_contig.csv"),   # Baseline
  read.csv("GC-WL-11570-LIVER_TCR_contig.csv"),         # Followup
  read.csv("GC-WL-10291-1-LIVER_TCR_contig.csv"),       # Baseline
  read.csv("GC-WL-11303-LIVER_TCR_contig.csv"),         # Followup
  read.csv("GC-WL-11040-LIVER_LIVER_TCR_contig.csv"),   # Baseline
  read.csv("GC-WL-11816-LIVER_TCR_contig.csv"),         # Followup
  read.csv("GC-WL-11183-LIVER_LIVER_TCR_contig.csv"),   # Baseline
  read.csv("GC-WL-11937-LIVER_TCR_contig.csv")          # Followup
)

sample_names <- c(
  "Baseline-10738-LIVER",
  "Followup-10738-LIVER",
  "Baseline-10291-1-LIVER",
  "Followup-10291-1-LIVER",
  "Baseline-11040-LIVER",
  "Followup-11040-LIVER",
  "Baseline-11183-LIVER",
  "Followup-11183-LIVER"
)

# Combine TCR data
combined.TCR <- combineTCR(
  contig_list,
  samples = sample_names,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)


# Load processed Seurat object

seuratObj <- readRDS(
  "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_LIVER.rds"
)

# Ensure barcode compatibility
# Keep original column order
orig_names <- colnames(seuratObj)

extract_core <- function(x) sub("_.*$", "", x)
seu_core <- extract_core(orig_names)

colnames(seuratObj) <- paste0(seuratObj$Sample, "_", seu_core)

tcr_to_seurat_map <- c(
  "Baseline-10738-LIVER"   = "GC-WL-10738-LIVER",
  "Followup-10738-LIVER"   = "GC-WL-11570-LIVER",
  
  "Baseline-10291-1-LIVER" = "GC-WL-10291-1-LIVER",
  "Followup-10291-1-LIVER" = "GC-WL-11303-LIVER",
  
  "Baseline-11040-LIVER"   = "GC-WL-11040-LIVER",
  "Followup-11040-LIVER"   = "GC-WL-11816-LIVER",
  
  "Baseline-11183-LIVER"   = "GC-WL-11183-LIVER",
  "Followup-11183-LIVER"   = "GC-WL-11937-LIVER"
)

combined.TCR.fixed <- mapply(
  function(df, oldname) {
    
    newname <- tcr_to_seurat_map[[oldname]]
    if (is.null(newname)) stop("No mapping for: ", oldname)
    
    # rewrite sample column
    df$sample <- newname
    
    # rewrite barcode prefix EXACTLY
    df$barcode <- sub(oldname, newname, df$barcode)
    
    df
  },
  combined.TCR,
  names(combined.TCR),
  SIMPLIFY = FALSE
)

names(combined.TCR.fixed) <- unname(tcr_to_seurat_map[names(combined.TCR)])

tcr_barcodes <- unlist(lapply(combined.TCR.fixed, `[[`, "barcode"))

length(tcr_barcodes)
head(tcr_barcodes)

match_rate <- mean(colnames(seuratObj) %in% tcr_barcodes)
match_rate

# Add TCR information to Seurat
seuratObj <- combineExpression(
  input.data = combined.TCR.fixed,
  sc.data    = seuratObj,
  cloneCall  = "gene",
  chain      = "both",
  group.by   = NULL
)

# Restore original column order
colnames(seuratObj) <- orig_names


# Subset for T cells (cells with TCR data)
Tcells <- subset(seuratObj, subset = !is.na(CTaa))

# Verify
cat("All T cells have CTaa?", all(!is.na(Tcells$CTaa)), "\n")

expanded_genes <- c("FGFBP2","GZMH","GNLY","ADGRG1","FCRL6","GZMB","CX3CR1","S1PR5","NKG7","KLRD1","PRF1","TBX21","PLEK","FCGR3A","ZEB2","CST7","CCL5","EFHD2","CTSW","SPON2","CD8A")

# Calculate average log-normalized expression per patient and visit
avg_exp <- AverageExpression(
  Tcells,
  features = expanded_genes,
  group.by = "VisitLabel",
  assays   = "RNA",
  slot     = "data"   # absolute log-normalized values
)

# Long format for ggplot
exp_df <- as.data.frame(avg_exp$RNA)
exp_df$Gene <- rownames(exp_df)

exp_long <- exp_df |>
  tidyr::pivot_longer(
    cols = -Gene,
    names_to = "VisitLabel",
    values_to = "LogExpression"
  )

visit_order <- c(
  "baseline 10738",  "followup 10738",
  "baseline 10291-1","followup 10291-1",
  "baseline 11040",  "followup 11040",
  "baseline 11183",  "followup 11183"
)

exp_long$VisitLabel <- factor(
  exp_long$VisitLabel,
  levels = visit_order
)

p_heatmap <- ggplot(
  exp_long,
  aes(x = VisitLabel, y = Gene, fill = LogExpression)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(LogExpression, 2)), size = 3) +
  scale_fill_gradient2(
    low = "coral",
    mid = "white",
    high = "red",
    midpoint = median(exp_long$LogExpression)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Expanded clone signature (absolute log-normalized expression)",
    x = "Visit",
    y = "Gene",
    fill = "log-normalized expression"
  )

print(p_heatmap)

pdf("/data/home/hdx044/plots/seurat/liver/FNA/ExpandedCloneGeneSignatureLogValuesinallLiverTcells.pdf",
  width = 6.5,
  height = 6.5
)
print(p_heatmap)
dev.off()

# End of the script
