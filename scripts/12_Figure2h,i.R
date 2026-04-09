############################################################
# Title: 12_Figure2h,i.R
#
# Purpose:
# This script identifies autoaggressive and protective CD8 T-cell populations
# based on RNA expression thresholds and quantifies the proportion of expanded
# versus non-expanded clonotypes across tissues, patients, and disease stages.
#
# Figures:
# Fig2h – Autoaggressive CD8 T cells (CD8A+ CD69+ CXCR6+)
# Fig2i – Protective CD8 T cells (CD8A+ CD69+ ITGAE−)
#
# Input:
# Seurat object containing T cells with clonality annotations
# File: TcellObj.rds
#
# Output:
# AutoaggressiveTcells_TissueStageWise.csv
# ProtectiveTcellsTissueStageWise.csv
#
# Notes:
# RNA expression is used for cell-type gating.
# Percentages are calculated at the cell level.
# Expanded clonotypes represent cells belonging to clonotypes observed in
# more than one cell.
# Output tables are exported for downstream plotting in Prism.
#
############################################################

library(Seurat)
library(patchwork)
library(tidyverse)

# Load T cells object
seurat_obj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/TcellObj.rds')

#### Fig2h ####
DefaultAssay(seurat_obj) <- "RNA"

# Identify autoaggressive CD8 T cells based on RNA expression
autoagressive_cells <- WhichCells(seurat_obj, expression = CD8A > 0.25 & CD69 > 0.25 & CXCR6 > 0.25)

# Subset autoaggressive T cells
AutoaggressiveTcells <- subset(
  seurat_obj,
 cells = autoagressive_cells
)

# Inspect autoaggressive T cell object
AutoaggressiveTcells

# Examine clonality distribution within autoaggressive T cells
table(AutoaggressiveTcells$Clonality)

# Count autoaggressive T cells by Stage, Patient, Tissue, and Clonality
auto_counts <- AutoaggressiveTcells@meta.data %>%
 count(Stage, Patient_ID, Tissue, Clonality, name = "n_auto")

# Count total T cells by Stage, Patient, Tissue, and Clonality
total_counts <- seurat_obj@meta.data %>%
 count(Stage, Patient_ID, Tissue, Clonality, name = "n_total")

# Calculate percentage of autoaggressive T cells
auto_percent <- auto_counts %>%
 left_join(total_counts,
           by = c("Stage", "Patient_ID", "Tissue", "Clonality")) %>%
 mutate(
  percent_autoaggressive = 100 * n_auto / n_total
 ) %>%
 arrange(Stage, Patient_ID, Tissue, Clonality)

# View long-format autoaggressive percentages
auto_percent

# Convert to wide format (Expanded vs Non-expanded)
auto_percent_wide <- auto_percent %>%
 select(Stage, Patient_ID, Tissue, Clonality, percent_autoaggressive) %>%
 pivot_wider(
  names_from = Clonality,
  values_from = percent_autoaggressive
 )

# Final autoaggressive summary table
auto_percent_wide
# Data is plotted in Prism
write_csv(auto_percent_wide, "/data/home/hdx044/files/screpertoire/AutoagressiveTcells_TissueStageWise.csv")


#### Fig2i ####
DefaultAssay(seurat_obj) <- "RNA"

# Define protective cells 
expr <- FetchData(
  seurat_obj,
  vars = c("CD8A", "CD69", "ITGAE"),
  assay = "RNA"
)
# Define protective cells
protective_cells <- rownames(expr)[
  expr$CD8A  > 0.25 &
    expr$CD69  > 0.25 &
    expr$ITGAE < 0.25
]

length(protective_cells)   # 9417

# Add Protective metadata
seurat_obj$Protective <- "No"
seurat_obj$Protective[protective_cells] <- "Yes"

table(seurat_obj$Protective)

# Check clonality within protective cells
table(
  seurat_obj$Protective[protective_cells],
  seurat_obj$Clonality[protective_cells]
)

# Count protective cells by Stage / Patient / Tissue / Clonality
protective_counts <- seurat_obj@meta.data %>%
  filter(Protective == "Yes") %>%
  count(Stage, Patient_ID, Tissue, Clonality, name = "n_protective")

# Count total cells
total_counts <- seurat_obj@meta.data %>%
  count(Stage, Patient_ID, Tissue, Clonality, name = "n_total")

# percentage of protective cells
protective_percent <- protective_counts %>%
  left_join(
    total_counts,
    by = c("Stage", "Patient_ID", "Tissue", "Clonality")
  ) %>%
  mutate(
    percent_protective = 100 * n_protective / n_total
  ) %>%
  arrange(Stage, Patient_ID, Tissue, Clonality)

protective_percent

# Create wide table
protective_percent_wide <- protective_percent %>%
  select(Stage, Patient_ID, Tissue, Clonality, percent_protective) %>%
  pivot_wider(
    names_from  = Clonality,
    values_from = percent_protective
  )

protective_percent_wide

# Data is plotted in Prism
write_csv(auto_percent_wide, "/data/home/hdx044/files/screpertoire/ProtectiveTcellsTissueStageWise.csv")

# End of the script
