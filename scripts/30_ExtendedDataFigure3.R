################################################################################
# Script name: 30_ExtendedDataFigure3.R
#
# Quantification of tissue‑specific T cell subsets
# from a Seurat single‑cell RNA‑seq object
#
# (a) PBMC CXCR3⁺ TBX21⁺ T cells
#     – counts and percentage per patient
#     – stratified by disease stage
#
# (b) Liver CD8⁺ CD69⁺ CXCR6⁺ T cells
#     – counts and percentage per patient
#     – stratified by disease stage
#
# Definitions:
#   Positive cells defined as gene expression > 0.25
#   Percentages calculated relative to total cells
#     within the corresponding tissue
#
# Input:
#   Seurat object: SeuObjx.rds
#
# Output:
#   CSV files for plotting in GraphPad Prism 11
#
################################################################################

# Load packages
library(Seurat)
library(tidyverse)

#Load seurat object
SeuObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

DefaultAssay(SeuObj) <- "RNA"

# Subset to PBMC
pbmc_obj <- subset(
  SeuObj,
  subset = Tissue == "PBMC"
)

# Fetch expression + metadata
expr_df <- FetchData(
  pbmc_obj,
  vars = c("CXCR3", "TBX21", "Patient_ID", "Stage")
)

# Define CXCR3+ TBX21+ cells (expression > 0)
expr_df <- expr_df %>%
  mutate(
    CXCR3_TBX21_pos = CXCR3 > 0.25 & TBX21 > 0.25
  )

# Count + calculate percentage
cxcr3_tbx21_summary <- expr_df %>%
  group_by(Patient_ID, Stage) %>%
  summarise(
    CXCR3_TBX21_Pos_Cells = sum(CXCR3_TBX21_pos),
    Total_PBMC_Cells     = n(),
    Percent_CXCR3_TBX21  = round(
      CXCR3_TBX21_Pos_Cells / Total_PBMC_Cells * 100, 2
    ),
    .groups = "drop"
  )

print(cxcr3_tbx21_summary)

## Extended Data Figure 3a ##
write.csv(
  cxcr3_tbx21_summary,
  "/data/home/hdx044/MASLD_fibrosis/files/PBMC_CXCR3_TBX21_counts_and_percentage_per_patient.csv",
  row.names = FALSE
)

# Data were plotted in GraphPad Prism 11.

# Subset to LIVER cells
liver_obj <- subset(
  SeuObj,
  subset = Tissue == "LIVER"
)

# Fetch expression + metadata
expr_df <- FetchData(
  liver_obj,
  vars = c("CD8A", "CD69", "CXCR6", "Patient_ID", "Stage")
)

# Define CD8+ CD69+ CXCR6+ cells
expr_df <- expr_df %>%
  mutate(
    CD8_CD69_CXCR6_pos =
      CD8A > 0.25 & CD69 > 0.25 & CXCR6 > 0.25
  )

# Count + calculate percentage
liver_cd8_cd69_cxcr6_summary <- expr_df %>%
  group_by(Patient_ID, Stage) %>%
  summarise(
    CD8_CD69_CXCR6_Pos_Cells = sum(CD8_CD69_CXCR6_pos),
    Total_LIVER_Cells       = n(),
    Percent_CD8_CD69_CXCR6  = round(
      CD8_CD69_CXCR6_Pos_Cells / Total_LIVER_Cells * 100, 2
    ),
    .groups = "drop"
  )

print(liver_cd8_cd69_cxcr6_summary)

## Extended Data Figure 3b ##
write.csv(
  liver_cd8_cd69_cxcr6_summary,
  "/data/home/hdx044/MASLD_fibrosis/files/LIVER_CD8_CD69_CXCR6_counts_and_percentage_per_patient.csv",
  row.names = FALSE
)
# Data were plotted in GraphPad Prism 11.


# End of the script