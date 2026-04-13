###############################################################################
# Script name: 11_Figure2e,f.R
#
# Description:
# This script performs multi‑tissue T‑cell receptor (TCR) repertoire analysis
# and integrates clonotype expansion information with single‑cell RNA‑seq data.
# TCR contig annotation files from SAT, VAT, PBMC, and liver samples across
# Healthy, No Fibrosis, and Fibrosis conditions are combined using scRepertoire.
#
# Expanded clonotypes are defined based on shared CDR3 amino‑acid sequences
# (CTaa) detected in more than one cell. Canonical 10x barcodes are extracted
# from TCR contig annotations and mapped back to the Seurat object to label
# individual cells as expanded or non‑expanded.
#
# Key steps performed in this script include:
#  - Loading and combining multi‑tissue TCR contig annotation CSV files
#  - Annotating samples with disease status and tissue identity
#  - Quantifying total T cells and unique TCR clonotypes
#  - Defining expanded, hyperexpanded, and singleton clonotypes
#  - Mapping expanded clonotypes to single cells via barcode matching
#  - Making T cell object
#  - Differential expression analysis of expanded vs non‑expanded T cells
#  - Visualization of expanded clone–associated RNA and ADT markers
#
# All analyses focus on characterising transcriptional and protein‑level
# differences associated with TCR clonal expansion across tissues in MASLD.
#
# Outputs:
#  - T cell object annotated with clonality metadata
#  - Differential expression results for expanded vs non‑expanded T cells
#  - Figures for RNA and ADT marker expression
#
# Purpose:
# To identify and characterise expanded T cell clonotypes and their molecular
# signatures across tissues and fibrosis stages.
#
###############################################################################

# Load required packages
library(Seurat)        
library(tidyverse)     
library(scRepertoire)  
library(patchwork)
library(SeuratExtend)

setwd ('/data/home/hdx044/files/screpertoire/demux_contig/TCR')

#Read files
#Healthy sample
s1 = read.csv("GC-WL-10742-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s2 = read.csv("GC-WL-10742-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s3 = read.csv("GC-WL-10742-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s4 = read.csv("GC-WL-10742-LIVER_LIVER_TCR_contig.csv")
#F0 samples
s5 = read.csv("GC-WL-9961-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s6 = read.csv("GC-WL-9961-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s7 = read.csv("GC-WL-9961-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s8 = read.csv("GC-WL-9961-LIVER_LIVER_TCR_contig.csv")
s9 = read.csv("GC-WL-9999-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s10 = read.csv("GC-WL-9999-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s11 = read.csv("GC-WL-9999-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s12 = read.csv("GC-WL-9999-LIVER_LIVER_TCR_contig.csv")
s13 = read.csv("GC-WL-10113-2-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s14 = read.csv("GC-WL-10113-2-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s15 = read.csv("GC-WL-10113-2-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s16 = read.csv("GC-WL-10113-2-LIVER_LIVER_TCR_contig.csv")
#F1 samples
s17 = read.csv("GC-WL-9680-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s18 = read.csv("GC-WL-9680-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s19 = read.csv("GC-WL-9680-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s20 = read.csv("GC-WL-9680-LIVER_LIVER_TCR_contig.csv")
s21 = read.csv("GC-WL-10203-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s22 = read.csv("GC-WL-10203-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s23 = read.csv("GC-WL-10203-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s24 = read.csv("GC-WL-10203-LIVER_LIVER_TCR_contig.csv")
s25 = read.csv("GC-WL-10113-1-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s26 = read.csv("GC-WL-10113-1-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s27 = read.csv("GC-WL-10113-1-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s28 = read.csv("GC-WL-10113-1-LIVER_LIVER_TCR_contig.csv")
s29 = read.csv("GC-WL-10380-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s30 = read.csv("GC-WL-10380-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s31 = read.csv("GC-WL-10380-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s32 = read.csv("GC-WL-10380-LIVER_LIVER_TCR_contig.csv")
s33 = read.csv("GC-WL-10291-2-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s34 = read.csv("GC-WL-10291-2-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s35 = read.csv("GC-WL-10291-2-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s36 = read.csv("GC-WL-10291-2-LIVER_LIVER_TCR_contig.csv")
s37 = read.csv("GC-WL-10202-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s38 = read.csv("GC-WL-10202-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s39 = read.csv("GC-WL-10202-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s40 = read.csv("GC-WL-10202-LIVER_LIVER_TCR_contig.csv")
s41 = read.csv("GC-WL-10634-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s42 = read.csv("GC-WL-10634-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s43 = read.csv("GC-WL-10634-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s44 = read.csv("GC-WL-10634-LIVER_LIVER_TCR_contig.csv")
s45 = read.csv("GC-WL-10738-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s46 = read.csv("GC-WL-10738-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s47 = read.csv("GC-WL-10738-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s48 = read.csv("GC-WL-10738-LIVER_LIVER_TCR_contig.csv")
#F2 samples
s49 = read.csv("GC-WL-10205-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s50 = read.csv("GC-WL-10205-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s51 = read.csv("GC-WL-10205-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s52 = read.csv("GC-WL-10205-LIVER_LIVER_TCR_contig.csv")
s53 = read.csv("GC-WL-9991-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s54 = read.csv("GC-WL-9991-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s55 = read.csv("GC-WL-9991-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s56 = read.csv("GC-WL-9991-LIVER_LIVER_TCR_contig.csv")
s57 = read.csv("GC-WL-9932-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s58 = read.csv("GC-WL-9932-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s59 = read.csv("GC-WL-9932-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s60 = read.csv("GC-WL-9932-LIVER_LIVER_TCR_contig.csv")
s61 = read.csv("GC-WL-11040-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s62 = read.csv("GC-WL-11040-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s63 = read.csv("GC-WL-11040-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s64 = read.csv("GC-WL-11040-LIVER_LIVER_TCR_contig.csv")
#F3 samples
s65 = read.csv("GC-WL-11051-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s66 = read.csv("GC-WL-11051-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s67 = read.csv("GC-WL-11051-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s68 = read.csv("GC-WL-11051-LIVER_LIVER_TCR_contig.csv")
s69 = read.csv("GC-WL-11183-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s70 = read.csv("GC-WL-11183-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s71 = read.csv("GC-WL-11183-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s72 = read.csv("GC-WL-11183-LIVER_LIVER_TCR_contig.csv")
s73 = read.csv("GC-WL-11471-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s74 = read.csv("GC-WL-11471-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s75 = read.csv("GC-WL-11471-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s76 = read.csv("GC-WL-11471-LIVER_LIVER_TCR_contig.csv")

#F1 samples
s77 = read.csv("GC-WL-11327-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s78 = read.csv("GC-WL-11327-SAT-VAT-PBMC_VAT_TCR_contig.csv")
s79 = read.csv("GC-WL-11327-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s80 = read.csv("GC-WL-11327-LIVER_LIVER_TCR_contig.csv")

#list
contig_list = list(s1,s2,s3,s4, s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45,s46,s47,s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,s60,s61,s62,s63,s64,s65,s66,s67,s68,s69,s70,s71,s72,s73,s74,s75,s76,s77,s78,s79,s80)

combined.TCR = combineTCR(contig_list, samples = c("GC-WL-10742-SAT-VAT-PBMC-SAT", "GC-WL-10742-SAT-VAT-PBMC-VAT", "GC-WL-10742-SAT-VAT-PBMC-PBMC", "GC-WL-10742-LIVER-LIVER",
                                                   "GC-WL-9961-SAT-VAT-PBMC-SAT", "GC-WL-9961-SAT-VAT-PBMC-VAT", "GC-WL-9961-SAT-VAT-PBMC-PBMC", "GC-WL-9961-LIVER-LIVER",
                                                   "GC-WL-9999-SAT-VAT-PBMC-SAT", "GC-WL-9999-SAT-VAT-PBMC-VAT", "GC-WL-9999-SAT-VAT-PBMC-PBMC", "GC-WL-9999-LIVER-LIVER",
                                                   "GC-WL-10113-2-SAT-VAT-PBMC-SAT", "GC-WL-10113-2-SAT-VAT-PBMC-VAT", "GC-WL-10113-2-SAT-VAT-PBMC-PBMC", "GC-WL-10113-2-LIVER-LIVER",
                                                   "GC-WL-9680-SAT-VAT-PBMC-SAT", "GC-WL-9680-SAT-VAT-PBMC-VAT", "GC-WL-9680-SAT-VAT-PBMC-PBMC", "GC-WL-9680-LIVER-LIVER",
                                                   "GC-WL-10203-SAT-VAT-PBMC-SAT", "GC-WL-10203-SAT-VAT-PBMC-VAT", "GC-WL-10203-SAT-VAT-PBMC-PBMC", "GC-WL-10203-LIVER-LIVER",
                                                   "GC-WL-10113-1-SAT-VAT-PBMC-SAT", "GC-WL-10113-1-SAT-VAT-PBMC-VAT", "GC-WL-10113-1-SAT-VAT-PBMC-PBMC", "GC-WL-10113-1-LIVER-LIVER",
                                                   "GC-WL-10380-SAT-VAT-PBMC-SAT", "GC-WL-10380-SAT-VAT-PBMC-VAT", "GC-WL-10380-SAT-VAT-PBM-PBMC", "GC-WL-10380-LIVER-LIVER",
                                                   "GC-WL-10291-2-SAT-VAT-PBMC-SAT", "GC-WL-10291-2-SAT-VAT-PBMC-VAT", "GC-WL-10291-2-SAT-VAT-PBMC-PBMC", "GC-WL-10291-2-LIVER-LIVER",
                                                   "GC-WL-10202-SAT-VAT-PBMC-SAT", "GC-WL-10202-SAT-VAT-PBMC-VAT", "GC-WL-10202-SAT-VAT-PBMC-PBMC", "GC-WL-10202-LIVER-LIVER",
                                                   "GC-WL-10634-SAT-VAT-PBMC-SAT","GC-WL-10634-SAT-VAT-PBMC-VAT","GC-WL-10634-SAT-VAT-PBMC-PBMC", "GC-WL-10634-LIVER-LIVER", 
                                                   "GC-WL-10738-SAT-VAT-PBMC-SAT","GC-WL-10738-SAT-VAT-PBMC-VAT","GC-WL-10738-SAT-VAT-PBMC-PBMC","GC-WL-10738-LIVER-LIVER", 
                                                   "GC-WL-10205-SAT-VAT-PBMC-SAT", "GC-WL-10205-SAT-VAT-PBMC-VAT", "GC-WL-10205-SAT-VAT-PBMC-PBMC", "GC-WL-10205-LIVER-LIVER",
                                                   "GC-WL-9991-SAT-VAT-PBMC-SAT", "GC-WL-9991-SAT-VAT-PBMC-VAT", "GC-WL-9991-SAT-VAT-PBMC-PBMC", "GC-WL-9991-LIVER-LIVER",
                                                   "GC-WL-9932-SAT-VAT-PBMC-SAT", "GC-WL-9932-SAT-VAT-PBMC-VAT", "GC-WL-9932-SAT-VAT-PBMC-PBMC", "GC-WL-9932-LIVER-LIVER",
                                                   "GC-WL-11040-SAT-VAT-PBMC-SAT", "GC-WL-11040-SAT-VAT-PBMC-VAT", "GC-WL-11040-SAT-VAT-PBMC-PBMC", "GC-WL-11040-LIVER-LIVER",
                                                   "GC-WL-11051-SAT-VAT-PBMC-SAT", "GC-WL-11051-SAT-VAT-PBMC-VAT", "GC-WL-11051-SAT-VAT-PBMC-PBMC", "GC-WL-11051-LIVER-LIVER",
                                                   "GC-WL-11183-SAT-VAT-PBMC-SAT", "GC-WL-11183-SAT-VAT-PBMC-VAT", "GC-WL-11183-SAT-VAT-PBMC-PBMC","GC-WL-11183-LIVER-LIVER",
                                                   "GC-WL-11471-SAT-VAT-PBMC-SAT", "GC-WL-11471-SAT-VAT-PBMC-VAT", "GC-WL-11471-SAT-VAT-PBMC-PBMC","GC-WL-11471-LIVER-LIVER",
                                                   "GC-WL-11327-SAT-VAT-PBMC-SAT", "GC-WL-11327-SAT-VAT-PBMC-VAT", "GC-WL-11327-SAT-VAT-PBMC-PBMC","GC-WL-11327-LIVER-LIVER"),
                          removeNA = FALSE, 
                          removeMulti = FALSE,
                          filterMulti = FALSE)

combined.TCR = addVariable(combined.TCR, 
                           variable.name = "Type", 
                           variables = c("Healthy", "Healthy", "Healthy", "Healthy",
                                         "No_fibrosis","No_fibrosis","No_fibrosis","No_fibrosis",
                                         "No_fibrosis","No_fibrosis","No_fibrosis","No_fibrosis",
                                         "No_fibrosis","No_fibrosis","No_fibrosis","No_fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis"))

head(combined.TCR[[1]])

# Count total number of T cells across all samples
total_tcells <- sum(sapply(combined.TCR, nrow))
cat("Total number of T cells:", total_tcells, "\n")

# Add Tissue and Sample information to combined TCR
combined.TCR <- combined.TCR %>%
 bind_rows(.id = "Sample") %>%
 mutate(Tissue = str_extract(Sample, "(?<=-)[A-Za-z]+$"))

# Combine all samples into one data frame
all_tcells <- dplyr::bind_rows(combined.TCR, .id = "Sample")

# Count unique clonotypes based on CDR3 AA sequence (clone identifier)
n_unique_clonotypes <- length(unique(all_tcells$CTaa))
cat("Number of unique clonotypes:", n_unique_clonotypes, "\n")

# Count cells per clonotype
clone_counts <- all_tcells %>%
 group_by(CTaa) %>%
 summarise(n_cells = n(), .groups = "drop")

# Identify expanded clonotypes (>1 cell)
hyperexpanded_clones <- clone_counts %>% filter(n_cells > 20)
expanded_clones <- clone_counts %>% filter(n_cells > 1)
singleton_clones <- clone_counts %>% filter(n_cells == 1)

# Calculate summary statistics
n_hyperexpanded_clonotypes <- nrow(hyperexpanded_clones)
n_hyperexpanded_cells <- sum(hyperexpanded_clones$n_cells)

n_expanded_clonotypes <- nrow(expanded_clones)
n_expanded_cells <- sum(expanded_clones$n_cells)

n_singleton_clonotypes <- nrow(singleton_clones)
n_singleton_cells <- sum(singleton_clones$n_cells)

#### Print clone counts ####
cat("HyperExpanded clonotypes:", n_hyperexpanded_clonotypes, "\n")
cat("HyperExpanded T cells:", n_hyperexpanded_cells, "\n")
cat("Expanded clonotypes:", n_expanded_clonotypes, "\n")
cat("Expanded T cells:", n_expanded_cells, "\n")
cat("Singleton clonotypes:", n_singleton_clonotypes, "\n")
cat("Singleton T cells:", n_singleton_cells, "\n")

# Sanity check (expanded + singleton should equal total)
cat("Check total cells (expanded + singleton):", n_expanded_cells + n_singleton_cells, "\n")

# Make sure column names match
colnames(combined.TCR)  # look for cell barcode column (often "cell" or "barcode")

# Count number of cells per clonotype (CTaa)
clone_counts <- combined.TCR %>%
 group_by(CTaa) %>%
 summarise(n_cells = n(), .groups = "drop")

# Define expanded clones
expanded_clones <- clone_counts$CTaa[ clone_counts$n_cells > 1 ]  # >1 cell = expanded

# Extract barcodes of expanded clones
expanded_tcr_barcodes <- combined.TCR$barcode[ combined.TCR$CTaa %in% expanded_clones ]
length(expanded_tcr_barcodes)
# 15564

# Load seurat object
seurat_obj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

DefaultAssay(seurat_obj) <- "RNA"

# Combine all TCR files
tcr_df <- dplyr::bind_rows(combined.TCR)

# Extract cell barcode from TCR barcode
tcr_df$cell_barcode <- sub(".*_", "", tcr_df$barcode)
tcr_df$tissue <- sub(".*-", "", sub("_.*", "", tcr_df$barcode))

head(tcr_df$tissue)

tcr_df$barcode_tissue <- paste(
  tcr_df$tissue,
  tcr_df$cell_barcode,
  sep = "_"
)

head(tcr_df$barcode_tissue)

tcr_meta <- tcr_df %>%
  dplyr::select(
    barcode_tissue,
    CTaa
  ) %>%
  distinct(barcode_tissue, .keep_all = TRUE)

# Prepare Seurat barcodes
seurat_obj@meta.data$cell_barcode <- sub(
  "_\\d+$", "", rownames(seurat_obj@meta.data)
)

seurat_obj@meta.data$tissue_short <- sub(
  " .*", "", seurat_obj@meta.data$Tissue
)

seurat_obj@meta.data$barcode_tissue <- paste(
  seurat_obj@meta.data$tissue_short,
  seurat_obj@meta.data$cell_barcode,
  sep = "_"
)

sum(tcr_df$barcode_tissue %in% seurat_obj@meta.data$barcode_tissue)

intersect(
  tcr_df$barcode_tissue,
  seurat_obj@meta.data$barcode_tissue
)[1:10]

matched_barcodes <- intersect(
  tcr_df$barcode_tissue,
  seurat_obj@meta.data$barcode_tissue
)

length(matched_barcodes)

cells_use <- rownames(seurat_obj@meta.data)[
  seurat_obj@meta.data$barcode_tissue %in% matched_barcodes &
    !duplicated(seurat_obj@meta.data$barcode_tissue)
]

seurat_Tcells <- subset(seurat_obj, cells = cells_use)

ncol(seurat_Tcells)  

# Clean Seurat barcode
seurat_Tcells@meta.data$cell_barcode <- sub("_\\d+$", "", rownames(seurat_Tcells@meta.data))

# Extract short tissue name from Seurat Tissue column
seurat_Tcells@meta.data$tissue_short <- sub(" .*", "", seurat_Tcells@meta.data$Tissue)

# Create Seurat tissue-aware barcode
seurat_Tcells@meta.data$barcode_tissue <- paste(
  seurat_Tcells@meta.data$tissue_short,
  seurat_Tcells@meta.data$cell_barcode,
  sep = "_"
)

tcr_meta <- tcr_meta[
  match(
    seurat_Tcells@meta.data$barcode_tissue,
    tcr_meta$barcode_tissue
  ),
]

seurat_Tcells <- AddMetaData(
  seurat_Tcells,
  metadata = tcr_meta$CTaa,
  col.name = "CTaa"
)

expanded_barcode_tissue <- paste(
  sub(".*-", "", sub("_.*", "", expanded_tcr_barcodes)),  # tissue
  sub(".*_", "", expanded_tcr_barcodes),                  # cell barcode
  sep = "_"
)

seurat_Tcells$Clonality <- ifelse(
  seurat_Tcells@meta.data$barcode_tissue %in% expanded_barcode_tissue,
  "Expanded",
  "Non-expanded"
)

saveRDS(seurat_Tcells, '/data/Blizard-AlazawiLab/rk/seurat/TcellObj.rds')

# Load T cells object
seurat_obj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/TcellObj.rds')

#### Fig2e #####
# Ensure the factor is set correctly
seurat_obj_qc$Clonality <- factor(seurat_obj$Clonality, levels = c("Non-expanded", "Expanded"))

# Set the default assay
DefaultAssay(seurat_obj) <- "RNA"

# Perform differential expression
de_markers <- FindMarkers(
  object = seurat_obj,
  ident.1 = "Expanded",
  ident.2 = "Non-expanded",
  group.by = "Clonality",
  logfc.threshold = 0.25,   
  min.pct = 0.25    
)

# View top genes sorted by log fold change
de_markers <- de_markers[order(de_markers$avg_log2FC, decreasing = TRUE), ]
head(de_markers, n=10)

# Extract significantly upregulated genes (optional)
sig_up_genes <- subset(de_markers, avg_log2FC > 0 & p_val_adj < 0.05)

expanded_genes <- c(
  "FGFBP2",
  "GZMH",
  "GNLY",
  "ADGRG1",
  "FCRL6",
  "GZMB",
  "CX3CR1",
  "S1PR5",
  "NKG7",
  "KLRD1",
  "PRF1",
  "TBX21",
  "PLEK",
  "FCGR3A",
  "ZEB2",
  "CST7",
  "CCL5",
  "EFHD2",
  "CTSW",
  "SPON2",
  "CD8A",
  "CD4"
)

# DotPlot2
ggp <- DotPlot2(
  seurat_obj,
  features = expanded_genes,
  group.by = "Tissue",
  split.by = "Clonality",
  legend_order = c("fill", "color", "size")  # Your preferred order+
) +
  scale_fill_gradient(low = "white", high = "blue") +
  scale_colour_manual(
    values = c(
      "Expanded" = "#084594",     # dark blue
      "Non-expanded" = "#9ecae1"  # light blue
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Expanded vs Non-expanded T cells")

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("ExpandedNonexpandedTcellsDifferentialGenes.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()

#### Fig2.f #####
DefaultAssay(seurat_obj) <- "ADTonly"

# Perform differential expression
de_markers <- FindMarkers(
  object = seurat_obj,
  ident.1 = "Expanded",
  ident.2 = "Non-expanded",
  group.by = "Clonality",
  logfc.threshold = 0.25,   
  min.pct = 0.25    
)

# View top genes sorted by log fold change
de_markers <- de_markers[order(de_markers$avg_log2FC, decreasing = TRUE), ]
print(de_markers)

adt_features <- c(
  "Hu.CD57",
  "Hu.CD8",
  "Hu.KLRG1",
  "Hu.CD11a",
  "Hu.CD18",
  "Hu.CD196",
  "Hu.CD4-RPA.T4",
  "Hu.CD7",
  "Hu.CD27",
  "Hu.CD62L"
)

DefaultAssay(seurat_obj) <- "ADT"

# DotPlot2
ggp <- DotPlot2(
  seurat_obj,
  features = adt_features,
  group.by = "Tissue",
  split.by = "Clonality",
  legend_order = c("fill", "color", "size")  # Your preferred order+
) +
  scale_fill_gradient(low = "white", high = "blue") +
  scale_colour_manual(
    values = c(
      "Expanded" = "#084594",     # dark blue
      "Non-expanded" = "#9ecae1"  # light blue
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Expanded vs Non-expanded T cells")

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("ExpandedNonexpandedTcellsDifferentialADTs.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()

# End of the script
