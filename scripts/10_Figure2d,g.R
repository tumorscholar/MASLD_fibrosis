###############################################################################
# Script name: 10_Figure2d,g.R
#
# Description:
# This script performs multi‑tissue T‑cell receptor (TCR) repertoire analysis
# using scRepertoire and integrates clonotype information with a single‑cell
# Seurat object. The analysis includes SAT, VAT, PBMC, and liver samples across
# Healthy, No Fibrosis, and Fibrosis conditions.
#
# The script reads per‑sample TCR contig annotation files generated from
# single‑cell V(D)J sequencing, combines them into a unified TCR repertoire
# object, and annotates samples with disease status and tissue identity.
#
# Key steps performed in this script include:
#  - Loading TCR contig annotation CSV files from multiple tissues and patients
#  - Combining TCR repertoires across all samples using scRepertoire::combineTCR
#  - Assigning disease status (Healthy, No_fibrosis, Fibrosis) to each sample
#  - Assigning tissue labels (SAT, VAT, PBMC, LIVER) to individual T cells
#  - Calculating total T‑cell counts and unique TCR clonotypes
#  - Exporting unique TRA/TRB CDR3 amino‑acid clonotypes for external annotation
#  - Mapping TCR clonotypes back to a Seurat single‑cell object via cell barcodes
#  - Quantifying clonal expansion and clonal occupancy across tissues and clusters
#
# Clonal expansion states are defined using scRepertoire and used to characterise
# singleton and expanded TCR clonotypes across tissues.
# Exported tables and figures are used for downstream visualisation and analysis.
#
# Outputs:
#  - Combined TCR repertoire object
#  - CSV file of unique TRA/TRB clonotype sequences
#  - CSV files describing expanded clonotypes and clonal occupancy
#
# Purpose:
# To characterise TCR clonotype expansion, and tissue distribution,
# and to integrate TCR repertoire information with single‑cell transcriptomic
# data in MASLD.
#

###############################################################################

# Load packages
library(tidyverse)     
library(scRepertoire)  
library(Seurat)        
library(patchwork)     

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


# 1️⃣ Define tissue order per patient
tissues <- c("SAT", "VAT", "PBMC", "LIVER")

# 2️⃣ Repeat for 20 patients
tissue_vector <- rep(tissues, times = 20)

# 3️⃣ Check length matches number of samples
length(tissue_vector)   # 80
length(combined.TCR)    # 80, matches ✅

# 4️⃣ Assign tissue info to each sample in the list
for (i in seq_along(combined.TCR)) {
 combined.TCR[[i]]$Tissue <- tissue_vector[i]
}

head(combined.TCR[[1]])

# Count total number of T cells across all samples
total_tcells <- sum(sapply(combined.TCR, nrow))
cat("Total number of T cells:", total_tcells, "\n")
# Total number of T cells: 49363 

# Combine all samples into one data frame
all_tcells <- dplyr::bind_rows(combined.TCR, .id = "Sample")

# Count unique clonotypes based on CDR3 AA sequence (clone identifier)
n_unique_clonotypes <- length(unique(all_tcells$CTaa))
cat("Number of unique clonotypes:", n_unique_clonotypes, "\n")
# 36897

unique_TRA_TRB_CTaa <- all_tcells %>%
  dplyr::filter(!is.na(CTaa), CTaa != "") %>%
  dplyr::distinct(CTaa)

nrow(unique_TRA_TRB_CTaa)

# Write to CSV
write_csv(unique_TRA_TRB_CTaa, "/data/home/hdx044/files/screpertoire/unique_TRA_TRB_CTaa.csv")

#### Fig2.d ####
# sequence from unique_TRA_TRB_CTaa.csv were annotated at TCRdb2.0 and TCRMatch of Immune Epitope Database Analysis Resource (IEDB) and plotted in PRISM.

#### Map VDJ seq identified T cell in single cell obj ####

# Load seurat object
seurat_obj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

# Extract and manipulate barcodes names
# Extract metadata
metadata <- seurat_obj@meta.data

# Manipulate row names (remove the numeric suffix after the underscore)
change_barcodes <- gsub("(_[0-9]+)", "", rownames(metadata))
new_barcodes <- paste0(seurat_obj$Sample, '-', seurat_obj$Tissue, '_', change_barcodes)

# Ensure that the new barcodes are unique
if (length(unique(new_barcodes)) != length(new_barcodes)) {
  stop("The new barcode IDs are not unique after manipulation!")
}

# Ensure consistency between metadata row names and Seurat object cell names
if (all(rownames(metadata) == colnames(seurat_obj))) {
  # Update the row names of the metadata
  rownames(metadata) <- new_barcodes
  
  # Update the cell names in the Seurat object
  colnames(seurat_obj) <- new_barcodes
  
  # Assign the manipulated metadata back to the Seurat object
  seurat_obj@meta.data <- metadata
} else {
  stop("Metadata row names and Seurat object cell names do not match!")
}

# Verify the changes
head(seurat_obj@meta.data)
tail(seurat_obj@meta.data)

write.csv(seurat_obj@meta.data, file="/data/home/hdx044/files/seurat/allTissue/renamed_barcode.csv")

#### Add T cells barcodes to single cell obj ####

# Load the CSV file with barcode and T cell information
metadata_df <- read.csv("/data/home/hdx044/files/seurat/allTissue/old/tcell_barcode.csv", header = TRUE)

# Ensure the barcodes are character
metadata_df$barcode <- as.character(metadata_df$barcode)

# Set the barcodes as row names
rownames(metadata_df) <- metadata_df$barcode

# Remove the barcode column as it's now the row names
metadata_df$barcode <- NULL

# Ensure the Seurat object and metadata data frame have matching barcodes
common_barcodes <- intersect(rownames(metadata_df), colnames(seurat_obj))
metadata_df <- metadata_df[common_barcodes, , drop = FALSE]

# Add the metadata to the Seurat object
seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata_df)

# Verify the added metadata
head(seurat_obj@meta.data)

#### Calculating clone Size all tissue ####
Idents(seurat_obj) <- "Tissue"

seurat_obj <- combineExpression(combined.TCR, 
                                seurat_obj, 
                                cloneCall="gene", 
                                group.by = "Type", 
                                proportion = TRUE)

# Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

seurat_obj <- combineExpression(combined.TCR, 
                                seurat_obj, 
                                cloneCall="gene", 
                                group.by = "Type", 
                                proportion = FALSE, 
                                cloneSize=c(Single=1, Small=2, Medium=6, Large=11, Hyperexpanded=21))

# clonal occupancy
clonalOccupy(seurat_obj, 
             x.axis = "seurat_clusters", exportTable = NULL)

#### Fig2.g #### 
metadata <- clonalOccupy(seurat_obj, 
             x.axis = "seurat_clusters", exportTable = T)

write.csv(metadata, "/data/home/hdx044/files/screpertoire/expandedCloneDetails.csv")

# exported metadata was plotted in Prism.

# End of the script
