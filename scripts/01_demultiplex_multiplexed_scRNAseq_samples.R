
###############################################################################
# Script name: 01_demultiplex_multiplexed_scRNAseq_samples.R
#
# Description:
# This script is designed to demultiplex multiplexed single‑cell RNA‑seq data
# generated using 10x Genomics technology and antibody‑based hashing.
# The primary purpose of this script is to assign correct tissue identities
# (SAT, VAT, or PBMC) to individual cells and store this information in the
# Seurat object metadata.
#
# The input data consist of Cell Ranger output files
# (filtered_feature_bc_matrix.h5) containing gene expression data and antibody
# capture counts. Antibody‑derived tags (ADT) are used to perform hashtag
# oligonucleotide (HTO)‑based demultiplexing of pooled samples.
#
# Key steps performed in this script include:
#  - Loading Cell Ranger HDF5 output files
#  - Creating Seurat objects for each sample
#  - Storing antibody capture data as a separate ADT assay
#  - Separating liver samples from multiplexed SAT‑VAT‑PBMC samples
#  - Normalizing HTO counts and performing demultiplexing using HTODemux
#  - Retaining only singlet cells from multiplexed samples
#  - Assigning tissue labels (SAT, VAT, PBMC) based on HTO classification
#  - Adding tissue and patient identifiers to Seurat metadata
#  - Merging all tissues into a single Seurat object
#
# Importantly, this script does NOT perform quality control (QC), filtering,
# normalization of gene expression data, or downstream analysis.
# All cells passing demultiplexing are retained.
# Quality control and downstream analyses are handled in separate scripts.
#
# Outputs:
#  - Seurat objects with tissue metadata assigned at the cell level
#  - Barcode‑level metadata tables linking cells to sample and tissue identity
#
# Purpose:
# To generate reliable tissue annotations for multiplexed scRNA‑seq data that
# can be used in downstream QC and analysis pipelines.

###############################################################################

# Load packages
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)

# Locate files 
dataDir <- '/data/Blizard-AlazawiLab/rk/cellranger'
files <- list.files(path = dataDir,
                    pattern = 'filtered_feature_bc_matrix.h5',
                    recursive = T,
                    full.names = T)

# Remove the files which are not multiplexed and not included
# Make sure to load only 40 files
files <- files[-c(27, 28)]

# Extract sample names from file names
pattern1 <- '/data/Blizard-AlazawiLab/rk/cellranger/'
samples <- gsub(pattern1,
                '',
                files)

pattern2 <- '/outs/filtered_feature_bc_matrix.h5'
samples <- gsub(pattern2,
                '',
                samples)

# Load all data 
rawDat <- lapply(files, function(x){
 Read10X_h5(x)
}) 

# Specify list element names 
names(rawDat) <- samples

# Make seurat objects in for loop
# Initialize an empty list to store Seurat objects
datObjs <- list()

# Loop through each sample
for (x in samples) {
 # Create Seurat object
 datObjs[[x]] <- CreateSeuratObject(counts = rawDat[[x]][['Gene Expression']], project = x)
 
 # Add antibody data to a new assay
 datObjs[[x]][['ADT']] <- CreateAssayObject(counts = rawDat[[x]][['Antibody Capture']])
 
 # Add sample information to metadata
 datObjs[[x]]@meta.data$sample <- x
 
}

setwd("~/seurat")
saveRDS(datObjs, 'SeuObjContig.rds')

# Load data list
mylist <- readRDS("~/seurat/SeuObjContig.rds")

# Split list into LIVER and SVP objects
names(mylist)

# Subset LIVER
grep('LIVER', names(mylist))
LIVERlist <- mylist[grep('LIVER', names(mylist))]

# Subset SVP
grep('SAT-VAT-PBMC', names(mylist))
SVPlist <- mylist[grep('SAT-VAT-PBMC', names(mylist))]


#### Add LIVER metadata ####
# Name LIVER list
names(LIVERlist)
LIVER <- names(LIVERlist)[1]
LIVERlist[[LIVER]]

for (x in names(LIVERlist)) {
  
  # Create metadata columns for Tissue type and Patient_ID 
  LIVERlist[[x]]$Tissue <- 'LIVER'
  id <- gsub('-LIVER', '', x)
  id <- gsub('GC-WL-', '', id)
  LIVERlist[[x]]$Patient_ID <- id
  
}

#### Add SVP metadata ####
# Name SVP list
names(SVPlist)
SVP <- names(SVPlist)[1]
SVPlist[[SVP]]

# Check ADT list
rownames(SVPlist$`GC-WL-10113-1-SAT-VAT-PBMC`@assays$ADT$counts)

for (x in names(SVPlist)) {
  # Split ADT and HTO data into separate assays
  SVPlist[[x]][['HTO']] <- CreateAssayObject(counts = SVPlist[[x]]@assays$ADT$counts[c('C0251', 'C0252', 'C0253'), ])
  SVPlist[[x]][['ADTonly']] <- CreateAssayObject(counts = SVPlist[[x]]@assays$ADT$counts[c(1:130), ])
  
  # Normalise HTO data
  SVPlist[[x]] <- NormalizeData(SVPlist[[x]], assay = "HTO", normalization.method = "CLR", margin = 2)
  
  # Demultiplex cells based on HTO enrichment
  SVPlist[[x]] <- HTODemux(SVPlist[[x]], assay = "HTO", positive.quantile = 0.95)
  
  # Keep only singlets
  SVPlist[[x]] <- subset(SVPlist[[x]], subset = HTO_classification.global == 'Singlet')
  
  # Demux results are saved in HTO classification column under meta data
  # Create tissue vector using HTO classification
  tissue_vector <- c(rep(NA, length(SVPlist[[x]]$HTO_classification)))
  tissue_vector[SVPlist[[x]]$HTO_classification == 'C0251'] <- 'SAT'
  tissue_vector[SVPlist[[x]]$HTO_classification == 'C0252'] <- 'VAT'
  tissue_vector[SVPlist[[x]]$HTO_classification == 'C0253'] <- 'PBMC'
  
  # Create metadata columns for Tissue type and Patient_ID and add tissue vector in tissue type
  SVPlist[[x]]$Tissue <- tissue_vector
  id <- gsub('-SAT-VAT-PBMC', '', x)
  id <- gsub('GC-WL-', '', id)
  SVPlist[[x]]$Patient_ID <- id
  
}

# Merge the liver and SAT, VAT and PBMC list
all.tissue.list <- c(LIVERlist, SVPlist)
SeuObjForDemux <- merge(x = all.tissue.list[[1]], 
                        y = all.tissue.list[2:40])

# Extract metadata
metadata <- SeuObjForDemux@meta.data
write.csv(metadata, "/data/Blizard-AlazawiLab/rk/seurat/demuxBarcodes.csv")

saveRDS(SeuObjForDemux, 'SeuObjForDemux.rds')

# Load obj
SeuObjForDemux <- readRDS("~/seurat/SeuObjForDemux.rds")

Idents(SeuObjForDemux) <- "sample"

# Set working directory
setwd("/data/Blizard-AlazawiLab/rk/barcodes")

# Subset cells where sample name includes "-SAT-VAT-PBMC"
SeuObj_SAT_VAT_PBMC <- subset(
  SeuObjForDemux,
  subset = sample %in% grep("-SAT-VAT-PBMC$", SeuObjForDemux$sample, value = TRUE)
)

table(SeuObj_SAT_VAT_PBMC$sample)
head(SeuObj_SAT_VAT_PBMC)

# Extract metadata
meta <- SeuObj_SAT_VAT_PBMC@meta.data

# Ensure cell names are included
meta$cell_barcode <- rownames(meta)

# Create output directory (optional)
output_dir <- "/data/Blizard-AlazawiLab/rk/barcodes/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Split by sample and write each as CSV
samples <- unique(meta$sample)
for (s in samples) {
  meta_subset <- meta[meta$sample == s, c("cell_barcode", "sample", "Tissue")]
  out_file <- file.path(output_dir, paste0(s, "_metadata.csv"))
  write.csv(meta_subset, file = out_file, row.names = FALSE)
}

# Set the path where metadata files are saved
input_dir <- "/data/Blizard-AlazawiLab/rk/barcodes/"

# List all CSV files in the directory
metadata_files <- list.files(input_dir, pattern = "_metadata\\.csv$", full.names = TRUE)

# Read and combine all metadata files
metadata_list <- lapply(metadata_files, read.csv)

# Combine into one data frame
combined_metadata <- do.call(rbind, metadata_list)

# Check the result
head(combined_metadata)
tail(combined_metadata)

# Remove sample suffix
combined_metadata$cell_barcode <- sub("_[2-4][0-9]$", "", combined_metadata$cell_barcode)
head(combined_metadata)
tail(combined_metadata)

# Check unique samples
unique_samples <- unique(combined_metadata$sample)

# Create an output directory
output_dir <- "/data/Blizard-AlazawiLab/rk/barcodes"
dir.create(output_dir, showWarnings = FALSE)

# Split and save per sample
for (samp in unique_samples) {
  sample_df <- subset(combined_metadata, sample == samp)
  
  output_file <- file.path(output_dir, paste0(samp, "_metadata.csv"))
  
  write.csv(sample_df[, c("cell_barcode", "sample", "Tissue")],
            file = output_file,
            row.names = FALSE)
  
  cat("Saved:", output_file, "\n")
}

# End of the script
