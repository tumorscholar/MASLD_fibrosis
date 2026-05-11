############################################################
# Script name: Export_Tissue_Expanded_Tcells_to_h5ad.R
#
# Purpose:
#   - Load Seurat object containing expanded T cells
#   - Remove PBMC-derived cells
#   - Retain tissue-expanded T cells only
#   - Export the subset to AnnData (.h5ad) format
#
# Output:
#   TissueExpandedTcells.h5ad
#
############################################################

# Load required libraries
library(Seurat)
library(reticulate)
library(scCustomize)

# Configure Python environment
reticulate::use_python("~/env/py_venv_hpc/bin/python", required = TRUE)

# Confirm packages are available
reticulate::py_module_available("anndata")  # Should return TRUE
reticulate::py_module_available("numpy")    # TRUE
reticulate::py_module_available("scanpy")   # TRUE

# Load expanded T cell
seurat_obj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinal.rds')
table(seurat_obj$expanded_shared_class)

# Subset for tissue expanded T cells
seurat_obj_subset <- subset(
  seurat_obj,
  subset = expanded_shared_class != "a_PBMC"
)

table(seurat_obj_subset$expanded_shared_class)

# Take metadata
meta_df <- seurat_obj_subset@meta.data

# Extract core barcode (AAAC...-1)
core <- sub("_.*$", "", rownames(meta_df))

# Build AnnData-style barcode
meta_df$barcode_match <- paste0(meta_df$Sample, "_", core)

# Put it FIRST column
meta_df <- meta_df[, c("barcode_match", setdiff(colnames(meta_df), "barcode_match"))]

# Write FULL metadata
write.csv(
  meta_df,
  "/data/Blizard-AlazawiLab/rk/python/BarcodesTcellsPython.csv",
  row.names = FALSE,
  quote = FALSE
)

# Export to AnnData (.h5ad)
as.anndata(x = seurat_obj, file_path = "/data/Blizard-AlazawiLab/rk/seurat/", file_name = "TissueExpandedTcells.h5ad")

# End of the script
