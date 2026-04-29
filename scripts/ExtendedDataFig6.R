##############################################################
# Title: Extended Data Figure 6.R
# 
# This script performs:
#  - Identification of non-expanded CD4+ and CD8+ T cells
#  - Quantification of CD4/CD8 proportions across tissues
#  - Assessment of tissue residency marker expression
#  - Reactome pathway enrichment of non-expanded T cells
#    using FGSEA across liver, PBMC, SAT, and VAT
#
# Input:
#  - Seurat object containing T cell scRNA-seq data
#
# Output:
#  plots and CSV tables for Extended Data
#
# Reproducibility:
#  - Fixed random seeds are used for FGSEA analyses
############################################################

library(Seurat)
library(tidyverse)
library(msigdbr)
library(fgsea)
library(viridis)
library(Matrix)

# Load T cells object
SeuObj  <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/TcellObj.rds')

# Change stage name writing style
SeuObj$Stage <- gsub(
  "No_fibrosis",
  "No_Fibrosis",
  SeuObj$Stage
)

# Extract expression values
cd4_expr <- FetchData(SeuObj, vars = "CD4", layer = "data")
cd8_expr <- FetchData(SeuObj, vars = c("CD8A", "CD8B"), layer = "data")

# Gate cells
SeuObj$T_cell_type <- case_when(
  cd4_expr$CD4 > 0.25 ~ "CD4",
  (cd8_expr$CD8A > 0.25 | cd8_expr$CD8B > 0.25) ~ "CD8",
  TRUE ~ "Other"
)

# Subset Non-expanded + CD4 or CD8 only
meta <- SeuObj@meta.data %>%
  filter(Clonality == "Non-expanded",
         T_cell_type %in% c("CD4", "CD8"))

# Calculate percentage per Tissue + Stage + T_cell_type
plot_data <- meta %>%
  group_by(Tissue, Stage, T_cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Tissue, Stage) %>%
  mutate(pct = n / sum(n) * 100)

# Set factor levels for correct ordering
plot_data$Stage <- factor(plot_data$Stage, levels = c("Healthy", "No_Fibrosis", "Fibrosis"))
plot_data$Tissue <- factor(plot_data$Tissue, levels = c("LIVER", "PBMC", "SAT", "VAT"))

#### Extended data Fig 6a ####
ggp <- ggplot(plot_data, aes(x = Tissue, y = pct, fill = T_cell_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Stage ~ ., switch = "y") +          # switch="y" moves strip to the left
  scale_fill_manual(values = c("CD4" = "#4472C4", "CD8" = "#E41A1C"),
                    name = "T Cell Type") +
  scale_y_continuous(name = "Percentage (%)") +  # adds y-axis label as "Percentage (%)"
  theme_bw() +
  theme(
    strip.text.y.left = element_text(angle = 90, size = 11),  # stage label on LEFT side, rotated
    strip.placement = "outside",                               # moves strip outside the axis
    strip.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12),
    panel.border = element_rect(color = "black"),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("NonExpandedTcellsCD4CD8.svg", width = 4, height = 4)   # width/height in cm by default for svg()
print(ggp)
dev.off()

#### Extended data Fig 6b ####

trm_genes <- c(
  "NR4A1",
  "BHLHE40",
  "RUNX3",
  "KLF2",
  "S1PR1",
  "SELL",
  "CD44",
  "ITGB1",
  "CD101",
  "CXCR6",
  "ITGA1",
  "ITGAE",
  "CD69"
)

# Make tissue stage column
SeuObj$Tissue.Stage <- paste(SeuObj$Tissue, SeuObj$Stage, sep = ".")

# Reset label
SeuObj$Tissue.Stage <- factor(
  SeuObj$Tissue.Stage,
  levels = c(
    "LIVER.Healthy", "LIVER.No_Fibrosis", "LIVER.Fibrosis",
    "PBMC.Healthy",  "PBMC.No_Fibrosis",  "PBMC.Fibrosis",
    "SAT.Healthy",   "SAT.No_Fibrosis",   "SAT.Fibrosis",
    "VAT.Healthy",   "VAT.No_Fibrosis",   "VAT.Fibrosis"
  )
)


ggp <- DotPlot2(SeuObj, features = trm_genes, group.by = "Tissue.Stage")

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("NonExpandedTcellsTRM.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()


#### Extended data Fig 6c ####

counts <- GetAssayData(
  SeuObj,
  assay = "RNA",
  layer = "counts"
)

# Subset: No Fibrosis + Fibrosis
cells_keep <- rownames(SeuObj@meta.data)[
  SeuObj$Stage %in% c("No_Fibrosis", "Fibrosis")
]

counts_sub <- counts[, cells_keep, drop = FALSE]

SeuObj_sub <- CreateSeuratObject(
  counts    = counts_sub,
  meta.data = SeuObj@meta.data[cells_keep, ]
)

# Subset further: Non‑expanded T cells
SeuObj_nonexp <- subset(
  SeuObj_sub,
  subset = Clonality == "Non-expanded"
)

table(SeuObj_nonexp$Clonality)

# Sanity check
table(SeuObj_nonexp$Stage)
table(SeuObj_nonexp$Clonality)

# Load Reactome pathways
reactome_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:REACTOME"
) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Define tissues
tissues <- c("LIVER", "PBMC", "SAT", "VAT")
fgsea_all <- list()

# Loop over tissues
for (tissue in tissues) {
  
  message("Processing: ", tissue)
  
  # Reproducibility
  set.seed(123 + which(tissues == tissue))
  
  # Cells for this tissue
  cells_tissue <- rownames(SeuObj_nonexp@meta.data)[
    grepl(paste0("^", tissue), SeuObj_nonexp$Tissue)
  ]
  
  if (length(cells_tissue) < 50) next
  
  obj_t <- subset(SeuObj_nonexp, cells = cells_tissue)
  
  counts <- GetAssayData(obj_t, assay = "RNA", slot = "counts")
  
  # Average expression across cells
  avg_expr <- Matrix::rowMeans(counts)
  
  # Rank genes
  ranks <- sort(avg_expr, decreasing = TRUE)
  
  # Jitter ties (reproducible)
  ranks <- ranks + runif(length(ranks), -1e-6, 1e-6)
  ranks <- sort(ranks, decreasing = TRUE)
  
  # FGSEA
  fg <- fgsea(
    pathways    = reactome_sets,
    stats       = ranks,
    minSize     = 15,
    maxSize     = 500,
    nPermSimple = 10000
  )
  
  fg$Tissue <- tissue
  fgsea_all[[tissue]] <- fg
}

# Combine results
fgsea_all <- bind_rows(fgsea_all)

top_pathway_per_tissue <- fgsea_all %>%
  dplyr::filter(!is.na(padj)) %>%       # drop NA pathways
  dplyr::group_by(Tissue) %>%
  dplyr::slice_min(order_by = padj, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

print(top_pathway_per_tissue, n = 40)

results_csv <- top_pathway_per_tissue %>%
  dplyr::select(
    Tissue,
    pathway,
    NES,
    padj,
    pval,
    size
  )

write.csv(
  results_csv,
  "/data/home/hdx044/files/fgsea/top10_pathways_nonexpanded_Tcells.csv",
  row.names = FALSE
)

selected_pathways <- c(
  "REACTOME_TRANSLATION",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "REACTOME_RRNA_PROCESSING",
  "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
  "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
  "REACTOME_SELENOAMINO_ACID_METABOLISM",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
  "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
  "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",
  "REACTOME_INFLUENZA_INFECTION",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  "REACTOME_FORMATION_OF_THE_40S_RIBOSOMAL_SUBUNIT",
  "REACTOME_FORMATION_OF_THE_60S_RIBOSOMAL_SUBUNIT",
  "REACTOME_RIBOSOME_BIOGENESIS",
  "REACTOME_MRNA_SPLICING",
  "REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_MRNA",
  "REACTOME_TRNA_AMINOACYLATION",
  "REACTOME_AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE",
  "REACTOME_CELLULAR_RESPONSES_TO_EXTERNAL_STIMULI",
  "REACTOME_PROTEIN_STABILIZATION"
)

# Subset FGSEA results 
plot_df <- fgsea_all %>%
  dplyr::filter(
    pathway %in% selected_pathways,
    !is.na(padj)
  ) %>%
  mutate(
    pathway = factor(pathway, levels = selected_pathways)
  )

table(plot_df$Tissue)

# Clean pathway names for the figure
plot_df$pathway <- gsub("REACTOME_", "", plot_df$pathway)
plot

p <- ggplot(
  plot_df,
  aes(
    x = Tissue,
    y = pathway,
    size = -log10(padj),
    color = NES
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_viridis(
    option = "viridis",
    direction = 1
  ) +
  scale_size(range = c(2, 6)) +
  labs(
    x = "Tissue",
    y = "Reactome pathway",
    color = "NES",
    size = expression(-log[10]~adj.~p)
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p

ggsave(
  filename = "/data/home/hdx044/plots/fgsea/Nonexpanded_Tcells_translation_dotplot.png",
  plot = p,
  width = 8,
  height = 5,
  units = "in"
)

# End of the script
