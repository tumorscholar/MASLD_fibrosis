#########################################################################
# Script name: 13_Figure3a-c.R
#
# Purpose:
# This script processes TCR contig data using scRepertoire to quantify
# expanded and hyperexpanded clonotypes, summarize clone statistics,
# and generate tissue- and stage-specific Venn diagrams.
#
# Workflow Overview:
# 1) Read TCR contig files and combine using scRepertoire.
# 2) Annotate samples with disease type and tissue labels.
# 3) Compute clone sizes and classify singleton, expanded, and
#    hyperexpanded clonotypes.
# 4) Generate Venn diagrams for expanded and hyperexpanded clonotypes
#    stratified by disease stage.
# 5) Export tabular summaries and save figures for downstream use.
#
# Input:
# - TCR contig CSV files (Cell Ranger output).
#
# Output:
# - SVG Venn diagrams for Expanded and Hyperexpanded clonotypes.
# - CSV tables summarizing clone counts and metadata.
#
#
##########################################################################

library(Seurat)
library(scRepertoire)
library(tidyverse)      
library(ggVennDiagram)  
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

#### Fig3a ####
# Venn diagram expanded clonetypes (3098)
# Combine all list elements into one dataframe
combined_TCR_df <- bind_rows(combined.TCR, .id = "sample_id")

# Check
dim(combined_TCR_df)
head(combined_TCR_df)
length(unique(combined_TCR_df$sample_id))
table(combined_TCR_df$sample_id)

# Ensure CTaa column exists in both dataframes
stopifnot("CTaa" %in% colnames(expanded_clones))
stopifnot("CTaa" %in% colnames(combined_TCR_df))

# Join to get tissue info for expanded clonotypes
expanded_combined <- combined_TCR_df %>%
  dplyr::select(CTaa, Tissue, Type) %>%
  distinct() %>%
  inner_join(expanded_clones, by = "CTaa") %>%
  group_by(CTaa) %>%
  summarise(
    tissues = paste(unique(Tissue), collapse = ","),
    Type = paste(unique(Type), collapse = ","),
    n_cells = sum(n_cells, na.rm = TRUE),
    .groups = "drop"
  )

head(expanded_combined)

expanded_tissue_summary <- expanded_combined %>%
  separate_rows(tissues, sep = ",") %>%
  group_by(tissues) %>%
  summarise(n_expanded_clonotypes = n_distinct(CTaa), .groups = "drop")

expanded_tissue_summary

# Expand tissue column so each clonotype has one row per tissue
expanded_clonotypes_long <- expanded_combined %>%
  separate_rows(tissues, sep = ",") %>%
  distinct(CTaa, tissues)

# Check total unique expanded clonotypes
length(unique(expanded_clonotypes_long$CTaa))  # should return 3098

# Optional sanity check: sum per tissue
expanded_tissue_summary <- expanded_clonotypes_long %>%
  group_by(tissues) %>%
  summarise(n_expanded_clonotypes = n_distinct(CTaa), .groups = "drop")

sum(expanded_tissue_summary$n_expanded_clonotypes)  # returns 4989

# Total unique expanded clonotypes across all tissues
length(unique(expanded_clonotypes_long$CTaa))  # should return 3098

# Sum per tissue (will be >3098 due to overlaps)
sum(expanded_tissue_summary$n_expanded_clonotypes)  # returns 4989

# Prepare long-format clonotype table
expanded_long <- expanded_combined %>%
  rename(Stage = Type) %>%  # Rename to "Stage" for clarity
  separate_rows(Stage, sep = ",") %>%
  mutate(Stage = trimws(Stage)) %>%
  separate_rows(tissues, sep = ",") %>%
  mutate(tissues = trimws(tissues)) %>%
  distinct(CTaa, Stage, tissues) %>%
  filter(tissues != "Unknown")

# Define Venn plotting function
plot_stage_venn <- function(stage_label, color = "#1f78b4") {
  
  # Filter the stage-specific clonotypes
  stage_data <- expanded_long %>%
    filter(Stage == stage_label)
  
  if (nrow(stage_data) == 0) {
    message("⚠️ No expanded clonotypes found for stage: ", stage_label)
    return(NULL)
  }
  
  # Count unique clonotypes for title
  n_stage_clonotypes <- stage_data %>%
    distinct(CTaa) %>%
    nrow()
  
  # Build tissue-wise list of clonotypes
  tissue_lists <- stage_data %>%
    group_by(tissues) %>%
    summarise(clones = list(unique(CTaa)), .groups = "drop") %>%
    deframe()
  
  # Print per-tissue summary
  tissue_counts <- stage_data %>%
    group_by(tissues) %>%
    summarise(n_clones = n_distinct(CTaa), .groups = "drop")
  print(tissue_counts)
  
  # Plot Venn diagram
  ggp <- ggVennDiagram(
    tissue_lists,
    label_alpha = 0,
    edge_size = 1.2,
    set_size = 6
  ) +
    scale_fill_gradient(low = "#FFFFFF", high = color) +
    theme(
      text = element_text(family = "Arial", face = "bold", size = 14),
      plot.background = element_rect(fill = "white"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = paste0(
      "Expanded clonotypes across – ", stage_label,
      " (n = ", n_stage_clonotypes, ")"
    ))
  
  print(ggp)
  return(ggp)
}

# Generate Venn diagrams per stage
venn_fib     <- plot_stage_venn("Fibrosis", color = "#e41a1c")
venn_no_fib  <- plot_stage_venn("No_fibrosis", color = "#377eb8")
venn_healthy <- plot_stage_venn("Healthy", color = "#4daf4a")

# Save plots
out_dir <- "/data/home/hdx044/plots/screpertoire"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "Venn_Fibrosis.svg"), plot = venn_fib, width = 6, height = 6)
ggsave(file.path(out_dir, "Venn_No_Fibrosis.svg"), plot = venn_no_fib, width = 6, height = 6)
ggsave(file.path(out_dir, "Venn_Healthy.svg"), plot = venn_healthy, width = 6, height = 6)

cat("✅ Saved 3 SVG plots in:", out_dir, "\n")


#### Fig3b ####
# Count number of cells per clonotype (clone size)
clone_counts <- all_tcells %>%
  group_by(CTaa) %>%
  summarise(n_cells = n(), .groups = "drop")

# Add clone size (n_cells) back to each cell
all_tcells <- all_tcells %>%
  left_join(clone_counts, by = "CTaa")

# Inspect barcode, clonotype, and clone size
head(all_tcells[, c("barcode", "CTaa", "n_cells")])

# Subset cells belonging to hyperexpanded clonotypes (>20 cells)
hyperexpanded_cells <- all_tcells %>%
  filter(n_cells > 20)

# Extract unique hyperexpanded clonotypes (one row per CTaa)
unique_CTaa <- hyperexpanded_cells %>%
  distinct(CTaa)

# View unique hyperexpanded clonotypes
unique_CTaa

# Summarise metadata for each hyperexpanded clonotype
CTaa_metadata <- hyperexpanded_cells %>%
  group_by(CTaa) %>%
  summarise(
    clone_size        = first(n_cells),                          # total clone size
    n_cells           = n(),                                      # number of rows/cells
    n_patients        = n_distinct(Patient_ID),                  # number of patients
    n_tissues         = n_distinct(Tissue),                      # number of tissues
    dominant_tissue   = names(sort(table(Tissue), decreasing = TRUE))[1],
    dominant_patient  = names(sort(table(Patient_ID), decreasing = TRUE))[1],
    Type              = first(Type),                              # disease type
    .groups = "drop"
  )

# Attach metadata to unique hyperexpanded clonotypes
unique_CTaa <- unique_CTaa %>%
  left_join(CTaa_metadata, by = "CTaa")

# Inspect final hyperexpanded clonotype table
unique_CTaa

# Export hyperexpanded clonotype metadata table
write_csv(
  unique_CTaa,
  "/data/home/hdx044/files/screpertoire/HyperexpandedClonesDetails.csv"
)

# Data was plotted in Prism.

#### Fig3c ####
# Venn diagram hyperexpanded clonotypes (95)
# Prepare long-format table for hyperexpanded clones
hyperexpanded_clones <- clone_counts %>% filter(n_cells > 20)
hyperexpanded_long <- expanded_long %>%
  filter(CTaa %in% hyperexpanded_clones$CTaa) %>%
  select(CTaa, Stage, tissues) %>%
  distinct()  # one row per clone per tissue per stage


# Venn plotting function per stage

plot_stage_venn <- function(stage_label, color="#1f78b4") {
  
  # Filter for stage
  stage_data <- hyperexpanded_long %>%
    filter(Stage == stage_label)
  
  if (nrow(stage_data) == 0) {
    message("⚠️ No hyperexpanded clones found for stage: ", stage_label)
    return(NULL)
  }
  
  # Count unique clones per tissue for plotting
  tissue_lists <- stage_data %>%
    group_by(tissues) %>%
    summarise(clones = list(unique(CTaa)), .groups = "drop") %>%
    deframe()
  
  n_stage_clonotypes <- length(unique(stage_data$CTaa))
  
  cat("\nStage:", stage_label, "(n hyperexpanded clones =", n_stage_clonotypes, ")\n")
  print(stage_data %>%
          group_by(tissues) %>%
          summarise(n_clones = n_distinct(CTaa)))
  
  # Plot Venn diagram
  ggp <- ggVennDiagram(
    tissue_lists,
    label_alpha = 0,
    edge_size = 1.2,
    set_size = 6
  ) +
    scale_fill_gradient(low = "#FFFFFF", high = color) +
    theme(
      text = element_text(family = "Arial", face = "bold", size = 14),
      plot.background = element_rect(fill = "white"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = paste0(
      "Hyperexpanded clonotypes – ",
      stage_label,
      " (n=", n_stage_clonotypes, ")"
    ))
  
  print(ggp)
  return(ggp)
}


# Plot Venn diagrams for each stage

venn_fib     <- plot_stage_venn("Fibrosis", color="#e41a1c")
venn_no_fib  <- plot_stage_venn("No_fibrosis", color="#377eb8")

# Step 5: Save plots
out_dir <- "/data/home/hdx044/plots/screpertoire"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "Hyperexpanded_Venn_Fibrosis.svg"), plot = venn_fib, width = 6, height = 6)
ggsave(file.path(out_dir, "Hyperexpanded_Venn_No_Fibrosis.svg"), plot = venn_no_fib, width = 6, height = 6)

cat("✅ Saved 3 SVG plots in:", out_dir, "\n")

# End of the script