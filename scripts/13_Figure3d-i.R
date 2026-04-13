###############################################################################
# Title: 13_Figure3-d-i.R
#
# Description:
#   Analysis of single-cell TCR repertoires across LIVER, PBMC, SAT, and VAT
#   tissues, comparing Healthy, No Fibrosis, and Fibrosis conditions.
#   The script identifies expanded clonotypes, classifies tissue-sharing
#   patterns, integrates clonality with a Seurat T-cell object, and performs
#   differential gene and ADT analyses with downstream visualisation.
#
# Key outputs:
#   - Expanded and hyperexpanded clonotype statistics
#   - Tissue- and stage-specific expansion profiles
#   - Seurat-based DotPlots for migration, activation, effector, and fibrosis
#     associated gene and ADT signatures
#
###############################################################################

library(Seurat)
library(scRepertoire)
library(SeuratExtend)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(patchwork)
library(ggupset)
library(ComplexUpset)
library(circlize)

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


# ️Define tissue order per patient
tissues <- c("SAT", "VAT", "PBMC", "LIVER")

# Repeat for 20 patients
tissue_vector <- rep(tissues, times = 20)

# Check length matches number of samples
length(tissue_vector)   # 80
length(combined.TCR)    # 80, matches ✅

# ️Assign tissue info to each sample in the list
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

#### Percent of expanded clones in tissue and stage ###

# Calculate expanded cells breakdown by Tissue
expanded_by_tissue <- all_tcells %>%
  mutate(is_expanded = CTaa %in% expanded_clones$CTaa) %>%
  group_by(Tissue) %>%
  summarise(
    total_cells = n(),
    expanded_cells = sum(is_expanded),
    percent_expanded = (expanded_cells / total_cells) * 100,
    contribution_to_total_expanded = (expanded_cells / n_expanded_cells) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(contribution_to_total_expanded))

cat("\n=== EXPANDED CELLS BY TISSUE ===\n")
print(expanded_by_tissue)

# Calculate expanded cells breakdown by Stage (Type)
expanded_by_stage <- all_tcells %>%
  mutate(is_expanded = CTaa %in% expanded_clones$CTaa) %>%
  group_by(Type) %>%
  summarise(
    total_cells = n(),
    expanded_cells = sum(is_expanded),
    percent_expanded = (expanded_cells / total_cells) * 100,
    contribution_to_total_expanded = (expanded_cells / n_expanded_cells) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(contribution_to_total_expanded))

cat("\n=== EXPANDED CELLS BY STAGE (Type) ===\n")
print(expanded_by_stage)

# Combined breakdown: Tissue × Stage
expanded_by_tissue_stage <- all_tcells %>%
  mutate(is_expanded = CTaa %in% expanded_clones$CTaa) %>%
  group_by(Tissue, Type) %>%
  summarise(
    total_cells = n(),
    expanded_cells = sum(is_expanded),
    percent_expanded = (expanded_cells / total_cells) * 100,
    contribution_to_total_expanded = (expanded_cells / n_expanded_cells) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(contribution_to_total_expanded))

cat("\n=== EXPANDED CELLS BY TISSUE × STAGE ===\n")
print(expanded_by_tissue_stage)

# Add column showing contribution to the overall 31.5% expanded population
expanded_by_tissue_stage <- expanded_by_tissue_stage %>%
  mutate(
    contribution_to_31.5_percent = (expanded_cells / total_tcells) * 100
  )

cat("\n=== EXPANDED CELLS BY TISSUE × STAGE (with contribution to 31.5%) ===\n")
print(expanded_by_tissue_stage)

# Verification: this should sum to 31.5%
cat("\nVerification: Sum of contributions should equal 31.5%\n")
cat("Total:", sum(expanded_by_tissue_stage$contribution_to_31.5_percent), "%\n")

# Verification
cat("\nVerification: Sum of all contributions should equal 100%\n")
cat("By Tissue:", sum(expanded_by_tissue$contribution_to_total_expanded), "%\n")
cat("By Stage:", sum(expanded_by_stage$contribution_to_total_expanded), "%\n")

# Extract the clonotype sequences (CTaa) of expanded clones
expanded_clone_list <- expanded_clones %>%
  select(CTaa, n_cells) %>%
  arrange(desc(n_cells))  # optional: sort by clone size

# Expanded clonotype details

# Combine samples
combined_df <- do.call(rbind, combined.TCR)

# Keep only needed columns
combined_df <- combined_df[, c("CTaa", "Tissue", "Type", "sample")]

# Add cell count per clonotype
cell_counts <- combined_df %>%
  filter(!is.na(CTaa)) %>%
  group_by(CTaa) %>%
  summarise(n_cells = n(), .groups = "drop")

# Filter for expanded clonotypes
expanded_df <- combined_df %>%
  filter(CTaa %in% expanded_clone_list$CTaa) %>%
  filter(!is.na(CTaa)) %>%
  distinct()

# Collapse tissues for CTaa × Type × Sample
expanded_summary <- expanded_df %>%
  group_by(CTaa, Type, sample) %>%
  summarise(
    Tissues = paste(sort(unique(Tissue)), collapse = ","),
    .groups = "drop"
  )

# Final summary (with sample, type, tissues)
expanded_summary_single <- expanded_summary %>%
  group_by(CTaa) %>%
  summarise(
    Type = paste(sort(unique(Type)), collapse = ","),
    sample = paste(sort(unique(sample)), collapse = ","),
    Tissues = paste(sort(unique(Tissues)), collapse = ","),
    .groups = "drop"
  )

# Add number of cells per clonotype
expanded_summary_single <- expanded_summary_single %>%
  left_join(cell_counts, by = "CTaa") %>%
  arrange(desc(n_cells))

# Output
nrow(expanded_summary_single)
head(expanded_summary_single)

# Add single and multi tissue expansion info
# ️Standardize tissue names and classify shared_class
expanded_summary_single <- expanded_summary_single %>%
  mutate(
    Tissues = toupper(Tissues),
    tissue_list = str_split(Tissues, ","),
    shared_class = case_when(
      sapply(tissue_list, function(x) all(x == "PBMC")) ~ "a_PBMC",
      sapply(tissue_list, function(x) {
        solid <- intersect(x, c("LIVER", "SAT", "VAT"))
        length(solid) == 1 && !("PBMC" %in% x)
      }) ~ "b_One_solid_tissue",
      sapply(tissue_list, function(x) {
        solid <- intersect(x, c("LIVER", "SAT", "VAT"))
        length(solid) >= 2 && !("PBMC" %in% x)
      }) ~ "c_Shared_between_solid_tissues",
      TRUE ~ NA_character_
    )
  )

expanded_summary_single2 <- expanded_summary_single %>%
  left_join(
    expanded_clone_list %>% select(CTaa, n_cells_clone = n_cells), 
    by = "CTaa"
  ) %>%
  mutate(
    n_cells = ifelse(is.na(n_cells_clone), 0, n_cells_clone)
  ) %>%
  select(-n_cells_clone)  # optional, to remove the temporary column


# Standardize tissue list and stages, handle NA class
expanded_summary_single2 <- expanded_summary_single2 %>%
  rowwise() %>%
  mutate(
    tissue_list = list(str_split(Tissues, ",")[[1]]),
    Stage_list = list(str_split(Type, ",")[[1]]),
    solid_count = sum(unlist(tissue_list) %in% c("LIVER", "SAT", "VAT")),
    pbmc_present = "PBMC" %in% unlist(tissue_list),
    shared_class = coalesce(shared_class, "d_Other")
  ) %>%
  ungroup()

# Summarize per shared_class
expanded_summary_final <- expanded_summary_single2 %>%
  group_by(shared_class) %>%
  summarise(
    n_clonotypes = n_distinct(CTaa),                     # unique clonotypes
    n_cells_total = sum(n_cells),
    n_LIVER = sum(sapply(tissue_list, function(x) "LIVER" %in% x)),
    n_cells_LIVER = sum(sapply(1:n(), function(i) if("LIVER" %in% tissue_list[[i]]) n_cells[i] else 0)),
    n_SAT = sum(sapply(tissue_list, function(x) "SAT" %in% x)),
    n_cells_SAT = sum(sapply(1:n(), function(i) if("SAT" %in% tissue_list[[i]]) n_cells[i] else 0)),
    n_VAT = sum(sapply(tissue_list, function(x) "VAT" %in% x)),
    n_cells_VAT = sum(sapply(1:n(), function(i) if("VAT" %in% tissue_list[[i]]) n_cells[i] else 0)),
    n_PBMC = sum(sapply(tissue_list, function(x) "PBMC" %in% x)),
    n_cells_PBMC = sum(sapply(1:n(), function(i) if("PBMC" %in% tissue_list[[i]]) n_cells[i] else 0)),
    stages = paste(unique(unlist(Stage_list)), collapse = ", "),
    .groups = "drop"
  )

# Check result
expanded_summary_final

expanded_lookup <- expanded_summary_single2 %>%
  select(CTaa, shared_class) %>%
  distinct()

#### Load seurat object ####
seurat_obj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/TcellObj.rds')

# Initialize column
seurat_obj$expanded_shared_class <- NA_character_

# Map only for expanded clones
expanded_idx <- which(seurat_obj$Clonality == "Expanded")
seurat_obj$expanded_shared_class[expanded_idx] <- expanded_lookup$shared_class[
  match(seurat_obj$CTaa[expanded_idx], expanded_lookup$CTaa)
]

# Check
table(seurat_obj$expanded_shared_class, useNA = "ifany")

# Add barcodes from Seurat
seurat_meta <- seurat_obj@meta.data %>%
  mutate(barcode = rownames(seurat_obj@meta.data))

# Filter for expanded clones — Tissue already exists in metadata
expanded_cells <- seurat_meta %>%
  filter(Clonality == "Expanded") %>%
  select(CTaa, expanded_shared_class, barcode, Tissue) %>%  # Tissue is already here
  mutate(Tissue = toupper(Tissue))

# Sanity check
cat("Tissue NAs:", sum(is.na(expanded_cells$Tissue)), "out of", nrow(expanded_cells), "\n")
head(expanded_cells)

expanded_summary_tissue <- expanded_cells %>%
  group_by(expanded_shared_class, Tissue) %>%
  summarise(
    n_unique_clonotypes = n_distinct(CTaa),
    n_cells = n(),                  # total cells in this tissue and class
    .groups = "drop"
  ) %>%
  arrange(expanded_shared_class, Tissue)

expanded_summary_tissue

Idents(seurat_obj) <- "expanded_shared_class"

# Expanded shared class 
# Subset Seurat object to keep only selected shared_class categories
selected_classes <- c(
  "a_PBMC",
  "b_One_solid_tissue",
  "c_Shared_between_solid_tissues"
)

seurat_obj_subset <- subset(
  seurat_obj,
  subset = expanded_shared_class %in% selected_classes
)

# Quick check
table(seurat_obj_subset$expanded_shared_class, useNA = "ifany")

# Extract metadata
meta <- seurat_obj_subset@meta.data %>%
  mutate(
    Tissue = toupper(Tissue)  # standardize tissue names
  )

# Summarize
expanded_summary_table <- meta %>%
  group_by(expanded_shared_class, Tissue) %>%
  summarise(
    n_unique_clonotypes = n_distinct(CTaa),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(expanded_shared_class, Tissue)

expanded_summary_table

# Standardize tissue names
meta <- seurat_obj_subset@meta.data %>%
  mutate(Tissue = toupper(Tissue))

# Filter for shared_between_solid_tissues
shared_cells <- meta %>%
  filter(expanded_shared_class == "c_Shared_between_solid_tissues") %>%
  select(CTaa, Tissue)

clonotype_tissues <- shared_cells %>%
  group_by(CTaa) %>%
  summarise(
    tissues = list(sort(unique(Tissue))),
    .groups = "drop"
  ) %>%
  mutate(
    tissue_combination = sapply(tissues, function(x) paste(x, collapse = " & "))
  )

clonotype_tissue_summary <- clonotype_tissues %>%
  group_by(tissue_combination) %>%
  summarise(
    n_clonotypes = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_clonotypes))

clonotype_tissue_summary

#### Unique clonotypes and T cells sharing details ####
meta %>%
  filter(expanded_shared_class %in% selected_classes) %>%
  group_by(expanded_shared_class) %>%
  summarise(
    n_unique_clonotypes = n_distinct(CTaa),
    n_cells = n(),
    .groups = "drop"
  )

# Save expadned obj
saveRDS(seurat_obj_subset, '/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinal.rds')

#### Load expanded obj ####
seurat_obj_subset <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinal.rds')

# Extract metadata
meta <- seurat_obj_subset@meta.data

# Summarize number of unique clonotypes and total cells
expanded_summary_table_calc <- meta %>%
  group_by(Stage, Tissue) %>%
  summarise(
    n_unique_clonotypes = n_distinct(CTaa),
    n_cells = n()
  ) %>%
  arrange(Stage, Tissue)

expanded_summary_table_calc


# Summarize number of unique clonotypes and total cells
expanded_summary_table_calc <- meta %>%
  group_by(expanded_shared_class, Tissue) %>%
  summarise(
    n_unique_clonotypes = n_distinct(CTaa),
    n_cells = n()
  ) %>%
  arrange(expanded_shared_class, Tissue)

expanded_summary_table_calc

exact_summary <- meta %>%
  filter(expanded_shared_class %in% c("b_One_solid_tissue", "c_Shared_between_solid_tissues")) %>%
  group_by(expanded_shared_class) %>%
  summarise(
    total_unique_clonotypes = n_distinct(CTaa),
    total_cells = n(),
    .groups = "drop"
  )

exact_summary

discriminative_genes <- c(
  "CCR7", "SELL", "S1PR1",      # PBMC-enriched
  "ITGAL", "ITGB2", "ITGA4", "ITGB7", "CD44", "CXCR3",  # One solid tissue
  "CXCR6", "CCR5", "ITGAE",  "CX3CR1", # Shared solid tissues
  "ICAM1", "VCAM1",  "ITGB1" , "ADGRE5"          # additional discriminators
)

seurat_obj_subset$expanded_shared_class <- dplyr::recode(
  seurat_obj_subset$expanded_shared_class,
  "a_PBMC" = "PBMC",
  "b_One_solid_tissue" = "One tissue",
  "c_Shared_between_solid_tissues" = "Shared tissue"
)

# Convert to factor with desired order
seurat_obj_subset$expanded_shared_class <- factor(
  seurat_obj_subset$expanded_shared_class,
  levels = c("Shared tissue", "One tissue", "PBMC")
)

# Check
table(seurat_obj_subset$expanded_shared_class)

#### fig3d ####
DefaultAssay(seurat_obj_subset) <- "RNA"

ggp <- DotPlot2(
  seurat_obj_subset,
  features = discriminative_genes,
  group.by = "expanded_shared_class"
) + 
  # Expression: white → blue
  scale_fill_gradient(low = "white", high = "blue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Migration and adhesion associated genes")+
  coord_flip()

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("MigrationAdhesionAssocitedGenes.svg", width = 8, height = 4)   # width/height in cm by default for svg()
print(ggp)
dev.off()

#### fig3e ####

## Migration-associated ADTs
migration_ADTs <- c(
  "Hu.CD11a",        # Integrin αL, adhesion / migration
  "Hu.CD11b",        # Integrin αM, adhesion / migration
  "Hu.CD11c",        # Integrin αX, adhesion / migration
  "Hu.CD18",         # Integrin β2, adhesion / migration
  "Hu.CD29",         # Integrin β1, cell adhesion/migration
  "Hu.CD49a",        # Integrin α1, tissue homing
  "Hu.CD49b",        # Integrin α2, migration
  "Hu.CD49d",        # Integrin α4, tissue trafficking
  "Hu.CX3CR1",       # Chemokine receptor, migration/homing
  "Hu.CD62L"         # L-selectin, lymph node homing
)

DefaultAssay(seurat_obj_subset) <- "ADTonly"

ggp <- DotPlot2(
  seurat_obj_subset,
  features = migration_ADTs,
  group.by = "expanded_shared_class",
  legend_order = c("fill", "size") 
) + 
  # Expression: white → blue
  scale_fill_gradient(low = "white", high = "blue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Migration and adhesion associated ADTs")+
  coord_flip()

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("MigrationAdhesionAssocitedADTs.svg", width = 8, height = 4)   # width/height in cm by default for svg()
print(ggp)
dev.off()


#### fig3f ####
effector_memory_activation_genes <- c(
  "GZMH", "GNLY", "FGFBP2", "GZMB", "PRF1",
  "CCL5", "NKG7", "CTSW", "CST7",
  "TBX21", "KLF2", "HCST", "HOPX", "EFHD2",
  "RASSF1", "PLEK"
)

seurat_obj_subset_tissue <- subset(
  seurat_obj_subset,
  Tissue %in% c("LIVER", "SAT", "VAT")
)

seurat_obj_subset_tissue$expanded_shared_class <- factor(
  seurat_obj_subset_tissue$expanded_shared_class,
  levels = c("Shared tissue", "One tissue")
)

# Check
table(seurat_obj_subset_tissue$expanded_shared_class)

Idents(seurat_obj_subset_tissue) <- "expanded_shared_class"

DefaultAssay(seurat_obj_subset_tissue) <- 'RNA'

ggp <- DotPlot2(
  seurat_obj_subset_tissue,
  features = effector_memory_activation_genes,
  group.by = "expanded_shared_class",
  legend_order = c("fill", "size")
) + 
  scale_fill_gradient(low = "white", high = "blue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Effector Memory Activation Genes") +
  coord_flip()

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("OneVsSharedTcellGenes.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()

#### fig3g ####
innate_il_genes <- c(
  "TLR2", "TLR4", "NLRP3", "IL1B", "IL6", "IL10", 
  "IL12A", "IL12B", "IL18", "IFNB1", "CXCL8", "CCL2",
  "CCL3", "STAT1", "IRF7"
)

ggp <- DotPlot2(
  seurat_obj_subset_tissue,
  features = innate_il_genes,
  group.by = "expanded_shared_class",
  legend_order = c("fill", "size")
) + 
  scale_fill_gradient(low = "white", high = "blue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Innate and Interleukin Genes") +
  coord_flip()

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("OneVsSharedTcellGenes1.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()

#### fig3h,i ####

# fibrosis specific genes  in shared tissue 

seurat_obj_subset_tissue_shared <- subset(
  seurat_obj_subset_tissue,
  expanded_shared_class %in% "Shared tissue"
)

table(seurat_obj_subset_tissue_shared$expanded_shared_class)

# Remove Healthy
# Remove all scale.data layers except the main one
layers_to_remove <- grep("scale\\.data\\.", 
                         names(seurat_obj_subset_tissue_shared@assays$RNA@layers), 
                         value = TRUE)

for (l in layers_to_remove) {
  seurat_obj_subset_tissue_shared@assays$RNA@layers[[l]] <- NULL
}

# Check layers are clean now
names(seurat_obj_subset_tissue_shared@assays$RNA@layers)

# Now subset
cells_keep <- rownames(seurat_obj_subset_tissue_shared@meta.data)[
  seurat_obj_subset_tissue_shared@meta.data$Patient_ID != "10742"
]

seurat_obj_subset_tissue_shared <- seurat_obj_subset_tissue_shared[, cells_keep]

# Check
table(seurat_obj_subset_tissue_shared$Patient_ID)
table(seurat_obj_subset_tissue_shared$Stage)

# Set factor order
seurat_obj_subset_tissue_shared$Stage <- factor(
  seurat_obj_subset_tissue_shared$Stage,
  levels = c("Fibrosis", "No_fibrosis")
)

# Check cell counts per group
table(seurat_obj_subset_tissue_shared$Stage)
# Fibrosis No_fibrosis 
#    1769        135

DefaultAssay(seurat_obj_subset_tissue_shared) <- "RNA"

markers <- FindMarkers(
  seurat_obj_subset_tissue_shared,
  ident.1 = "Fibrosis",
  ident.2 = "No_fibrosis",
  group.by = "Stage",
  test.use = "wilcox"  # default
)

# Filter for significant genes (adjusted p-value < 0.05)
sig_genes <- markers %>% filter(p_val_adj < 0.05)
print(sig_genes)

fibrosis_effector_genes <- c(
  "TRAV38-1",
  "CCL3L1",
  "CCL4L2",
  "XCL2",
  "CCL13",
  "JUN",
  "GAS5"
)

# Reorder Tissue factor in the Seurat object
seurat_obj_subset_tissue_shared$Tissue <- 
  factor(seurat_obj_subset_tissue_shared$Tissue, levels = c("LIVER", "SAT", "VAT"))

DefaultAssay(seurat_obj_subset_tissue_shared) <- "RNA"

# 4. Run DotPlot2
ggp <- DotPlot2(
  seurat_obj_subset_tissue_shared,
  features = fibrosis_effector_genes,
  group.by = "Tissue",
  split.by = "Stage",
  legend_order = c("fill", "color", "size")   # Correct legend order
) +
  scale_fill_gradient(low = "white", high = "blue") +
  scale_colour_manual(
    values = c(
      "Fibrosis" = "#084594",      # dark blue
      "No fibrosis" = "#9ecae1"    # light blue
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Fibrosis - No Fibrosis shared T cells")

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("fibrosisSharedGenes.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()

# fibrosis specific ADT in shared tissue 

DefaultAssay(seurat_obj_subset_tissue_shared) <- "ADTonly"

# Subset Seurat object to only include cells with ADT counts
adt_cells <- colnames(seurat_obj_subset_tissue_shared[["ADTonly"]])
seurat_obj_ADT <- subset(seurat_obj_subset_tissue_shared, cells = adt_cells)

# Set the Tissue factor order
seurat_obj_ADT$Tissue <- factor(seurat_obj_ADT$Tissue, levels = c("LIVER", "SAT", "VAT"))

# Make sure ADT assay is active
DefaultAssay(seurat_obj_ADT) <- "ADTonly"

markers <- FindMarkers(
  seurat_obj_ADT,
  ident.1 = "Fibrosis",
  ident.2 = "No_fibrosis",
  group.by = "Stage",
  test.use = "wilcox"  # default
)

# Filter for significant genes (adjusted p-value < 0.05)
sig_adt <- markers %>% filter(p_val_adj < 0.05)
print(sig_adt)

fibrosis_tcell_adts <- c(
  "Hu.CD7",        # Pan–T cell marker; relative downshift often reflects altered T-cell composition or activation in fibrosis
  "Hu.CD45RA",     # Typically naive, BUT in fibrosis often marks TEMRA CD8+ T cells (CD45RA+ CCR7−), linked to chronic tissue damage
  "Hu.TCR.Vd2",    # γδ T cells; expanded in chronic inflammation and tissue fibrosis, promoting fibroblast activation
  "Hu.CD25",       # IL-2Rα; fibrosis-associated activation and/or Treg expansion in chronic inflammatory niches
  "Hu.CD127",      # IL-7Rα; reduced relative to CD25 suggests effector/Treg skewing in fibrosis
  "Hu.CD154",      # CD40L; activated CD4 T cells driving macrophage–fibroblast crosstalk
  "Hu.CD40",       # Costimulatory receptor; reflects heightened immune activation in fibrotic tissue
  "Hu.CD152",      # CTLA-4; immune checkpoint upregulation due to chronic antigen stimulation in fibrosis
  "Hu.CD272",      # BTLA; inhibitory receptor marking exhausted/regulated T cells in chronic fibrotic inflammation
  "Hu.CD107a",     # Degranulation marker; increased cytotoxic T-cell activity contributing to epithelial/tissue injury
  "Hu.KLRG1",      # Terminal effector / senescent T cells; hallmark of chronic antigen exposure in fibrosis
  "Hu.CD49d",      # Integrin α4; enhanced tissue homing and retention of T cells in fibrotic organs
  "Hu.CD49b",      # Integrin α2; effector/memory T cells interacting with collagen-rich fibrotic matrix
  "Hu.CD71"        # Proliferation/metabolic activation; increased turnover of activated T cells in fibrosis
)

# 4. Run DotPlot2
ggp <- DotPlot2(
  seurat_obj_ADT,
  features = fibrosis_tcell_adts,
  group.by = "Tissue",
  split.by = "Stage",
  legend_order = c("fill", "color", "size")   # Correct legend order
) +
  scale_fill_gradient(low = "white", high = "blue") +
  scale_colour_manual(
    values = c(
      "Fibrosis" = "#084594",      # dark blue
      "No_fibrosis" = "#9ecae1"    # light blue
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Fibrosis - No Fibrosis shared T cells")

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("fibrosisSharedADTs.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()

# End of the script
