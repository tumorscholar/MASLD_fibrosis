library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(scRepertoire)
library(dplyr)
library(SeuratExtend)
library(ggupset)  
library(forcats)
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

# Write to CSV
write_csv(expanded_by_tissue_stage, "/data/home/hdx044/files/screpertoire/expanded_by_tissue_stage.csv")


# Extract the clonotype sequences (CTaa) of expanded clones
expanded_clone_list <- expanded_clones %>%
 select(CTaa, n_cells) %>%
 arrange(desc(n_cells))  # optional: sort by clone size

# Write to CSV
write_csv(expanded_clone_list, "/data/home/hdx044/files/screpertoire/expanded_clonotypes_3098.csv")

cat("✅ File 'expanded_clonotypes_3098.csv' written successfully.\n")

#### Expanded clonotype details ####

library(dplyr)
library(tidyr)

# 1️⃣ Combine samples
combined_df <- do.call(rbind, combined.TCR)

# 2️⃣ Keep only needed columns
combined_df <- combined_df[, c("CTaa", "Tissue", "Type", "sample")]

# 3️⃣ Add cell count per clonotype
cell_counts <- combined_df %>%
 filter(!is.na(CTaa)) %>%
 group_by(CTaa) %>%
 summarise(n_cells = n(), .groups = "drop")

# 4️⃣ Filter for expanded clonotypes
expanded_df <- combined_df %>%
 filter(CTaa %in% expanded_clone_list$CTaa) %>%
 filter(!is.na(CTaa)) %>%
 distinct()

# 5️⃣ Collapse tissues for CTaa × Type × Sample
expanded_summary <- expanded_df %>%
 group_by(CTaa, Type, sample) %>%
 summarise(
  Tissues = paste(sort(unique(Tissue)), collapse = ","),
  .groups = "drop"
 )

# 6️⃣ Final summary (with sample, type, tissues)
expanded_summary_single <- expanded_summary %>%
 group_by(CTaa) %>%
 summarise(
  Type = paste(sort(unique(Type)), collapse = ","),
  sample = paste(sort(unique(sample)), collapse = ","),
  Tissues = paste(sort(unique(Tissues)), collapse = ","),
  .groups = "drop"
 )

# 7️⃣ Add number of cells per clonotype
expanded_summary_single <- expanded_summary_single %>%
 left_join(cell_counts, by = "CTaa") %>%
 arrange(desc(n_cells))

# Output
nrow(expanded_summary_single)
head(expanded_summary_single)

write_csv(expanded_summary_single, "/data/home/hdx044/files/screpertoire/expanded_clonotypes_3098_details.csv")

#### Add single and multi tissue expansion info ####
library(dplyr)
library(stringr)

# 1️⃣ Standardize tissue names and classify shared_class
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


# 3️⃣ Standardize tissue list and stages, handle NA class
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

# 4️⃣ Summarize per shared_class
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

# 5️⃣ Check result
expanded_summary_final


write_csv(expanded_summary_final, "/data/home/hdx044/files/screpertoire/expanded_clonotypes_3098_details_TissueWise.csv")

#### Load seurat object ####
SeuObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

#### Check tolerance vs activated ####

tolerance_genes <- c("PDCD1LG2", "PDCD1", "CTLA4", "TIGIT", "LAG3", 
                     "FOXP3", "IL10", "IL10RA", "IDO1", "HAVCR2", "BTLA")

damage_genes <- c("GZMB", "GNLY", "PRF1", "NKG7", "CCL5",
                  "IFNG", "TBX21", "CXCL9", "CXCL10", "CD38", "ENTPD1")


# Add module scores with explicit naming
seurat_obj <- AddModuleScore(
 seurat_obj,
 features = list(tolerance_genes, damage_genes),
 name = c("Tolerance_gene_score", "Damage_gene_score")
)

# Check the metadata column names that were created
colnames(seurat_obj@meta.data)[grepl("Tolerance_gene_score|Damage_gene_score", colnames(seurat_obj@meta.data))]

library(dplyr)
library(stringr)

seurat_obj$ExpandedClones <- recode(
 seurat_obj$ExpandedClones,
 "Yes" = "Expanded",
 "No"  = "Non_expanded"
)

# Check
table(seurat_obj$ExpandedClones)

seurat_obj <- RenameCells(
 object = seurat_obj,
 new.names = gsub("_score1$", "_score", colnames(seurat_obj))
)

ggp <- VlnPlot2(seurat_obj, 
         features = c("Tolerance_gene_score", "Damage_gene_score"),  
         pt = FALSE, group.by = "Tissue", 
         split.by = "ExpandedClones",
         stat.method = "wilcox.test", 
         hide.ns = TRUE) +
 ggtitle("Tolerance and Damage Genes Scores") 

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("Tolerance_and_Damage_Genes_Scores.svg", width = 10, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()

####  Single vs Expanded genes #####

# Ensure the factor is set correctly
seurat_obj$ExpandedClones <- factor(seurat_obj$ExpandedClones, levels = c("Non_expanded", "Expanded"))

# Set the default assay (RNA or whichever contains gene expression)
DefaultAssay(seurat_obj) <- "RNA"

# Perform differential expression
de_markers <- FindMarkers(
 object = seurat_obj,
 ident.1 = "Expanded",
 ident.2 = "Non_expanded",
 group.by = "ExpandedClones",
 logfc.threshold = 0.25,   
 min.pct = 0.25    
)

# View top genes sorted by log fold change
de_markers <- de_markers[order(de_markers$avg_log2FC, decreasing = TRUE), ]
head(de_markers)

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

# Then run DotPlot2
ggp <- DotPlot2(
 seurat_obj,
 features = expanded_genes,
 group.by = "Tissue",
 split.by = "ExpandedClones",
 legend_order = c("fill", "color", "size")  # Your preferred order+
) +
 scale_fill_gradient(low = "white", high = "blue") +
 scale_colour_manual(
  values = c(
   "Expanded" = "#084594",     # dark blue
   "Non_expanded" = "#9ecae1"  # light blue
  )
 ) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 ggtitle("Expanded vs Non-expanded T cells")

ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("ExpandedNonexpandedTcellsDifferentialGenes.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()

# Set the default assay (RNA or whichever contains gene expression)
DefaultAssay(seurat_obj) <- "ADTonly"

# Perform differential expression
de_markers <- FindMarkers(
 object = seurat_obj,
 ident.1 = "Expanded",
 ident.2 = "Non_expanded",
 group.by = "ExpandedClones",
 logfc.threshold = 0.25,   
 min.pct = 0.25    
)

# View top genes sorted by log fold change
de_markers <- de_markers[order(de_markers$avg_log2FC, decreasing = TRUE), ]
print(de_markers)

expanded_ADTs <- c(
 "Hu.CD57",
 "Hu.CD8",
 "Hu.KLRG1",
 "Hu.CD5",
 "Hu.CD4-RPA.T4",
 "Hu.CD27",
 "Hu.CD62L",
 "Hu.CD7"
)

# 1. Set default assay to ADT
DefaultAssay(seurat_obj) <- "ADTonly"

# 2. Select all cells with ADT counts > 0
adt_cells <- colnames(seurat_obj)[Matrix::colSums(seurat_obj[["ADTonly"]]@counts) > 0]

# 3. Subset Seurat object
seurat_obj_ADT <- subset(seurat_obj, cells = adt_cells)

# 4. Check tissue composition
table(seurat_obj_ADT$Tissue)

# 5. Keep only features that exist in the ADT assay
valid_features <- expanded_ADTs[expanded_ADTs %in% rownames(seurat_obj_ADT[["ADTonly"]])]

# 6. Set Idents to ExpandedClones
Idents(seurat_obj_ADT) <- "ExpandedClones"

# 7. DotPlot using valid features only

ggp <- DotPlot2(
 seurat_obj_ADT,
 features = valid_features,
 group.by = "Tissue",
 split.by = "ExpandedClones"
) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 ggtitle("Expanded vs Non-expanded T cells (ADTs)") +
 coord_flip()


# 8. Print the plot
ggp

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("ExpandedNonexpandedTcellsDifferentialGenes.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()


library(dplyr)
library(purrr)
library(stringr)

# Merge all individual sample TCR data frames into one
combined_TCR_df <- combined.TCR %>%
 map_df(~as.data.frame(.x), .id = "sample_id")

# Quick check
colnames(combined_TCR_df)
# Should include: barcode, CTaa, Type, Tissue

# Count number of cells per clonotype
clone_counts <- combined_TCR_df %>%
 group_by(CTaa) %>%
 summarise(n_cells = n(), .groups = "drop")

# Expanded clones (>1 cell)
expanded_clones <- clone_counts$CTaa[clone_counts$n_cells > 1]

# Barcodes of expanded clones
expanded_barcodes <- combined_TCR_df$barcode[combined_TCR_df$CTaa %in% expanded_clones]

# Add Clonality metadata
SeuObj$Clonality <- ifelse(
 colnames(SeuObj) %in% expanded_barcodes,
 "Expanded",
 "Non-expanded"
)

# Quick check
table(SeuObj$Clonality)

# Example: use expanded_summary_single2 or expanded_summary_single
expanded_lookup <- expanded_summary_single2 %>%
 select(CTaa, shared_class) %>%
 distinct()

# Initialize column
seurat_obj$expanded_shared_class <- NA_character_

# Map only for expanded clones
expanded_idx <- which(seurat_obj$Clonality == "Expanded")
seurat_obj$expanded_shared_class[expanded_idx] <- expanded_lookup$shared_class[
 match(seurat_obj$CTaa[expanded_idx], expanded_lookup$CTaa)
]

# Check
table(seurat_obj$expanded_shared_class, useNA = "ifany")

library(dplyr)
library(tidyr)

# Add barcodes from Seurat
seurat_meta <- seurat_obj@meta.data %>%
 mutate(barcode = rownames(seurat_obj@meta.data))

# Filter for expanded clones
expanded_cells <- seurat_meta %>%
 filter(Clonality == "Expanded") %>%
 select(CTaa, expanded_shared_class, barcode)

expanded_cells <- expanded_cells %>%
 left_join(
  combined_TCR_df %>% select(barcode, Tissue),
  by = "barcode"
 ) %>%
 mutate(Tissue = toupper(Tissue))

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

#### Expanded shared class ####

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

library(dplyr)

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

library(dplyr)
library(tidyr)

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

# Save expadned obj
saveRDS(seurat_obj_subset, '/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinal.rds')


# Expanded T cells
ExpandedTcells <- subset(SeuObj, subset = clone_category %in% c("expanded","hyperexpanded", "singleton"))

ExpandedTcells$clone_group <- ifelse(
 ExpandedTcells$clone_category %in% c("expanded", "hyperexpanded"),
 "expanded+hyperexpanded",
 "singleton"
)

#### check autoaggressive cell count ####
DefaultAssay(ExpandedTcells) <- "RNA"
autoagressive_cells <- WhichCells(ExpandedTcells, expression = CD8A > 0.25 & CD69 > 0.25 & CXCR6 > 0.25)

AutoaggressiveTcells <- subset(
 ExpandedTcells,
 cells = autoagressive_cells
)

AutoaggressiveTcells
table(AutoaggressiveTcells$clone_category)

AutoaggressiveTcells$clone_group <- ifelse(
 AutoaggressiveTcells$clone_category %in% c("expanded", "hyperexpanded"),
 "expanded+hyperexpanded",
 "singleton"
)

table(AutoaggressiveTcells$clone_group)

library(dplyr)

auto_counts <- AutoaggressiveTcells@meta.data %>%
 count(Stage, Patient_ID, Tissue, clone_group, name = "n_auto")

total_counts <- ExpandedTcells@meta.data %>%
 count(Stage, Patient_ID, Tissue, clone_group, name = "n_total")

auto_percent <- auto_counts %>%
 left_join(total_counts,
           by = c("Stage", "Patient_ID", "Tissue", "clone_group")) %>%
 mutate(
  percent_autoaggressive = 100 * n_auto / n_total
 ) %>%
 arrange(Stage, Patient_ID, Tissue, clone_group)

auto_percent


library(tidyr)

auto_percent_wide <- auto_percent %>%
 select(Stage, Patient_ID, Tissue, clone_group, percent_autoaggressive) %>%
 pivot_wider(
  names_from = clone_group,
  values_from = percent_autoaggressive
 )

auto_percent_wide

write_csv(auto_percent_wide, "/data/home/hdx044/files/screpertoire/AutoagressiveTcells_details_TissueStageWise.csv")

#### check protective cell count ####

protective_cells <- WhichCells(ExpandedTcells, expression = CD8A > 0.25 & CD69 > 0.25 & ITGAE < 0.25)

protectiveTcells <- subset(
 ExpandedTcells,
 cells = protective_cells
)

table(protectiveTcells$clone_category)

protectiveTcells$clone_group <- ifelse(
 protectiveTcells$clone_category %in% c("expanded", "hyperexpanded"),
 "expanded+hyperexpanded",
 "singleton"
)

table(protectiveTcells$clone_group)

library(dplyr)

auto_counts <- protectiveTcells@meta.data %>%
 count(Stage, Patient_ID, Tissue, clone_group, name = "n_auto")

total_counts <- ExpandedTcells@meta.data %>%
 count(Stage, Patient_ID, Tissue, clone_group, name = "n_total")

auto_percent <- auto_counts %>%
 left_join(total_counts,
           by = c("Stage", "Patient_ID", "Tissue", "clone_group")) %>%
 mutate(
  protectiveTcells = 100 * n_auto / n_total
 ) %>%
 arrange(Stage, Tissue, clone_group)

auto_percent


library(tidyr)

auto_percent_wide <- auto_percent %>%
 select(Stage, Patient_ID, Tissue, clone_group, protectiveTcells) %>%
 pivot_wider(
  names_from = clone_group,
  values_from = protectiveTcells
 )

auto_percent_wide
write_csv(auto_percent_wide, "/data/home/hdx044/files/screpertoire/ProtectiveTcells_details_TissueStageWise.csv")

#### Load expadned obj ####
seurat_obj_subset <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinal.rds')

# Assuming your Seurat object is seurat_obj_subset
# and metadata columns are named 'Tissue' and 'expanded_shared_class'

# Extract metadata
meta <- seurat_obj_subset@meta.data

# Summarize number of unique clonotypes and total cells
expanded_summary_table_calc <- meta %>%
 group_by(Stage, Tissue) %>%
 summarise(
  n_unique_clonotypes = n_distinct(CTaa), # replace with your clonotype column name
  n_cells = n()
 ) %>%
 arrange(Stage, Tissue)

expanded_summary_table_calc


# Summarize number of unique clonotypes and total cells
expanded_summary_table_calc <- meta %>%
 group_by(expanded_shared_class, Tissue) %>%
 summarise(
  n_unique_clonotypes = n_distinct(CTaa), # replace with your clonotype column name
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

#### fig.3D dot plot ####
DefaultAssay(seurat_obj_subset_tissue) <- "RNA"

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


#### fig.4a Plot effector-memory genes ####
effector_memory_activation_genes <- c(
 "GZMH", "GNLY", "FGFBP2", "GZMB", "PRF1",
 "CCL5", "NKG7", "CTSW", "CST7",
 "TBX21", "KLF2", "HCST", "HOPX", "EFHD2",
 "RASSF1", "PLEK"
)

innate_il_genes <- c(
 "TLR2", "TLR4", "NLRP3", "IL1B", "IL6", "IL10", 
 "IL12A", "IL12B", "IL18", "IFNB1", "CXCL8", "CCL2",
 "CCL3", "STAT1", "IRF7"
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

DefaultAssay(seurat_obj_subset) <- "RNA"
Idents(seurat_obj_subset) <- seurat_obj_subset$expanded_shared_class

subset_obj <- subset(seurat_obj_subset, expanded_shared_class %in% c("Shared tissue", "One tissue"))
Idents(subset_obj) <- subset_obj$expanded_shared_class

adt_markers <- FindMarkers(
 object = subset_obj,
 ident.1 = "Shared tissue",
 ident.2 = "One tissue",
 assay = "ADTonly",
 logfc.threshold = 0.1,        # set threshold as needed
 min.pct = 0.25              # expressed in at least 10% of cells
)

# Order by log fold-change and filter for significance
adt_markers_sig <- adt_markers %>%
 rownames_to_column(var = "ADT_tag") %>%   # if rownames are protein names
 filter(p_val_adj < 0.05) %>%
 arrange(desc(avg_log2FC))

# Print results
print(adt_markers_sig)

effector_memory_ADTs <- c(
 "Hu.KLRG1",      # senescence / effector memory
 "Hu.CD16",       # cytotoxic effector
 "Hu.TCR.Vd2",    # γδ T cells, often effector
 "Hu.CD268",      # activating receptor (NK/T cells)
 "Hu.GPR56",      # effector/memory marker
 "Hu.CD49b",      # integrin, cytotoxic cells
 "Hu.CD163"      # tissue resident / macrophage effector (if relevant)
)

# Set default assay to ADT
DefaultAssay(seurat_obj_subset_tissue) <- "ADTonly"

# 1. Subset Seurat object to only include cells with ADT counts
adt_cells_tissue <- colnames(seurat_obj_subset_tissue[["ADTonly"]])
seurat_obj_ADT_tissue <- subset(seurat_obj_subset_tissue, cells = adt_cells_tissue)

# 2. Set the Tissue factor order
seurat_obj_ADT_tissue$Tissue <- factor(seurat_obj_ADT_tissue$Tissue, levels = c("LIVER", "SAT", "VAT"))

Idents(seurat_obj_ADT_tissue) <- "expanded_shared_class"

# 4. Run DotPlot2
DefaultAssay(seurat_obj_ADT_tissue) <- "ADTonly"

# Only keep features that exist
valid_features <- effector_memory_ADTs[effector_memory_ADTs %in% rownames(seurat_obj_ADT_tissue[["ADTonly"]])]

ggp <- DotPlot2(
 seurat_obj_ADT_tissue,
 features = valid_features,
 legend_order = c("fill", "size")
) +
 scale_fill_gradient(low = "white", high = "blue") +
 theme(axis.text.x = element_text(angle = 45, hjust = 1),
       legend.position = "right") +
 ggtitle("Effector Memory activation markers")+  coord_flip() 

print(ggp)

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("sharedADTs.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()

#### fig.4b fibrosis specific gens in shared tissue ####

seurat_obj_subset_tissue_shared <- subset(
 seurat_obj_subset_tissue,
 expanded_shared_class %in% "Shared tissue"
)

table(seurat_obj_subset_tissue_shared$expanded_shared_class)

# 1. Remove Healthy
seurat_obj_subset_tissue_shared <- subset(
 seurat_obj_subset_tissue_shared,
 Stage != "c_Healthy"
)

# 2. Recode the remaining groups
seurat_obj_subset_tissue_shared$Stage <- dplyr::recode(
 seurat_obj_subset_tissue_shared$Stage,
 "a_Fibrosis" = "Fibrosis",
 "b_No_Fibrosis" = "No fibrosis"
)

# 3. Set factor order (Fibrosis first, as typical in MASH/MASLD analysis)
seurat_obj_subset_tissue_shared$Stage <- factor(
 seurat_obj_subset_tissue_shared$Stage,
 levels = c("Fibrosis", "No fibrosis")
)

# Check cell counts per group
table(seurat_obj_subset_tissue_shared$Stage)
# Fibrosis No Fibrosis 
#    1798        136


DefaultAssay(seurat_obj_subset_tissue_shared) <- "ADTonly"
DefaultAssay(seurat_obj_subset_tissue_shared) <- "RNA"

markers <- FindMarkers(
 seurat_obj_subset_tissue_shared,
 ident.1 = "Fibrosis",
 ident.2 = "No Fibrosis",
 group.by = "Stage",
 test.use = "wilcox"  # default
)

# Filter for significant genes (adjusted p-value < 0.05)
sig_genes <- markers %>% filter(p_val_adj < 0.05)
print(sig_genes)

# Optional: filter for upregulated in Fibrosis
up_genes <- sig_genes %>% filter(avg_log2FC > 0)

# View top genes
print(up_genes)

fibrosis_effector_genes <- c(
 "TRAV38-1",
 "CCL3L1",
 "CCL4L2",
 "XCL2",
 "CCL13",
 "JUN",
 "GAS5"
)

fibrosis_effector_ADTs <- c(
 "Hu.CD267",    # ICOS, activation marker
 "Hu.CD183",    # CXCR3, Th1/memory migration
 "Hu.CD154",    # CD40L, T cell activation
 "Hu.CD23",     # FCER2, B/T cell interactions
 "Hu.CD58",     # LFA-3, T cell adhesion & activation
 "Hu.CD40",     # CD40, co-stimulation
 "Hu.CD124",    # IL-4Rα, cytokine receptor
 "Hu.CD112",    # Nectin-2, adhesion & immune activation
 "Hu.HLA.E",    # MHC class I, antigen presentation
 "Hu.CD185",    # CXCR5, migration to follicles / memory
 "Hu.CD223",    # LAG3, inhibitory checkpoint (often activated T cells)
 "Hu.CD319",    # SLAMF7, cytotoxic / NK interactions
 "Hu.CD105-43A3", # Endoglin, adhesion & signaling
 "Hu.CD64",     # FcγRI, phagocytosis / activation
 "Hu.CD83",     # Maturation / activation marker
 "Hu.CD95",     # Fas, activation-induced apoptosis
 "Hu.CD24"      # Adhesion / memory T cell maturation
)

library(forcats)

# Reorder Tissue factor in the Seurat object
seurat_obj_subset_tissue_shared$Tissue <- 
 factor(seurat_obj_subset_tissue_shared$Tissue, levels = c("LIVER", "SAT", "VAT"))

DefaultAssay(seurat_obj_subset_tissue_shared) <- "ADTonly"

# 1. Subset Seurat object to only include cells with ADT counts
adt_cells <- colnames(seurat_obj_subset_tissue_shared[["ADTonly"]])
seurat_obj_ADT <- subset(seurat_obj_subset_tissue_shared, cells = adt_cells)

# 2. Set the Tissue factor order
seurat_obj_ADT$Tissue <- factor(seurat_obj_ADT$Tissue, levels = c("LIVER", "SAT", "VAT"))

# 3. Make sure ADT assay is active
DefaultAssay(seurat_obj_ADT) <- "ADTonly"

# 4. Run DotPlot2
ggp <- DotPlot2(
 seurat_obj_ADT,
 features = fibrosis_effector_ADTs,
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
svg("fibrosisSharedADTs.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()

DefaultAssay(seurat_obj_subset_tissue_shared) <- "RNA"

# Then run DotPlot
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

setwd("/data/home/hdx044/plots/seurat/allTcells")
svg("fibrosisSharedGenes.svg", width = 5, height = 5)   # width/height in cm by default for svg()
print(ggp)
dev.off()


#### Characteristic of c_Shared_between_solid_tissues", "b_One_solid_tissue ####

Idents(seurat_obj_subset) <- seurat_obj_subset$expanded_shared_class

DefaultAssay(seurat_obj_subset) <- "RNA"

dge_results <- FindMarkers(
 seurat_obj_subset,
 ident.1 = "c_Shared_between_solid_tissues",
 ident.2 = "b_One_solid_tissue",
 logfc.threshold = 0.25,
 min.pct = 0.25,
 test.use = "wilcox"
)

write.csv(dge_results, file = "/data/home/hdx044/files/screpertoire/ExpandedClone/Differential_genes_SHARED_CLONES.csv")


dge_results <- FindMarkers(
 seurat_obj_subset,
 ident.1 = "c_Shared_between_solid_tissues",
 ident.2 = "b_One_solid_tissue")

write.csv(dge_results, file = "/data/home/hdx044/files/screpertoire/ExpandedClone/Differential_genes_SHARED_CLONES_relaxed.csv")

# Upregulated genes (higher in shared clones)
up_genes <- rownames(dge_results[dge_results$avg_log2FC > 0 & dge_results$p_val < 0.05, ])


#### Find genes of tissue shared clones in PBMC ####
# 1. load metadata
meta <- seurat_obj_subset@meta.data

# 2. quick check
table(meta$expanded_shared_class, useNA = "ifany")
# expect to see "c_Shared_between_solid_tissues" in the output

# Defensive exact-match for the class label
class_label <- "c_Shared_between_solid_tissues"

# Subset rows flagged with that label (use grepl if labels can have extra text)
rows_shared <- meta[grepl(class_label, meta$expanded_shared_class, fixed = TRUE), ]

# Build a data.frame with CTaa and Tissue (base-R)
shared_ctaa_df <- rows_shared[, c("CTaa", "Tissue"), drop = FALSE]

# Replace empty strings with NA for clarity
shared_ctaa_df$CTaa[shared_ctaa_df$CTaa == ""] <- NA

# Show top rows and counts
head(shared_ctaa_df, 20)
cat("Total rows in this class:", nrow(shared_ctaa_df), "\n")
cat("Unique CTaa in this class:", length(unique(shared_ctaa_df$CTaa[!is.na(shared_ctaa_df$CTaa)])), "\n")


#### Venn diagram expanded clonetypes (3098) ####

library(dplyr)
library(tidyr)
library(ggVennDiagram)

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

# Step 1: Expand tissue column so each clonotype has one row per tissue
expanded_clonotypes_long <- expanded_combined %>%
 separate_rows(tissues, sep = ",") %>%
 distinct(CTaa, tissues)

# Step 2: Check total unique expanded clonotypes
length(unique(expanded_clonotypes_long$CTaa))  # should return 3098

# Step 3: Optional sanity check: sum per tissue
expanded_tissue_summary <- expanded_clonotypes_long %>%
 group_by(tissues) %>%
 summarise(n_expanded_clonotypes = n_distinct(CTaa), .groups = "drop")

sum(expanded_tissue_summary$n_expanded_clonotypes)  # returns 4989

# Total unique expanded clonotypes across all tissues
length(unique(expanded_clonotypes_long$CTaa))  # should return 3098

# Sum per tissue (will be >3098 due to overlaps)
sum(expanded_tissue_summary$n_expanded_clonotypes)  # returns 4989


# Step 1: Prepare long-format clonotype table
expanded_long <- expanded_combined %>%
 rename(Stage = Type) %>%  # Rename to "Stage" for clarity
 separate_rows(Stage, sep = ",") %>%
 mutate(Stage = trimws(Stage)) %>%
 separate_rows(tissues, sep = ",") %>%
 mutate(tissues = trimws(tissues)) %>%
 distinct(CTaa, Stage, tissues) %>%
 filter(tissues != "Unknown")

# Step 2: Define Venn plotting function
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

# Step 3: Generate Venn diagrams per stage
venn_fib     <- plot_stage_venn("Fibrosis", color = "#e41a1c")
venn_no_fib  <- plot_stage_venn("No_fibrosis", color = "#377eb8")
venn_healthy <- plot_stage_venn("Healthy", color = "#4daf4a")

# Step 5: Save plots
out_dir <- "/data/home/hdx044/plots/screpertoire"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "Venn_Fibrosis.svg"), plot = venn_fib, width = 6, height = 6)
ggsave(file.path(out_dir, "Venn_No_Fibrosis.svg"), plot = venn_no_fib, width = 6, height = 6)
ggsave(file.path(out_dir, "Venn_Healthy.svg"), plot = venn_healthy, width = 6, height = 6)

cat("✅ Saved 3 SVG plots in:", out_dir, "\n")


#### Venn diagram  shared expanded clonetypes (1384) ####
library(dplyr)
library(tidyr)
library(ggVennDiagram)

# Prepare long-format table (tissue + stage info)
expanded_long <- expanded_combined %>%
 separate_rows(Type, sep = ",") %>%
 mutate(Type = trimws(Type)) %>%
 separate_rows(tissues, sep = ",") %>%
 mutate(tissues = trimws(tissues)) %>%
 distinct(CTaa, Type, tissues, n_cells) %>%
 filter(tissues != "Unknown")

clone_sharing <- expanded_long %>%
 group_by(CTaa) %>%
 summarise(
  n_tissues = n_distinct(tissues),
  tissues_present = paste(unique(tissues), collapse = ", "),
  .groups = "drop"
 ) %>%
 mutate(shared_status = ifelse(n_tissues > 1, "Shared_between_tissues", "Unique_to_one_tissue"))

expanded_long <- expanded_long %>%
 left_join(clone_sharing, by = "CTaa")

expanded_summary <- clone_sharing %>%
 summarise(
  total_clonotypes = n(),
  shared_clonotypes = sum(n_tissues > 1),
  unique_clonotypes = sum(n_tissues == 1)
 )
print(expanded_summary)

plot_stage_venn <- function(stage_label, color="#1f78b4") {
 
 # Filter for this stage
 stage_data <- expanded_long %>%
  filter(Type == stage_label & shared_status == "Shared_between_tissues")
 
 if (nrow(stage_data) == 0) {
  message("⚠️ No shared clonotypes found for stage: ", stage_label)
  return(NULL)
 }
 
 # Unique clonotypes per stage
 n_stage_clonotypes <- n_distinct(stage_data$CTaa)
 
 # Create list of clonotypes per tissue
 tissue_lists <- stage_data %>%
  group_by(tissues) %>%
  summarise(clones = list(unique(CTaa)), .groups = "drop") %>%
  deframe()
 
 # Print per-tissue clonotype counts
 cat("\nStage:", stage_label, "\n")
 print(stage_data %>% group_by(tissues) %>% summarise(n_clones = n_distinct(CTaa)))
 
 # Plot Venn
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
   "Shared expanded clonotypes – ",
   stage_label,
   " (n=", n_stage_clonotypes, ")"
  ))
 
 print(ggp)
 return(ggp)
}

venn_fib     <- plot_stage_venn("Fibrosis", color="#e41a1c")
venn_no_fib  <- plot_stage_venn("No_fibrosis", color="#377eb8")
venn_healthy <- plot_stage_venn("Healthy", color="#4daf4a")

# Step 5: Save plots
out_dir <- "/data/home/hdx044/plots/screpertoire"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "Shared_Venn_Fibrosis.svg"), plot = venn_fib, width = 6, height = 6)
ggsave(file.path(out_dir, "Shared_Venn_No_Fibrosis.svg"), plot = venn_no_fib, width = 6, height = 6)
ggsave(file.path(out_dir, "Shared_Venn_Healthy.svg"), plot = venn_healthy, width = 6, height = 6)

cat("✅ Saved 3 SVG plots in:", out_dir, "\n")

#### Venn diagram hyperexpanded clonotypes (95) ####

library(dplyr)
library(tidyr)
library(ggVennDiagram)


# Step 1: Filter hyperexpanded clones

hyperexpanded_clones <- clone_counts %>% filter(n_cells > 20)


# Step 2: Prepare long-format table for hyperexpanded clones

hyperexpanded_long <- expanded_long %>%
 filter(CTaa %in% hyperexpanded_clones$CTaa) %>%
 select(CTaa, Type, tissues) %>%
 distinct()  # one row per clone per tissue per stage


# Step 3: Venn plotting function per stage

plot_stage_venn <- function(stage_label, color="#1f78b4") {
 
 # Filter for stage
 stage_data <- hyperexpanded_long %>%
  filter(Type == stage_label)
 
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


# Step 4: Plot Venn diagrams for each stage

venn_fib     <- plot_stage_venn("Fibrosis", color="#e41a1c")
venn_no_fib  <- plot_stage_venn("No_fibrosis", color="#377eb8")

# Step 5: Save plots
out_dir <- "/data/home/hdx044/plots/screpertoire"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "Hyperexpanded_Venn_Fibrosis.svg"), plot = venn_fib, width = 6, height = 6)
ggsave(file.path(out_dir, "Hyperexpanded_Venn_No_Fibrosis.svg"), plot = venn_no_fib, width = 6, height = 6)

cat("✅ Saved 3 SVG plots in:", out_dir, "\n")

#### Find FNA clones in other tissues ####

# Check that combined.TCR is a named list of samples
names(combined.TCR)

# Check metadata for each sample
lapply(combined.TCR, function(df) colnames(df))

# Collapse the list into one long table
# Ensure combined.TCR is a list of data.frames
stopifnot(is.list(combined.TCR))
stopifnot(all(vapply(combined.TCR, is.data.frame, logical(1))))

# Bind all with the list name kept as SampleName
all_tcr <- bind_rows(
 imap(combined.TCR, ~ mutate(.x, SampleName = .y))
)

all_tcr <- all_tcr %>%
 mutate(
  Donor  = str_match(SampleName, "^GC-WL-([0-9A-Za-z\\-]+)")[,2],
  Tissue = str_extract(SampleName, "[A-Z]+$")    # last capital-word
 )

head(all_tcr %>% select(SampleName, Donor, Tissue, CTaa))


## Build liver_df from your pasted 3-column list
# Paste your entire list into txt, including the header line:
txt <- "
Donor\tTissue\tCTaa
10738\tLIVER\tNA_CATSRESGSHGEQFF
10738\tLIVER\tNA_CSAKPGPSGYYEQYF
10738\tLIVER\tCATDAYSSYKLIF_CASRTPGNEQFF
10738\tLIVER\tCAVGPPSGGSYIPTF_CAISESGTHVGYNEQFF
10738\tLIVER\tCARPSQYSGAGSYQLTF_CSAKPGPSGYYEQYF
10738\tLIVER\tCAAVGRGGGNKLTF_CASIPRDRGRKNYGYTF
10738\tLIVER\tCAGAPLDSNYQLIW_NA
10738\tLIVER\tCVVSSRMDSSYKLIF_CASSPGTGIYNEQFF
10738\tLIVER\tCILAGDMRF_CSAKPGPSGYYEQYF
10738\tLIVER\tCAGAPLDSNYQLIW_CATSRESGSHGEQFF
10738\tLIVER\tCILAGDMRF;CARPSQYSGAGSYQLTF_CSAKPGPSGYYEQYF
10738\tLIVER\tNA_CASRTPGNEQFF
10738\tLIVER\tNA_CSARDRLMGKNIQYF
10738\tLIVER\tNA_CASSPGTGIYNEQFF
10738\tLIVER\tCALSPNTGNQFYF_CASSWLSGTGNTIYF
10738\tLIVER\tCATDAYSSYKLIF_CASRRDGSPMDGYTF
10738\tLIVER\tCAVLSPDKLIF_NA
10738\tLIVER\tCAVMDSNYQLIW_CASRTPGNEQFF
10738\tLIVER\tCAVSLGGGKLIF_NA
10738\tLIVER\tNA_CAISESGTHVGYNEQFF
10738\tLIVER\tNA_CASGPGTSRGYTF
10738\tLIVER\tNA_CASRPPRGLPYEQYF
10738\tLIVER\tNA_CASRRDGSPMDGYTF
10738\tLIVER\tNA_CASSEPQDYLNYEQYF
10738\tLIVER\tNA_CASSFGGRGYEQYF
10738\tLIVER\tNA_CASSHRLAGRTDTQYF
10738\tLIVER\tNA_CASSLGESSYGYTF
10738\tLIVER\tNA_CASSYLDREGGYTF
10738\tLIVER\tNA_CASSYSKPRTRANEQFF
10738\tLIVER\tNA_CSARVYTEAFF
11040\tLIVER\tCLVGSNYGGSQGNLIF_CASSHPTRGSYNEQFF
11040\tLIVER\tCAGPNNAGNMLTF_CASSFRRSGLREKLFF
11040\tLIVER\tCAVMDSNYQLIW_CASSDSTGGGYTF
11040\tLIVER\tCAVNMWTGGFKTIF_CASSERPGRGGYTF
11040\tLIVER\tCAVRGPGAGSYQLTF_CASSIGRRASGYTF
11040\tLIVER\tNA_CSASGRAFGTEAFF
11040\tLIVER\tCAAMDSNYQLIW_CSARGVDDRALPNQETQYF
11040\tLIVER\tCAGRTPGGATNKLIF_CASSRDRVHKLFF
11040\tLIVER\tCALLDSNYQLIW_CASSTTPGQGAAEQYF
11040\tLIVER\tCALSDTMTTDSWGKLQF_CASSEALEAFF
11040\tLIVER\tCAMRETRSNDYKLSF;CALSVKGNNRLAF_CSATTGQGITEAFF
11040\tLIVER\tCAMSKSGGGADGLTF_CASSSRAVTNTQYF
11040\tLIVER\tCASRKQLTF;CAVRDSNYQLIW_CSARGEGDSQPQHF
11040\tLIVER\tCATDTGGGSNYKLTF_CASSQDRGLATTYNEQFF
11040\tLIVER\tCAVMDSNYQLIW_CASSDTSGSYNEQFF
11040\tLIVER\tCGTETGSWGKLQF_CASSLSGTKGDTGELFF
11040\tLIVER\tCIVRVAGGETSGSRLTF_CASTTRAYGGGPGYTF
11183\tLIVER\tCACRSTNAGKSTF_NA
11183\tLIVER\tCAARIQGAQKLVF_CAWSVQGDQPQHF
11183\tLIVER\tNA_CSARLGGGTDTQYF
11183\tLIVER\tCAARLQGAQKLVF_CAWSVQGDQPQHF
11183\tLIVER\tCAGSMEYGNKLVF_CASSIWTIRGGDTQYF
11183\tLIVER\tCAVMDSNYQLIW_CASTLAETSTDTQYF
11183\tLIVER\tCAARTQGAQKLVF_CAWSVQGDQPQHF
11183\tLIVER\tCAASEDYQLIW_CASSSMSGGGTIGYTF
11183\tLIVER\tCAASSGYALNF_CASSYSGTGNTGELFF
11183\tLIVER\tCALSGYNAGNMLTF_CASCTGGTVYGYTF
11183\tLIVER\tCAVLYGNNRLAF_CASSFHPGTNSPLHF
11183\tLIVER\tCAVTGNQFYF_CASSLGGLLYNEQFF
11183\tLIVER\tCAYRSPGSSNTGKLIF_CASSQNSYNSPLHF
11183\tLIVER\tCGFIGNQGGKLIF_CASSYGAALYEQFF
11183\tLIVER\tNA_CASSPGAREKLFF
11183\tLIVER\tNA_CASSYSGTGNTGELFF
11183\tLIVER\tCAVDTGGFKTIF_CASSLYRTTDTQYF
11183\tLIVER\tCAVIDSSYKLIF_CASSQEQGSGANVLTF
11183\tLIVER\tCAVKDSNYQLIW_CASSPVGESGELFF
11183\tLIVER\tCAVQAPGGKLIF_CASSTDGNYGYTF
11183\tLIVER\tCAVVDSNYQLIW_CSASPPSVTYNEQFF
11183\tLIVER\tCAYRSARGNYQLIW_CASSLLIWGAHEQYF
11183\tLIVER\tNA_CASSIGGQIYTEAFF
11183\tLIVER\tNA_CASSLGGLLYNEQFF
11183\tLIVER\tNA_CASSPARETNTGELFF
11183\tLIVER\tNA_CASSVTVRGSYEQYF
11183\tLIVER\tNA_CASTPGHKYEQYF
11183\tLIVER\tCACRSTNAGKSTF_CASSAGGSGANVLTF
11183\tLIVER\tNA_CAWSVQGDQPQHF
11183\tLIVER\tCALESSTLTF_CASSLGLAGTYEQYF
11183\tLIVER\tCATDEDSWGKLQF_CSARKERGGLNEQFF
11183\tLIVER\tNA_CASSAGGSGANVLTF
11183\tLIVER\tNA_CAWSRQGDQPQHF
11183\tLIVER\tCAVGFSDGQKLLF_CASSLVGATYEQYF
11183\tLIVER\tCAMSPISNFGNEKLTF_CASSFSGTWGGIPYEQYF
11183\tLIVER\tCATYNYGQNFVF_CASSVTVRGSYEQYF
11183\tLIVER\tCAVKDSNYQLIW_CASSHGQGGANVLTF
11183\tLIVER\tNA_CASSPRGSGDTDTQYF
11183\tLIVER\tCAAQTGNQFYF_CSVGQGGPADTQYF
11183\tLIVER\tCAASGGFNKFYF_CASSVGVADGYTF
11183\tLIVER\tCAAYNTDKLIF_CASSLGAGYTF
11183\tLIVER\tCAGQGEYGGSQGNLIF_CASSYGGAYGYTF
11183\tLIVER\tCAGYNQGGKLIF_NA
11183\tLIVER\tCALRRANSGYALNF_CASSTYSGTVIF
11183\tLIVER\tCASMDSNYQLIW_CAWNAGDTGELFF
11183\tLIVER\tCAVLDSNYQLIW_CASTDREDTEAFF
11183\tLIVER\tCAVLDSNYQLIW_NA
11183\tLIVER\tCAVMDSNYQLIW_CASSERDTGEQYF
11183\tLIVER\tCAVMDSNYQLIW_CASSFGGGHTDTQYF
11183\tLIVER\tCAVMDSNYQLIW_CSARGPGQEETQYF
11183\tLIVER\tCAVRDRNNYGQNFVF_NA
11183\tLIVER\tCAVRDSNYQLIW_CASSQDPGQGADEQYF
11183\tLIVER\tCAVRDSNYQLIW_NA
11183\tLIVER\tCAVRDTSNAGNMLTF_CATSSPSSEGLIGELFF
11183\tLIVER\tCAVVDSNYQLIW_CASSPRGSGDTDTQYF
11183\tLIVER\tCVPWRYKLSF_CASSFPGSSYNEQFF
11183\tLIVER\tCVVNGRGSSASKIIF_NA
11183\tLIVER\tNA_CASTLAETSTDTQYF
10291-1\tLIVER\tCAENLNTGFQKLVF_CASSLGYEQYF
10291-1\tLIVER\tCAYRTPSGSARQLTF_CSASPVAANSPLHF
10291-1\tLIVER\tCIVRVHTGNQFYF_CASSQTVLDGETQYF
10291-1\tLIVER\tNA_CASRRQETYNEQFF
10291-1\tLIVER\tCATLRNSGNTPLVF_CASTPFGTGELFF
10291-1\tLIVER\tCAASGSSNTGKLIF_CASSAVSRGSIGNTIYF
10291-1\tLIVER\tCAGEGSNFGNEKLTF_CASRDNYEQYF
10291-1\tLIVER\tCALDMRSNYQLIW_CASSLVRTNTGELFF
10291-1\tLIVER\tCALHSGGYQKVTF_CASSEKDRFDNEQFF
10291-1\tLIVER\tNA_CASSLVRTNTGELFF
10291-1\tLIVER\tCALANTNAGKSTF_CASRRSTSGNTIYF
10291-1\tLIVER\tCAVGFGNVLHC_CASRTGYYGYTF
10291-1\tLIVER\tNA_CSASPVAANSPLHF
10291-1\tLIVER\tCAANYGQNFVF_CASGSHGQGETDTQYF
10291-1\tLIVER\tCAFMKPGEGAQKLVF_CASHPLGFYSNQPQHF
10291-1\tLIVER\tCAGMDSNYQLIW_CASSDSSTDTQYF
10291-1\tLIVER\tCALSVMDSSYKLIF_CAIELAVAGELFF
10291-1\tLIVER\tCAVMDSNYQLIW_NA
10291-1\tLIVER\tCGTVNNDMRF_CASSVGTGSNTEAFF
10291-1\tLIVER\tNA_CASSLFGDAGTATQYF
10291-1\tLIVER\tCAANKLTF_CASSLVAGRETQYF
10291-1\tLIVER\tCAASITLGAQKLVF_CASSLDGATQYF
10291-1\tLIVER\tCAENLNTGFQKLVF_NA
10291-1\tLIVER\tCAGGRRALTF_CASSPGTSGSTDTQYF
10291-1\tLIVER\tCAVRDSNYQLIW_NA
10291-1\tLIVER\tCAEGTSYGKLTF_CASSKPPRGRNTIYF
10291-1\tLIVER\tCAGMDSNYQLIW_NA
10291-1\tLIVER\tCAMREPARHNNDMRF;CAAPGGSYIPTF_CATSATGRKDQPQHF
10291-1\tLIVER\tCAMSGAGGATNKLIF_CAISEGADGYTF
10291-1\tLIVER\tCATDARSGNTPLVF_NA
10291-1\tLIVER\tCAVDMNYGGSQGNLIF_CASSESGQGYEQYF
10291-1\tLIVER\tCAVMDSNYQLIW_CASRTGDTGELFF
10291-1\tLIVER\tCAYRSARRTPLVF_CASSLVRLAGGPTGELFF
10291-1\tLIVER\tCVVRGPGRRALTF_CASSQDTSLTNEKLFF
10291-1\tLIVER\tCATDALLSGATNKLIF_CASRRQETYNEQFF
10291-1\tLIVER\tCATDAWGSQGNLIF_CASRSSGQETQYF
10291-1\tLIVER\tCAESQEEYGNKLVF_CASSRRQGSNIQYF
10291-1\tLIVER\tCAMGGDNNNDMRF_CASSYRRGNTEAFF
10291-1\tLIVER\tCAVQGDNQGGKLIF_CAWSVVDNEQFF
10291-1\tLIVER\tCLVGGAGSYQLTF_CASSLFGDAGTATQYF
10291-1\tLIVER\tCAPLGDDKIIF_CASIRDKGTDTQYF
10291-1\tLIVER\tCATDARSGNTPLVF_CASSLRGQGVRSPLHF
10291-1\tLIVER\tCATDDWGSQGNLIF_CASRRAGGDEQYF
10291-1\tLIVER\tCAVFSNTGKLIF_CASSSLESYEQYF
10291-1\tLIVER\tNA_CASSLGTSTGELFF
10291-1\tLIVER\tCAVVDSNYQLIW_CASSPGQSTDTQYF
10291-1\tLIVER\tNA_CASSQTVLDGETQYF
10291-1\tLIVER\tNA_CASRDNYEQYF
10291-1\tLIVER\tNA_CASSLGYEQYF
10291-1\tLIVER\tNA_CASTPFGTGELFF
10291-1\tLIVER\tCAEGLVNDMRF_CASSLRSGAKNIQYF
10291-1\tLIVER\tNA_CASRSSGQETQYF
10291-1\tLIVER\tNA_CASSLDGATQYF
10291-1\tLIVER\tNA_CASSRRQGSNIQYF
10291-1\tLIVER\tCAALDSNYQLIW_CASSQSSEVHEQYF
10291-1\tLIVER\tCALDSGGSNYKLTF_CASTPFGTGELFF;CASSPSGAPDKYF
10291-1\tLIVER\tCAMSGGTGGSYIPTF_CASSPAAGEDQETQYF
10291-1\tLIVER\tCAPMDSNYQLIW_CASSDSGGAPYNEQFF
10291-1\tLIVER\tCAVLDSNYQLIW_CSAISPSGSNSLETQYF
10291-1\tLIVER\tCAVMDSNYQLIW_CSARQGENEQFF
10291-1\tLIVER\tCAVRDSNYQLIW_CASSDSSGANVLTF
10291-1\tLIVER\tCLVGDFDGYNKLIF_CASSPASGSYNSPLHF
10291-1\tLIVER\tCVVKAAGNKLTF_CASSYAESYGYTF
10291-1\tLIVER\tNA_CASADFSNQPQHF
10291-1\tLIVER\tNA_CASGSGEHQPQHF
10291-1\tLIVER\tNA_CASNRGTFTDTQYF
10291-1\tLIVER\tNA_CASRTGDTGELFF
10291-1\tLIVER\tNA_CASSFRRGSYEQYF
10291-1\tLIVER\tNA_CASSLDWGGGETQYF
10291-1\tLIVER\tNA_CASSLINEQFF
10291-1\tLIVER\tNA_CATSPRPGPGYNEQFF
10291-1\tLIVER\tNA_CSAEGLAGYYEQYF
"

# Parse tab-separated block -> data.frame
liver_df <- read.delim(textConnection(txt), header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")

# Basic cleanup
liver_df <- liver_df %>%
 mutate(
  Donor  = trimws(as.character(Donor)),
  Tissue = trimws(as.character(Tissue)),
  CTaa   = trimws(as.character(CTaa))
 )

# Remove any blank lines / rows with missing CTaa
liver_df <- liver_df %>% filter(!is.na(CTaa), CTaa != "")

# Vector of unique CTaa
liver_clones <- unique(liver_df$CTaa)

# Quick checks
nrow(liver_df)             # total rows pasted
length(liver_clones)       # unique CTaa; you expect 177
length(liver_clones) == 177

# Which of the 177 are in any tissue?
shared_clones <- all_tcr %>%
 filter(CTaa %in% liver_clones) %>%
 select(Donor, Tissue, SampleName, CTaa)

# Clone → tissues
clone_where <- shared_clones %>%
 group_by(CTaa) %>%
 summarise(
  Donors  = paste(sort(unique(Donor)), collapse = ", "),
  Tissues = paste(sort(unique(Tissue)), collapse = ", "),
  Samples = paste(sort(unique(SampleName)), collapse = "; "),
  .groups = "drop"
 )

# Liver-unique vs shared outside liver
clone_liver_only <- clone_where %>% filter(Tissues == "LIVER")
clone_shared     <- clone_where %>% filter(Tissues != "LIVER")

# Quick counts
cat("Unique CTaa provided: ", length(liver_clones), "\n")
cat("Found in data:        ", n_distinct(shared_clones$CTaa), "\n")
cat("Also outside liver:   ", n_distinct(shared_clones %>% filter(Tissue != "LIVER") %>% pull(CTaa)), "\n")
# Also outside liver:    53 

clone_outside_tissue_map <- shared_clones %>%
 filter(Tissue != "LIVER") %>%
 group_by(CTaa) %>%
 summarise(
  Donors  = paste(sort(unique(Donor)), collapse = ", "),
  Tissues = paste(sort(unique(Tissue)), collapse = ", "),
  Samples = paste(sort(unique(SampleName)), collapse = "; "),
  .groups = "drop"
 )

print (clone_outside_tissue_map, n=53)

write_csv(clone_outside_tissue_map, "/data/home/hdx044/files/screpertoire/FNA/liver/sharedExpandedClones.csv")















