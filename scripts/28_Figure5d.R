####################################################################################
# Script name: 28_Figure5d.R
#
#   Liver FNA TCR Repertoire Analysis
# 
# - scRepertoire-based TCR analysis
# - Baseline vs follow-up liver FNAs
# - Top-100 clonal dominance
# - Baseline-expanded clonotype tracking
#
#####################################################################################     

library(scRepertoire)  
library(dplyr)         
library(tidyr)         
library(ggplot2)  

setwd('/data/home/hdx044/files/screpertoire/demux_contig/TCR')

#Read TCR contig files
s3 = read.csv("GC-WL-10738-LIVER_LIVER_TCR_contig.csv")          # 10738 LIVER baseline
s4 = read.csv("GC-WL-11570-LIVER_TCR_contig.csv")                # 10738 LIVER followup

s6 = read.csv("GC-WL-10291-1-LIVER_TCR_contig.csv")              # 10291-1 LIVER baseline
s7 = read.csv("GC-WL-11303-LIVER_TCR_contig.csv")                # 10291-1 LIVER followup

s10 = read.csv("GC-WL-11040-LIVER_LIVER_TCR_contig.csv")         # 11040 LIVER baseline
s11 = read.csv("GC-WL-11816-LIVER_TCR_contig.csv")               # 11040 LIVER followup

s14 = read.csv("GC-WL-11183-LIVER_LIVER_TCR_contig.csv")         # 11183 LIVER baseline
s15 = read.csv("GC-WL-11937-LIVER_TCR_contig.csv")               # 11183 LIVER followup


# List of TCR data 
contig_list = list(s3,s4,s6,s7,s10,s11,s14,s15)

#Correct sample names (aligned with your real metadata)
sample_names <- c(
 "Baseline-10738-LIVER",
 "Followup-10738-LIVER",
 
 "Baseline-10291-1-LIVER",
 "Followup-10291-1-LIVER",
 
 "Baseline-11040-LIVER", 
 "Followup-11040-LIVER",
 
 "Baseline-11183-LIVER", 
 "Followup-11183-LIVER"
)

# Combine
combined.TCR <- combineTCR(
 contig_list,
 samples = sample_names,
 removeNA = FALSE,
 removeMulti = FALSE,
 filterMulti = FALSE
)

# Add Type correctly
Type_vector <- c(
 "Baseline-LIVER",  # 10738 liver baseline
 "Followup-LIVER",  # 10738 liver followup
 
 "Baseline-LIVER",  # 10291 liver baseline
 "Followup-LIVER",  # 10291 liver followup
 
 "Baseline-LIVER",  # 11040 liver baseline
 "Followup-LIVER",  # 11040 liver followup
 
 "Baseline-LIVER",  # 11183 liver baseline
 "Followup-LIVER"   # 11183 liver followup
)

combined.TCR <- addVariable(
 combined.TCR,
 variable.name = "Type",
 variables = Type_vector
)

head(combined.TCR[[1]])

# Order the Type factor for plotting
combined.TCR$Type <- factor(
 combined.TCR$Type,
 levels = c(
  "Baseline-LIVER",
  "Followup-LIVER"
 )
)

combined.TCR$Type <- NULL

for (i in seq_along(combined.TCR)) {
 combined.TCR[[i]]$Type <- Type_vector[i]
}

for (i in seq_along(combined.TCR)) {
 combined.TCR[[i]]$Type <- factor(
  combined.TCR[[i]]$Type,
  levels = c( "Baseline-LIVER", "Followup-LIVER")
 )
}

table(unlist(lapply(combined.TCR, function(x) x$Type)))

head(combined.TCR[[1]])

# Count total T cells
total_tcells <- sum(sapply(combined.TCR, nrow))
cat("Total number of T cells:", total_tcells, "\n")
# Total number of T cells: 8952

# Combine all samples into one dataframe
all_tcells <- dplyr::bind_rows(combined.TCR, .id = "Sample")

# Count unique clonotypes
n_unique_clonotypes <- length(unique(all_tcells$CTaa))
cat("Number of unique clonotypes:", n_unique_clonotypes, "\n")
# Number of unique clonotypes: 7160

# Count cells per clonotype
clone_counts <- all_tcells %>%
 group_by(CTaa) %>%
 summarise(n_cells = n(), .groups = "drop")

# Separate expanded and singleton clones
expanded_clones <- clone_counts %>% filter(n_cells > 1)
singleton_clones <- clone_counts %>% filter(n_cells == 1)

n_expanded_cells <- sum(expanded_clones$n_cells)
n_singleton_cells <- sum(singleton_clones$n_cells)

cat("Expanded T cells:", n_expanded_cells, "\n")
#Expanded T cells: 2305 

cat("Singleton T cells:", n_singleton_cells, "\n")
#Singleton T cells: 6647 

cat("Overall expansion percentage:", round((n_expanded_cells/total_tcells)*100, 2), "%\n")
#Overall expansion percentage: 25.75 %

# Repertoire occupancy
# Calculate Top 100 Clone Frequencies
top_clone_list <- vector("list", length(combined.TCR))

for (i in seq_along(combined.TCR)) {
  
  sample_data <- combined.TCR[[i]]
  sample_name <- names(combined.TCR)[i]
  
  # Clone frequency table
  clone_freq <- table(sample_data$CTstrict)
  clone_freq_sorted <- sort(clone_freq, decreasing = TRUE)
  
  total_cells <- sum(clone_freq)
  top100 <- (sum(clone_freq_sorted[
    1:min(100, length(clone_freq_sorted))
  ]) / total_cells) * 100
  
  # Metadata (LIVER only)
  top_clone_list[[i]] <- data.frame(
    Sample         = sample_name,
    Timepoint      = ifelse(grepl("Baseline", sample_name), "Baseline", "Followup"),
    Tissue         = "LIVER",
    Donor          = gsub("(Baseline-|Followup-)(.+)-LIVER", "\\2", sample_name),
    Top100_Percent = round(top100, 2),
    stringsAsFactors = FALSE
  )
}

# Combine results
top_clone_analysis <- bind_rows(top_clone_list)

# Safety checks
stopifnot(all(top_clone_analysis$Tissue == "LIVER"))
stopifnot(nrow(top_clone_analysis) == length(combined.TCR))

cat("\n=== TOP 100 PER SAMPLE ===\n")
print(top_clone_analysis)

# Baseline vs Followup summary
top100_comparison <- top_clone_analysis %>%
  group_by(Tissue, Timepoint) %>%
  summarise(
    Mean_Top100 = mean(Top100_Percent),
    SD_Top100   = sd(Top100_Percent),
    n_samples   = n(),
    .groups = "drop"
  )

cat("\n=== TOP 100 COMPARISON: Baseline vs Followup ===\n")
print(top100_comparison)

# Paired change per donor
top100_change <- top_clone_analysis %>%
  select(Donor, Tissue, Timepoint, Top100_Percent) %>%
  pivot_wider(names_from = Timepoint, values_from = Top100_Percent) %>%
  filter(!is.na(Baseline) & !is.na(Followup)) %>%
  mutate(
    Change = Followup - Baseline,
    Percent_Change = round(Change / Baseline * 100, 1),
    Direction = case_when(
      Change < -2 ~ "Improved ↓",
      Change >  2 ~ "Worsened ↑",
      TRUE        ~ "Stable"
    )
  )

cat("\n=== TOP 100 CHANGE BY DONOR ===\n")
print(top100_change)

# Statistical test (paired Wilcoxon)

cat("\n=== STATISTICAL TEST (LIVER, paired Wilcoxon) ===\n")

baseline_paired <- top100_change$Baseline
followup_paired <- top100_change$Followup

if (length(baseline_paired) > 1) {
  wtest <- wilcox.test(baseline_paired, followup_paired, paired = TRUE)
  
  cat(sprintf(
    "p = %.4f | Mean change = %.2f (%.1f%%)\n",
    wtest$p.value,
    mean(followup_paired - baseline_paired),
    (mean(followup_paired) - mean(baseline_paired)) /
      mean(baseline_paired) * 100
  ))
}

# Prepare data for plotting
top100_plot_data <- top_clone_analysis %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("Baseline", "Followup"))
  )

#### Fig5.d ####
p_lines <- ggplot(
  top100_plot_data,
  aes(
    x = Timepoint,
    y = Top100_Percent,
    group = Donor,
    color = Donor
  )
) +
  geom_line(size = 1.4, alpha = 0.85) +
  geom_point(size = 4, alpha = 0.95) +   # inherits same colour as line
  facet_wrap(~Tissue) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Top 100 TCR Clonal Dominance in Liver",
    subtitle = "Paired baseline vs follow‑up FNA samples",
    x = "Timepoint",
    y = "Top 100 Clones (% of repertoire)",
    color = "Donor"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(breaks = seq(0, 60, 5))

print(p_lines)

setwd("/data/home/hdx044/plots/screpertoire/liver/FNA")
svg("FNA_Top100_TCR_ClonalDominance_liver.svg", width = 6, height = 6)   
print(p_lines)
dev.off()

#### Baseline expanded clonotype at followup ####

# Identify paired LIVER samples
baseline_liver_samples <- names(combined.TCR)[grepl("^Baseline-.*-LIVER$", names(combined.TCR))]
followup_liver_samples <- names(combined.TCR)[grepl("^Followup-.*-LIVER$", names(combined.TCR))]

donors <- intersect(
  gsub("Baseline-|\\-LIVER", "", baseline_liver_samples),
  gsub("Followup-|\\-LIVER", "", followup_liver_samples)
)

# Build table per donor
liver_counts_table <- lapply(donors, function(donor) {
  
  baseline <- combined.TCR[[baseline_liver_samples[grepl(donor, baseline_liver_samples)]]]
  followup <- combined.TCR[[followup_liver_samples[grepl(donor, followup_liver_samples)]]]
  
  baseline_total <- nrow(baseline)
  followup_total <- nrow(followup)
  
  # Baseline expanded clonotypes (>=2 cells, SAME sample)
  baseline_counts <- baseline %>%
    group_by(CTaa) %>%
    summarise(Baseline_Count = n(), .groups = "drop") %>%
    filter(Baseline_Count >= 2)
  
  # Follow-up counts for same clonotypes
  followup_counts <- followup %>%
    filter(CTaa %in% baseline_counts$CTaa) %>%
    group_by(CTaa) %>%
    summarise(Followup_Count = n(), .groups = "drop")
  
  baseline_counts %>%
    left_join(followup_counts, by = "CTaa") %>%
    mutate(
      Followup_Count = replace_na(Followup_Count, 0),
      Baseline_Total_Tcells = baseline_total,
      Followup_Total_Tcells = followup_total,
      Donor = donor,
      Tissue = "LIVER"
    ) %>%
    select(
      Donor,
      Tissue,
      CTaa,
      Baseline_Count,
      Baseline_Total_Tcells,
      Followup_Count,
      Followup_Total_Tcells
    )
})

# Combine all donors
liver_baseline_followup_counts <- bind_rows(liver_counts_table)

# Sanity checks
nrow(liver_baseline_followup_counts)   # should be 177
head(liver_baseline_followup_counts)

print(liver_baseline_expanded_table, n=Inf)

write.csv(
  liver_baseline_followup_counts,
  "/data/home/hdx044/MASLD_fibrosis/files/LIVER_baseline_followup_clonotype_cell_counts.csv",
  row.names = FALSE
)

# End of the script
