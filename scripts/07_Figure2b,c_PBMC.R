###############################################################################
# Script name: 06_Figure2b,c_PBMC.R
#
# Description:
# This script performs T‑cell receptor (TCR) repertoire analysis for PBMC
# samples using scRepertoire. It integrates TCR contig annotation files from
# multiple patients across Healthy, No Fibrosis, and Fibrosis groups and
# generates clonotype‑level summaries and diversity metrics.
#
# The script reads per‑sample TCR contig CSV files derived from single‑cell
# V(D)J sequencing, combines them into a unified TCR repertoire object, and
# annotates each sample with disease status (Healthy, No Fibrosis, Fibrosis).
#
# Key steps performed in this script include:
#  - Loading PBMC‑derived TCR contig annotation files
#  - Combining TCR repertoires across samples using scRepertoire::combineTCR
#  - Annotating samples with fibrosis status
#  - Exporting clonotype information for downstream analysis
#  - Quantifying clonal expansion using clonalQuant
#  - Assessing TCR diversity using Shannon and inverse Simpson indices
#  - Generating publication‑ready figures summarizing clonotype proportions
#    and diversity across disease groups
#
# All analyses in this script are restricted to PBMC samples and focus on
# comparing TCR clonality and diversity across fibrosis stages.
#
# Outputs:
#  - CSV file containing exported TCR clonotypes
#  - SVG figures showing clonal expansion and diversity metrics
#
# Purpose:
# To characterise PBMC‑specific TCR repertoire structure and diversity across
# fibrosis stages in MASLD.
#

###############################################################################

library(scRepertoire)
library(patchwork)
library(tidyverse)

setwd ('/data/home/hdx044/files/screpertoire/demux_contig/TCR')

#Read files
#Healthy sample
s1 = read.csv("GC-WL-10742-SAT-VAT-PBMC_PBMC_TCR_contig.csv")

#F0 samples
s2 = read.csv("GC-WL-9961-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s3 = read.csv("GC-WL-9999-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s4 = read.csv("GC-WL-10113-2-SAT-VAT-PBMC_PBMC_TCR_contig.csv")

#F1 samples
s5 = read.csv("GC-WL-9680-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s6 = read.csv("GC-WL-10203-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s7 = read.csv("GC-WL-10113-1-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s8 = read.csv("GC-WL-10380-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s9 = read.csv("GC-WL-10291-2-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s10 = read.csv("GC-WL-10202-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s11 = read.csv("GC-WL-10634-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s12 = read.csv("GC-WL-10738-SAT-VAT-PBMC_PBMC_TCR_contig.csv")

#F2 samples
s13 = read.csv("GC-WL-10205-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s14 = read.csv("GC-WL-9991-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s15 = read.csv("GC-WL-9932-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s16 = read.csv("GC-WL-11040-SAT-VAT-PBMC_PBMC_TCR_contig.csv")

#F3 samples
s17 = read.csv("GC-WL-11051-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s18 = read.csv("GC-WL-11183-SAT-VAT-PBMC_PBMC_TCR_contig.csv")
s19 = read.csv("GC-WL-11471-SAT-VAT-PBMC_PBMC_TCR_contig.csv")

# F1 samples
s20 = read.csv("GC-WL-11327-SAT-VAT-PBMC_PBMC_TCR_contig.csv")


#list
contig_list = list(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20)

combined.TCR = combineTCR(contig_list, samples = c("10742",
                                                   "9961",
                                                   "9999",
                                                   "10113-2",
                                                   "9680",
                                                   "10203",
                                                   "10113-1",
                                                   "10380",
                                                   "10291-2",
                                                   "10202",
                                                   "10634", 
                                                   "10738", 
                                                   "10205",
                                                   "9991",
                                                   "9932",
                                                   "11040",
                                                   "11051",
                                                   "11183",
                                                   "11471",
                                                   "11327"),
                          removeNA = FALSE, 
                          removeMulti = FALSE,
                          filterMulti = FALSE)

combined.TCR = addVariable(combined.TCR, 
                           variable.name = "Type", 
                           variables = c("Healthy",
                                         "No_fibrosis",
                                         "No_fibrosis",
                                         "No_fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis"))

head(combined.TCR[[1]])

for (i in seq_along(combined.TCR)) {
 combined.TCR[[i]]$Type <- factor(combined.TCR[[i]]$Type,
                                  levels = c("Healthy", "No_fibrosis", "Fibrosis"))
}

setwd ("/data/home/hdx044/files/screpertoire")

exportClones(combined.TCR, 
             write.file = TRUE,
             file.name = "allclonesPBMC.csv")


#### Fig2.b ####
ggp <- clonalQuant(combined.TCR,
                   cloneCall = "aa",
                   chain = "both",
                   group.by = "Type",
                   scale = TRUE,
                   order.by = c("Healthy", "No_fibrosis", "Fibrosis")) +
 scale_fill_manual(values = c(
  "Healthy" = "#88CCEE",       # light blue
  "No_fibrosis" = "#DDCC77",   # yellow
  "Fibrosis" = "#CC6677"       # reddish pink
 )) +
 scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +  # Tick every 10
 ggtitle("PBMC") + 
 theme_classic(base_size = 12) +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.line = element_line(color = "black"),
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  legend.key.size = unit(1.5, "lines")
 )

ggp

ggsave(
 filename = "uniqueClonePercentPBMC.svg",
 plot = ggp,              
 device = "svg",
 path = "~/plots/screpertoire/PBMC",
 width = 3,                      
 height = 5,                     
 units = "in",
 dpi = 600                        
)

#### Fig2.c ####
# Force factor level order across combined.TCR list
combined.TCR <- lapply(combined.TCR, function(df) {
  df$Type <- factor(df$Type, levels = c("Healthy", "No_fibrosis", "Fibrosis"))
  return(df)
})

# Calculate Shannon diversity
p_shannon <- clonalDiversity(
  combined.TCR, 
  cloneCall = "aa",
  group.by = "Type",
  metric = "shannon"
) +
  scale_fill_manual(
    name = "Type",
    values = c(
      "Healthy" = "#88CCEE",
      "No_fibrosis" = "#DDCC77",
      "Fibrosis" = "#CC6677"
    ),
    breaks = c("Healthy", "No_fibrosis", "Fibrosis")
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  labs(title = "Shannon Index", y = "Shannon Diversity")

# Calculate Inverse Simpson diversity
p_inv_simpson <- clonalDiversity(
  combined.TCR, 
  cloneCall = "aa",
  group.by = "Type",
  metric = "inv.simpson"
) +
  scale_fill_manual(
    name = "Type",
    values = c(
      "Healthy" = "#88CCEE",
      "No_fibrosis" = "#DDCC77",
      "Fibrosis" = "#CC6677"
    ),
    breaks = c("Healthy", "No_fibrosis", "Fibrosis")
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.5, "lines")
  ) +
  labs(title = "Inverse Simpson Index", y = "Inverse Simpson Diversity")

# Combine side by side
combined_diversity <- p_shannon | p_inv_simpson

combined_diversity <- combined_diversity +
  plot_annotation(
    title = "PBMC TCR Diversity",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

print(combined_diversity)

ggsave(
 filename = "clonalLengthPBMC.svg",
 plot = combined_diversity,              
 device = "svg",
 path = "~/plots/screpertoire/PBMC",
 width = 5,                      
 height = 5,                     
 units = "in",
 dpi = 600                        
)


# End of the script
