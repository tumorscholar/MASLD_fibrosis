
library(Seurat)
library(scRepertoire)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)


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
                                         "No_Fibrosis","No_Fibrosis","No_Fibrosis","No_Fibrosis",
                                         "No_Fibrosis","No_Fibrosis","No_Fibrosis","No_Fibrosis",
                                         "No_Fibrosis","No_Fibrosis","No_Fibrosis","No_Fibrosis",
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
length(combined.TCR)    # 80, matches

# ️Assign tissue info to each sample in the list
for (i in seq_along(combined.TCR)) {
  combined.TCR[[i]]$Tissue <- tissue_vector[i]
}

head(combined.TCR[[1]])

# Get all unique patient IDs
patient_ids <- c("GC-WL-10742", "GC-WL-9961", "GC-WL-9999", "GC-WL-10113-2",
                 "GC-WL-9680", "GC-WL-10203", "GC-WL-10113-1", "GC-WL-10380",
                 "GC-WL-10291-2", "GC-WL-10202", "GC-WL-10634", "GC-WL-10738",
                 "GC-WL-10205", "GC-WL-9991", "GC-WL-9932", "GC-WL-11040",
                 "GC-WL-11051", "GC-WL-11183", "GC-WL-11471", "GC-WL-11327")

# Generate plots for all patients
plot_list <- lapply(patient_ids, function(patient_id) {
  
  # Get samples for this patient
  patient_samples <- unique(unlist(
    lapply(combined.TCR, function(df) {
      df$sample[grepl(patient_id, df$sample)]
    })
  ))
  
  # Get stage for this patient
  patient_stage <- combined.TCR[[which(sapply(combined.TCR, function(df) 
    any(grepl(patient_id, df$sample))))[1]]]$Type[1]
  
  # Generate alluvial plot
  p <- clonalCompare(
    combined.TCR,
    top.clones = 5,
    samples = patient_samples,
    cloneCall = "aa",
    graph = "alluvial"
  ) +
    labs(
      title = paste0(patient_id, " | ", patient_stage),
      fill = "Clonotypes"
    ) +
    scale_x_discrete(labels = c("LIVER", "PBMC", "SAT", "VAT")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.text.x = element_text(size = 9),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      panel.border = element_blank(),        # removes box around each individual plot
      plot.margin = unit(c(0.3, 0.8, 0.3, 0.3), "cm")  # top, RIGHT, bottom, left — increase right for more width
    )
  
  return(p)
})

# Layout settings
ncol_fig <- 4                                
nrow_fig <- ceiling(length(plot_list) / ncol_fig)

# Width per column in inches — increase this number for wider individual plots
width_per_col <- 6    
height_per_row <- 5   

# Arrange plots in grid
final_plot <- arrangeGrob(
  grobs = plot_list,
  ncol = ncol_fig,
  nrow = nrow_fig
)

# Save as TIFF
ggsave(
  filename = "/data/home/hdx044/plots/screpertoire/alluvial_top5_clones_all_patients.tiff",
  plot = final_plot,
  width = ncol_fig * 6,       
  height = nrow_fig * 6,  
  units = "in",
  dpi = 300,
  compression = "lzw" 
)

# Save as PDF
ggsave(
  filename = "/data/home/hdx044/plots/screpertoire/alluvial_top5_clones_all_patients.pdf",
  plot = final_plot,
  width = ncol_fig * 6,
  height = nrow_fig * 6,
  units = "in"
)

# End of the script
