################################################################################
# Script name: 31_ExtendedDataFigure4_SAT.R
#
# TCR Repertoire Gene Usage and Diversity in SAT
#
# Analyses:
#   (a) Percentage usage of TRB V and J genes
#       – stratified by disease stage
#       – Healthy, No fibrosis, Fibrosis
#
#   (b) TRB positional amino‑acid entropy
#       – comparison across disease stages
#
# Input:
#   • CellRanger TCR contig CSV files for SAT
#
# Output:
#   • Publication‑quality figures
#
###############################################################################

library(scRepertoire) 
library(ggplot2)      
library(dplyr)        
library(grid) 

setwd ('/data/home/hdx044/files/screpertoire/demux_contig/TCR')

#Read files
#Healthy sample
s1 = read.csv("GC-WL-10742-SAT-VAT-PBMC_SAT_TCR_contig.csv")

#F0 samples
s2 = read.csv("GC-WL-9961-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s3 = read.csv("GC-WL-9999-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s4 = read.csv("GC-WL-10113-2-SAT-VAT-PBMC_SAT_TCR_contig.csv")

#F1 samples
s5 = read.csv("GC-WL-9680-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s6 = read.csv("GC-WL-10203-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s7 = read.csv("GC-WL-10113-1-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s8 = read.csv("GC-WL-10380-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s9 = read.csv("GC-WL-10291-2-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s10 = read.csv("GC-WL-10202-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s11 = read.csv("GC-WL-10634-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s12 = read.csv("GC-WL-10738-SAT-VAT-PBMC_SAT_TCR_contig.csv")

#F2 samples
s13 = read.csv("GC-WL-10205-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s14 = read.csv("GC-WL-9991-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s15 = read.csv("GC-WL-9932-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s16 = read.csv("GC-WL-11040-SAT-VAT-PBMC_SAT_TCR_contig.csv")

#F3 samples
s17 = read.csv("GC-WL-11051-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s18 = read.csv("GC-WL-11183-SAT-VAT-PBMC_SAT_TCR_contig.csv")
s19 = read.csv("GC-WL-11471-SAT-VAT-PBMC_SAT_TCR_contig.csv")

# F1 samples
s20 = read.csv("GC-WL-11327-SAT-VAT-PBMC_SAT_TCR_contig.csv")


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

#### Extended Data Figure 4a (Percent V genes) ####
ggp <- percentGenes(combined.TCR, 
                    chain = "TRB",
                    gene = "Vgene",
                    group.by = "Type",
                    order.by = c("Healthy", "No_fibrosis", "Fibrosis"))+
  ggtitle("SAT") +
  ylab("Genes")+
  xlab("Stage")+
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.5, "lines")
  )

ggp

ggsave(
  filename = "percentVGenesSATStagewise.pdf",
  plot = ggp,
  device = pdf,
  path = "~/plots/screpertoire/SAT",
  width = 3,
  height = 12,
  units = "in",
  dpi = 300
)


#### Extended Data Figure 4a (Percent J genes) ####
ggp <- percentGenes(combined.TCR, 
                    chain = "TRB",
                    gene = "Jgene",
                    group.by = "Type",
                    order.by = c("Healthy", "No_fibrosis", "Fibrosis"))+
  ggtitle("SAT") +
  ylab("Genes")+
  xlab("Stage")+
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.5, "lines")
  )

ggp

ggsave(
  filename = "percentJGenesSATStagewise.pdf",
  plot = ggp,
  device = pdf,
  path = "~/plots/screpertoire/SAT",
  width = 3,
  height = 5,
  units = "in",
  dpi = 300
)

#### Extended Data Figure 4b ####
ggp <- positionalEntropy(combined.TCR, 
                         chain = "TRB",
                         group.by = "Type",
                         order.by = c("Healthy", "No_fibrosis", "Fibrosis"),
                         aa.length = 20) + 
  scale_colour_manual(name = "Type",
                      values = c(
                        "Healthy" = "#88CCEE",       # light blue
                        "No_fibrosis" = "#DDCC77",   # yellow
                        "Fibrosis" = "#CC6677"       # reddish pink
                      )) +
  ggtitle("SAT") +
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
  filename = "positionalEntropySAT.pdf",
  plot = ggp,              
  device = "pdf",
  path = "~/plots/screpertoire/SAT",
  width = 5,                      
  height = 5,                     
  units = "in",
  dpi = 600                        
)

# End of the script