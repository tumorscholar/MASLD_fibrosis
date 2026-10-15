###############################################################################
# 26 - Extended Data Figure 8 - clonal distribution of top clonotypes across tissues (alluvial)
#      Clonal distribution of the top T cell clonotypes across tissues in
#      individual patients (alluvial plot).
#
# For each patient with matched liver, PBMC, SAT and VAT (n = 20; 10291-1 is
# liver-only and is excluded), the N most abundant clonotypes (by pooled size
# across the patient's four tissues) are shown as an alluvial plot: each clonotype
# is a ribbon whose height at each tissue is its proportion of that tissue's
# repertoire. One panel per patient, each with its own clonotype legend, combined
# into a single figure.
#
# Source: v7_tcr_TR captures pooled by combineTCR with the SVP library
# demultiplexed to SAT/VAT/PBMC (identical loader to the Figure 2/3 scripts).
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
suppressMessages({
  library(scRepertoire); library(tidyverse); library(ggalluvial)
  library(patchwork); library(viridisLite)
})
source("00_nature_style.R")

tcr_dir     <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
barcode_dir <- "/data/Blizard-AlazawiLab/rk/barcodes_v7"
fig_dir     <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/ExtendedData"
tab_dir     <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/ExtendedData"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

N_TOP <- 15                                              # top clonotypes per patient to display

# ---- cohort (20 patients with liver + PBMC + SAT + VAT) + stage -------------
liver_pids <- c("9680","9932","9961","9991","9999","10113-1","10113-2","10202","10203","10205",
                "10291-1","10291-2","10380","10634","10738","10742","11040","11051","11183","11327","11471")
matched_pids <- setdiff(liver_pids, "10291-1")           # 20 matched-tissue patients
stage_map  <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No_Fibrosis","9991"="Fibrosis","9999"="No_Fibrosis",
                "10113-1"="Fibrosis","10113-2"="No_Fibrosis","10202"="Fibrosis","10203"="Fibrosis","10205"="Fibrosis",
                "10291-1"="Fibrosis","10291-2"="Fibrosis","10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis",
                "10742"="Healthy","11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis")

bfiles <- list.files(barcode_dir, pattern = "_metadata\\.csv$", full.names = TRUE)
gex_map <- do.call(rbind, lapply(bfiles, read.csv, stringsAsFactors = FALSE))
gex_map <- data.frame(Sample = gex_map$sample, Tissue = gex_map$Tissue, bc = gex_map$cell_barcode)
gex_map <- gex_map[!is.na(gex_map$Tissue), ]
capture_of <- function(tissue, pid) if (tissue=="LIVER") paste0("GC-WL-",pid,"-LIVER") else paste0("GC-WL-",pid,"-SVP")
tcr_file_for <- function(cap) { c1 <- file.path(tcr_dir, paste0(cap,"-TCR"),"outs","filtered_contig_annotations.csv")
  c2 <- file.path(tcr_dir, paste0(cap,"-TCR"),"filtered_contig_annotations.csv"); if (file.exists(c1)) c1 else c2 }
read_contig <- function(tissue, pid) { d <- read.csv(tcr_file_for(capture_of(tissue,pid)), stringsAsFactors=FALSE)
  if (tissue=="LIVER") return(d); keep <- gex_map$bc[gex_map$Sample==capture_of(tissue,pid) & gex_map$Tissue==tissue]; d[d$barcode %in% keep,] }
build <- function(tissue, pids) { ct <- combineTCR(lapply(pids, function(p) read_contig(tissue,p)), samples=pids,
                   removeNA=FALSE, removeMulti=FALSE, filterMulti=FALSE); df <- bind_rows(ct); df$Tissue<-tissue; df$Patient_ID<-df$sample; df }
tis_lev <- c("LIVER","PBMC","SAT","VAT")
all_tcells <- bind_rows(lapply(tis_lev, function(t) build(t, matched_pids))) %>% filter(!is.na(CTaa), CTaa!="")

# ---- per patient: top-N clonotypes x tissue proportions ---------------------
build_pt <- function(d) {
  tot <- sort(table(d$CTaa), decreasing = TRUE)
  top <- names(tot)[seq_len(min(N_TOP, length(tot)))]
  tis_tot <- table(factor(d$Tissue, levels = tis_lev))
  cc <- as.data.frame(table(CTaa = factor(d$CTaa, levels = top),
                            Tissue = factor(d$Tissue, levels = tis_lev)), stringsAsFactors = FALSE)
  cc <- cc[!is.na(cc$CTaa), ]
  cc$proportion <- cc$Freq / as.numeric(tis_tot[cc$Tissue])
  cc$proportion[!is.finite(cc$proportion)] <- 0
  cc$CTaa   <- factor(cc$CTaa, levels = top)             # rank 1 (largest) first
  cc$Tissue <- factor(cc$Tissue, levels = tis_lev)
  cc
}
per_pt <- lapply(matched_pids, function(pid) {
  d <- all_tcells[all_tcells$Patient_ID == pid, c("CTaa","Tissue")]
  if (!nrow(d)) return(NULL)
  x <- build_pt(d); x$Patient_ID <- pid; x
})
names(per_pt) <- matched_pids
write.csv(do.call(rbind, per_pt),
          file.path(tab_dir, "EDFig7_topClonotypes_alluvial_v7.csv"), row.names = FALSE)

# ---- one alluvial panel per patient (each with its own clonotype legend) ----
panel <- function(pid) {
  cc <- per_pt[[pid]]; if (is.null(cc)) return(NULL)
  top <- levels(cc$CTaa)
  pal <- setNames(viridisLite::magma(length(top), begin = 0.05, end = 0.95, direction = -1), top)
  ggplot(cc, aes(x = Tissue, y = proportion, stratum = CTaa, alluvium = CTaa, fill = CTaa)) +
    geom_alluvium(alpha = 0.75, decreasing = FALSE, linewidth = 0.05, colour = NA) +
    geom_stratum(width = 0.42, linewidth = 0.08, decreasing = FALSE, colour = "grey20") +
    scale_fill_manual(values = pal, name = "Clonotypes") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(title = paste0(pid, " | ", stage_map[pid]), x = NULL, y = "Proportion") +
    theme_nature(9) +
    theme(legend.text = element_text(size = 6), legend.title = element_text(size = 8),
          legend.key.size = unit(3, "mm"), legend.position = "right",
          axis.text = element_text(size = 8), axis.title = element_text(size = 8),
          plot.title = element_text(size = 10, hjust = 0, face = "bold"))
}
plots <- Filter(Negate(is.null), lapply(matched_pids, panel))
fig <- wrap_plots(plots, ncol = 4)
# larger canvas so the bigger fonts and per-panel clonotype legends stay legible
save_nature(fig, file.path(fig_dir, "EDFig7_topClonotypes_alluvial_v7"), 52, 64)

cat("Done Extended Data Figure 7 (alluvial, top", N_TOP, "clonotypes, 20 matched patients).\n")
# End of the script
