###############################################################################
# 29 - Extended Table 1 - public hyperexpanded TCR clonotypes shared across patients
#        patients.
#
# Uses the SAME source as the Fig 2/3 figures and the manuscript counts:
# the v7_tcr_TR captures pooled by combineTCR, with the SVP library demultiplexed
# to SAT/VAT/PBMC via the barcode map (identical loader to 11_Figure3a-c.R /
# 31_recompute_all_numbers_v7.R). This keeps the table consistent with the "two
# public clones" reported in the text (a clonotype is public if detected in >=2
# patients; hyperexpanded if pooled clone size >20 cells).
#
# Columns: Ctaa-a_Ctaa-b, Epitope Annotation, n_patients, Tissues, Stage, n_cells
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
library(scRepertoire); library(tidyverse)

tcr_dir     <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
barcode_dir <- "/data/Blizard-AlazawiLab/rk/barcodes_v7"
ann_dir     <- "/data/Blizard-AlazawiLab/rk/Tcell/CloneAnnotation"
tab_dir     <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/ExtendedData"
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---- cohort + stage (identical to 11 / 31) ----------------------------------
liver_pids <- c("9680","9932","9961","9991","9999","10113-1","10113-2","10202","10203","10205",
                "10291-1","10291-2","10380","10634","10738","10742","11040","11051","11183","11327","11471")
svp_pids   <- setdiff(liver_pids, "10291-1")
stage_map  <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No Fibrosis","9991"="Fibrosis","9999"="No Fibrosis",
                "10113-1"="Fibrosis","10113-2"="No Fibrosis","10202"="Fibrosis","10203"="Fibrosis","10205"="Fibrosis",
                "10291-1"="Fibrosis","10291-2"="Fibrosis","10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis",
                "10742"="Healthy","11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis")

# ---- barcode -> tissue map (pre-QC demux) -----------------------------------
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
tissues <- list(LIVER=liver_pids, PBMC=svp_pids, SAT=svp_pids, VAT=svp_pids)
all_tcells <- bind_rows(Map(build, names(tissues), tissues)) %>% filter(!is.na(CTaa), CTaa!="")
all_tcells$Stage <- stage_map[all_tcells$Patient_ID]

# ---- per-clonotype aggregation ----------------------------------------------
cc  <- table(all_tcells$CTaa)
hyper <- names(cc)[cc > 20]
cat("Hyperexpanded clonotypes (>20 cells):", length(hyper), "\n")
meta <- all_tcells %>% filter(CTaa %in% hyper) %>% group_by(CTaa) %>%
  summarise(n_cells = n(), n_patients = n_distinct(Patient_ID),
            Tissues = paste(sort(unique(Tissue)), collapse = ","),
            Stage   = paste(sort(unique(Stage)), collapse = ", "), .groups = "drop")

# ---- epitope annotation (TCRdb condition), joined by CTaa -------------------
ann <- readr::read_csv(file.path(ann_dir, "v7_all_clones_annotation.csv"), show_col_types = FALSE) %>%
  mutate(TCRdb_Condition = ifelse(is.na(TCRdb_Condition), "", TCRdb_Condition)) %>%
  group_by(CTaa) %>%
  summarise(`Epitope Annotation` = { v <- unique(TCRdb_Condition[TCRdb_Condition != "" & TCRdb_Condition != "#N/A"])
  if (length(v)) paste(v, collapse = ",") else "#N/A" }, .groups = "drop")

tbl <- meta %>% filter(n_patients >= 2) %>%          # public
  left_join(ann, by = "CTaa") %>%
  mutate(`Epitope Annotation` = ifelse(is.na(`Epitope Annotation`), "#N/A", `Epitope Annotation`)) %>%
  arrange(desc(n_cells)) %>%
  transmute(`Ctaa-a_Ctaa-b` = CTaa, `Epitope Annotation`, n_patients, Tissues, Stage, n_cells)

cat("Public hyperexpanded (>=2 patients):", nrow(tbl), "\n")
print(as.data.frame(tbl), row.names = FALSE)
write.csv(tbl, file.path(tab_dir, "ExtendedTable1_public_hyperexpanded_clones_v7.csv"), row.names = FALSE)
message("Saved ExtendedTable1_public_hyperexpanded_clones_v7.csv (v7_tcr_TR source; matches 'two public clones' in text)")
# End of the script