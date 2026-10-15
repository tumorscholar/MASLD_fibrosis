###############################################################################
# 05b - VDJ-based T-cell % per tissue x STAGE (the numbers for the Results text)
#
# Identical computation to 05_Figure2a.R (a cell is a T cell if its barcode has a
# productive TRA/TRB contig; denominator = QC-passed cells), but summarised per
# Tissue x Stage so the manuscript can quote the same VDJ definition shown in
# Fig 2a. Reports BOTH statistics:
#   mean_pct  = mean of per-patient percentages  <- this is what the Fig 2a bar shows
#   pooled_pct= 100 * sum(T cells) / sum(cells)  <- cohort pooled
# Use mean_pct in the text so it matches the figure. Read-only.
###############################################################################
library(Seurat); library(tidyverse)

gex_rds <- "/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds"
tcr_dir <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure2"
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
productive_only <- TRUE

tcr_file_for <- function(sample) {
 cands <- c(file.path(tcr_dir, paste0(sample, "-TCR"), "outs", "filtered_contig_annotations.csv"),
            file.path(tcr_dir, paste0(sample, "-TCR"),          "filtered_contig_annotations.csv"),
            file.path(tcr_dir, sample, "outs",                  "filtered_contig_annotations.csv"),
            file.path(tcr_dir, sample,                          "filtered_contig_annotations.csv"))
 hit <- cands[file.exists(cands)]; if (length(hit)) hit[1] else NA_character_
}
stage_map <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No_fibrosis","9991"="Fibrosis",
               "9999"="No_fibrosis","10113-1"="Fibrosis","10113-2"="No_fibrosis","10202"="Fibrosis",
               "10203"="Fibrosis","10205"="Fibrosis","10291-1"="Fibrosis","10291-2"="Fibrosis",
               "10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis","10742"="Healthy",
               "11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis")

obj <- readRDS(gex_rds)
gex <- obj@meta.data %>% mutate(bc = sub("^.*_", "", rownames(.))) %>%
 transmute(Sample, Patient_ID, Tissue, bc)

captures <- unique(gex$Sample)
tcr_bc <- lapply(captures, function(s) {
 f <- tcr_file_for(s); if (is.na(f)) { message("No TCR file for capture: ", s); return(NULL) }
 d <- read.csv(f, stringsAsFactors = FALSE)
 if (productive_only && "productive" %in% names(d)) d <- d[tolower(as.character(d$productive)) %in% c("true","1"), , drop = FALSE]
 if ("chain" %in% names(d)) d <- d[d$chain %in% c("TRA","TRB"), , drop = FALSE]
 if (!nrow(d)) return(NULL)
 data.frame(Sample = s, bc = unique(d$barcode), stringsAsFactors = FALSE)
}) %>% bind_rows()

gex$hasTCR <- paste(gex$Sample, gex$bc) %in% paste(tcr_bc$Sample, tcr_bc$bc)

frac <- gex %>% group_by(Patient_ID, Tissue) %>%
 summarise(nTotal = n(), nT = sum(hasTCR), .groups = "drop") %>%
 mutate(pct = 100 * nT / nTotal,
        Stage  = factor(stage_map[Patient_ID], levels = c("Healthy","No_fibrosis","Fibrosis")),
        Tissue = factor(Tissue, levels = c("LIVER","PBMC","SAT","VAT")))

perstage <- frac %>% group_by(Tissue, Stage) %>%
 summarise(n_patients = n(),
           mean_pct   = round(mean(pct), 1),                       # <- Fig 2a bar value
           pooled_pct = round(100 * sum(nT) / sum(nTotal), 1),
           .groups = "drop") %>% arrange(Tissue, Stage)
cat("\n=== VDJ T-cell % per tissue x stage (mean_pct = Fig 2a bar; use this in the text) ===\n")
print(as.data.frame(perstage), row.names = FALSE)
write.csv(perstage, file.path(tab_dir, "Fig2a_Tcell_pct_perTissueStage_V7.csv"), row.names = FALSE)

cat("\n>>> For para 62 (LIVER), quote mean_pct: ",
    "Fibrosis / No_fibrosis / Healthy = ",
    paste(perstage$mean_pct[perstage$Tissue=="LIVER"][match(c("Fibrosis","No_fibrosis","Healthy"),
                                                            perstage$Stage[perstage$Tissue=="LIVER"])], collapse = " / "), "\n")
# End of the script