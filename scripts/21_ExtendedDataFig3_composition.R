###############################################################################
# 21 - Extended Data Figure 3 - per-donor cell-type composition (Fig 1 supplementary)
#
# For every cell type, in every tissue, the proportion of that cell type is
# computed per donor, plotted as box plots by disease stage, and tested Fibrosis
# vs No Fibrosis (Wilcoxon, BH-adjusted; Healthy n = 1 shown for reference only).
# Provides the per-donor, per-cell-type composition comparison across tissues and
# disease stages that accompanies the Figure 1 atlas.
###############################################################################
suppressMessages({library(Seurat); library(dplyr); library(tidyr); library(ggplot2)})
source("00_nature_style.R")

seu <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds")
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure1"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure1"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

md <- seu@meta.data
ct <- intersect(c("cell_type","celltype"), colnames(md))[1]

# Seurat meta.data can carry list-columns / factors, and 'count' may be masked by
# another loaded package (e.g. plyr) - both throw "Argument 'x' is not a vector:
# list". Pull the needed columns into a CLEAN atomic data.frame first. as.vector
# unwraps any 1-column list, then unlist forces atomic; then coerce to character.
flat <- function(x) { x <- if (is.list(x)) vapply(x, function(v) paste(unlist(v), collapse=";"), "") else x
as.character(x) }
cat("meta column classes:  Patient_ID:", class(md$Patient_ID),
    " Tissue:", class(md$Tissue), " ", ct, ":", class(md[[ct]]),
    " Stage:", class(md$Stage), "\n")
df <- data.frame(
 Patient_ID = flat(md$Patient_ID),
 Tissue     = flat(md$Tissue),
 cell_type  = flat(md[[ct]]),
 Stage      = gsub("No_Fibrosis", "No_fibrosis", flat(md$Stage)),
 stringsAsFactors = FALSE)
stage_lev <- c("Healthy","No_fibrosis","Fibrosis")

# ---- per-donor proportion of each cell type within each tissue --------------
comp <- df %>%
 dplyr::count(Patient_ID, Tissue, cell_type, name = "n") %>%
 group_by(Patient_ID, Tissue) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
 tidyr::complete(tidyr::nesting(Patient_ID, Tissue), cell_type, fill = list(n = 0, pct = 0))
# attach stage (one per patient)
pstage <- df %>% distinct(Patient_ID, Stage)
comp <- comp %>% left_join(pstage, by = "Patient_ID") %>%
 mutate(Stage = factor(Stage, levels = stage_lev))

# ---- statistics: Fibrosis vs No Fibrosis per (Tissue, cell_type) ------------
stat <- comp %>% filter(Stage %in% c("No_fibrosis","Fibrosis")) %>%
 group_by(Tissue, cell_type) %>%
 summarise(p = tryCatch(wilcox.test(pct ~ Stage)$p.value, error = function(e) NA_real_),
           med_Nofib = median(pct[Stage=="No_fibrosis"]),
           med_Fib   = median(pct[Stage=="Fibrosis"]), .groups = "drop") %>%
 group_by(Tissue) %>% mutate(padj = p.adjust(p, "BH")) %>% ungroup()
write.csv(stat, file.path(tab_dir, "Fig1_composition_stats_FibVsNoFib_v7.csv"), row.names = FALSE)
cat("Significant (padj<0.05) cell-type shifts Fibrosis vs No Fibrosis:\n")
print(as.data.frame(stat %>% filter(padj < 0.05) %>% arrange(Tissue, padj)), row.names = FALSE)

star <- function(p) ifelse(is.na(p),"", ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*",""))))

# ---- one box-plot figure per tissue (facet by cell type) --------------------
for (tis in c("LIVER","PBMC","SAT","VAT")) {
 d <- comp %>% filter(Tissue == tis)
 s <- stat %>% filter(Tissue == tis) %>%
  mutate(lab = star(padj), y = tapply(d$pct, d$cell_type, max)[cell_type] * 1.05)
 p <- ggplot(d, aes(Stage, pct, fill = Stage)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.2) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.6) +
  geom_text(data = s, aes(x = 2.5, y = y, label = lab), inherit.aes = FALSE, size = 2) +
  facet_wrap(~ cell_type, scales = "free_y") +
  scale_fill_manual(values = stage_cols, name = "Stage") +
  labs(x = NULL, y = paste0("% of ", tis, " cells (per donor)"),
       title = paste0(tis, " - cell-type composition by stage (Fib vs No-Fib: BH-adj Wilcoxon)")) +
  theme_nature(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
        strip.text = element_text(size = 5.5))
 save_nature(p, file.path(fig_dir, paste0("Fig1_composition_box_", tolower(tis), "_v7")), 22, 20)
}
cat("\nDone Fig 1 composition box plots (per-donor, with stats). Healthy shown for reference (n=1).\n")
# End of the script