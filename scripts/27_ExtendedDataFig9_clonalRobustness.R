###############################################################################
# 27 - Extended Data Figure 9 - clonal-expansion robustness (Fig 2 supplementary)
#
# Demonstrates that the clonal expansion is not driven by any single donor or any
# single clone. Uses the same v7_tcr_TR pooled repertoire as the Figure 2/3
# scripts. Three analyses:
#   A. Per-donor expanded fraction (expansion defined within each donor x tissue),
#      box plots by tissue x stage with per-donor points and per-donor Wilcoxon.
#   B. Leave-one-donor-out: each donor is dropped in turn and the pooled percentage
#      expanded (overall and per tissue) recomputed; the range across the 21 runs
#      quantifies the influence of any single donor.
#   C. Drop-largest-clone: the percentage expanded is recomputed after removing the
#      largest clone and the top five clones.
# Also reports per-donor normalised Shannon diversity by tissue x stage.
#
# Read-only on the data; writes figures and tables only.
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
suppressMessages({library(scRepertoire); library(tidyverse)})
source("00_nature_style.R")
if (!exists("stage_cols")) stage_cols <- c("No_fibrosis" = "#8E44AD", "Fibrosis" = "#C2185B", "Healthy" = "#7F8C8D")

tcr_dir     <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
barcode_dir <- "/data/Blizard-AlazawiLab/rk/barcodes_v7"
fig_dir     <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure2"
tab_dir     <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure2"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- pooled TCR repertoire (identical loader to 11 / 31) --------------------
liver_pids <- c("9680","9932","9961","9991","9999","10113-1","10113-2","10202","10203","10205",
                "10291-1","10291-2","10380","10634","10738","10742","11040","11051","11183","11327","11471")
svp_pids   <- setdiff(liver_pids, "10291-1")
stage_map  <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No_fibrosis","9991"="Fibrosis","9999"="No_fibrosis",
                "10113-1"="Fibrosis","10113-2"="No_fibrosis","10202"="Fibrosis","10203"="Fibrosis","10205"="Fibrosis",
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
tissues <- list(LIVER=liver_pids, PBMC=svp_pids, SAT=svp_pids, VAT=svp_pids)
all_tcells <- bind_rows(Map(build, names(tissues), tissues)) %>% filter(!is.na(CTaa), CTaa!="")
all_tcells$Stage <- unname(stage_map[all_tcells$Patient_ID])

tis_lev   <- c("LIVER","PBMC","SAT","VAT")
stage_lev <- c("No_fibrosis","Fibrosis")               # stats: MASLD only (Healthy n=1)
all_tcells$Tissue <- factor(all_tcells$Tissue, levels = tis_lev)
star <- function(p) ifelse(is.na(p),"", ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*","ns"))))

# pooled expansion (manuscript definition: pooled clone size > 1)
cc      <- table(all_tcells$CTaa)
exp_all <- names(cc)[cc > 1]
n_tot   <- nrow(all_tcells)
cat(sprintf("Pooled repertoire: %d T cells, %d clonotypes, %.1f%% expanded (>1 cell).\n",
            n_tot, length(cc), 100*sum(cc[cc>1])/n_tot))

# =====================================================================
# A. PER-DONOR expanded fraction (expansion defined WITHIN donor x tissue)
# =====================================================================
MINCELL <- 20                                          # ignore tiny donor x tissue samples
# clone size WITHIN each donor x tissue via add_count (robust; avoids ave()-in-summarise)
perdonor <- all_tcells %>%
  add_count(Patient_ID, Tissue, Stage, CTaa, name = "clone_n") %>%
  group_by(Patient_ID, Tissue, Stage) %>%
  filter(dplyr::n() >= MINCELL) %>%
  summarise(n_cells = dplyr::n(),
            expanded_frac = 100 * sum(clone_n > 1) / dplyr::n(), .groups = "drop")
write.csv(perdonor, file.path(tab_dir, "Fig2_perdonor_expandedFraction_v7.csv"), row.names = FALSE)

pd <- perdonor %>% filter(Stage %in% stage_lev) %>%
  mutate(Stage = factor(Stage, levels = stage_lev))
statA <- pd %>% group_by(Tissue) %>%
  summarise(p = tryCatch(wilcox.test(expanded_frac ~ Stage)$p.value, error=function(e) NA_real_),
            n_Nofib = sum(Stage=="No_fibrosis"), n_Fib = sum(Stage=="Fibrosis"),
            y = max(expanded_frac)*1.05, .groups="drop") %>% mutate(lab = star(p))
write.csv(statA, file.path(tab_dir, "Fig2_perdonor_expandedFraction_stats_v7.csv"), row.names = FALSE)
cat("\nA. Per-donor expanded fraction, No-fib vs Fibrosis (per-donor Wilcoxon):\n")
print(as.data.frame(statA), row.names = FALSE)

pA <- ggplot(pd, aes(Stage, expanded_frac, fill = Stage)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.8) +
  geom_text(data = statA, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE, size = 2.2) +
  facet_wrap(~ Tissue, nrow = 1) +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_x_discrete(labels = c("No_fibrosis" = "No fibrosis", "Fibrosis" = "Fibrosis")) +
  labs(x = NULL, y = "% of donor's T cells in expanded clones",
       title = "Per-donor clonal expansion (expansion defined within donor; each point = one donor)") +
  theme_nature(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_nature(pA, file.path(fig_dir, "Fig2_perdonor_expandedFraction_v7"), 18, 7)

# =====================================================================
# B. LEAVE-ONE-DONOR-OUT on the pooled headline % expanded
# =====================================================================
donors <- sort(unique(all_tcells$Patient_ID))
overall_full <- 100 * sum(cc[cc>1]) / n_tot
# fraction of a set's cells in clones expanded WITHIN each tissue.
# Pure base R (ave + tapply) so a masked dplyr verb (scRepertoire loads plyr) can
# never break it. Returns a data.frame(Tissue, pct).
within_tis <- function(df) {
  csz <- ave(rep(1L, nrow(df)), df$Tissue, df$CTaa, FUN = length)   # clone size within tissue
  agg <- tapply(csz > 1, as.character(df$Tissue), function(v) 100 * mean(v))
  data.frame(Tissue = names(agg), pct = as.numeric(agg), stringsAsFactors = FALSE)
}
tis_full <- within_tis(all_tcells); names(tis_full)[names(tis_full) == "pct"] <- "full"

loo <- lapply(donors, function(d) {
  sub <- all_tcells %>% filter(Patient_ID != d)
  ccx <- table(sub$CTaa)
  data.frame(dropped = d,
             overall = 100 * sum(ccx[ccx>1]) / nrow(sub),
             within_tis(sub) %>% tidyr::pivot_wider(names_from = Tissue, values_from = pct))
})
loo <- bind_rows(loo)
write.csv(loo, file.path(tab_dir, "Fig2_leaveOneDonorOut_expanded_v7.csv"), row.names = FALSE)
loo_rng <- loo %>% summarise(across(where(is.numeric), list(min = min, max = max)))
cat(sprintf("\nB. Leave-one-donor-out: overall %% expanded = %.1f%% (full); LOO range %.1f-%.1f%% over %d donors.\n",
            overall_full, min(loo$overall), max(loo$overall), length(donors)))
cat("   Per-tissue full vs LOO range:\n"); print(tis_full)

# range plot: full value (point) with LOO min-max whiskers, per tissue + overall (base R)
num <- loo[, setdiff(names(loo), "dropped"), drop = FALSE]
loo_long <- data.frame(grp = names(num),
                       lo = sapply(num, min), hi = sapply(num, max),
                       stringsAsFactors = FALSE, row.names = NULL)
full_vals <- rbind(data.frame(grp = "overall", full = overall_full, stringsAsFactors = FALSE),
                   data.frame(grp = as.character(tis_full$Tissue), full = tis_full$full, stringsAsFactors = FALSE))
loo_long <- merge(loo_long, full_vals, by = "grp", all.x = TRUE)
loo_long$grp <- factor(loo_long$grp, levels = c("overall", tis_lev))
pB <- ggplot(loo_long, aes(grp, full)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.4, colour = "grey40") +
  geom_point(size = 1.6, colour = "#C2185B") +
  labs(x = NULL, y = "% T cells in expanded clones",
       title = "Leave-one-donor-out stability (point = full data; whisker = range dropping any one donor)") +
  theme_nature(7)
save_nature(pB, file.path(fig_dir, "Fig2_leaveOneDonorOut_expanded_v7"), 12, 7)

# =====================================================================
# C. DROP-LARGEST-CLONE(S): headline stability to removing top clones
# =====================================================================
ord <- sort(cc, decreasing = TRUE)
drop_top <- function(k) {
  drop_ct <- names(ord)[seq_len(k)]
  sub <- all_tcells %>% filter(!CTaa %in% drop_ct)
  ccx <- table(sub$CTaa)
  data.frame(removed = paste0("top ", k), n_cells_removed = sum(cc[drop_ct]),
             overall_pct = 100 * sum(ccx[ccx>1]) / nrow(sub),
             liver_pct = { lv <- sub[as.character(sub$Tissue)=="LIVER", ]; 100*mean(ave(rep(1L,nrow(lv)), lv$CTaa, FUN=length) > 1) })
}
cloneC <- bind_rows(data.frame(removed="none", n_cells_removed=0,
                               overall_pct = overall_full,
                               liver_pct = tis_full$full[tis_full$Tissue=="LIVER"]),
                    drop_top(1), drop_top(5))
write.csv(cloneC, file.path(tab_dir, "Fig2_dropTopClones_expanded_v7.csv"), row.names = FALSE)
cat("\nC. Drop-largest-clone(s) - headline % expanded barely moves:\n")
print(cloneC, row.names = FALSE)

# ---- panel C figure: PER-DONOR stability, faceted BY TISSUE -----------------
# (ED Fig 6c.) The pooled table above collapses everyone together; here each
# donor is shown separately, within each tissue, so it is visible whether the
# expansion is broad across donors or driven by a few - and how each donor's
# expanded fraction moves when that donor's OWN largest clone / top-5 clones are
# removed. Each line = one donor; colour = disease stage.
drop_perdonor <- function() {
  keys <- all_tcells %>% dplyr::count(Patient_ID, Tissue, Stage, name = "n") %>%
    dplyr::filter(n >= MINCELL)
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    sub  <- all_tcells[all_tcells$Patient_ID == keys$Patient_ID[i] &
                         as.character(all_tcells$Tissue) == as.character(keys$Tissue[i]), ]
    ctab <- sort(table(sub$CTaa), decreasing = TRUE)
    do.call(rbind, lapply(c(0, 1, 5), function(k) {
      drop <- if (k > 0) names(ctab)[seq_len(min(k, length(ctab)))] else character(0)
      s2   <- sub[!sub$CTaa %in% drop, ]
      cs   <- ave(rep(1L, nrow(s2)), s2$CTaa, FUN = length)
      data.frame(Patient_ID = keys$Patient_ID[i], Tissue = as.character(keys$Tissue[i]),
                 Stage = keys$Stage[i],
                 removed = c("none", "top 1", "top 5")[match(k, c(0, 1, 5))],
                 expanded_frac = 100 * mean(cs > 1), stringsAsFactors = FALSE)
    }))
  }))
}
ppC <- drop_perdonor()
ppC$removed <- factor(ppC$removed, levels = c("none", "top 1", "top 5"))
ppC$Tissue  <- factor(ppC$Tissue,  levels = tis_lev)
# order the Stage legend Healthy -> No fibrosis -> Fibrosis (matches other figures)
ppC$Stage   <- factor(ppC$Stage,   levels = c("Healthy", "No_fibrosis", "Fibrosis"))
write.csv(ppC, file.path(tab_dir, "Fig2_dropTopClones_perDonor_v7.csv"), row.names = FALSE)

pC <- ggplot(ppC, aes(removed, expanded_frac, group = Patient_ID)) +
  geom_line(alpha = 0.5, linewidth = 0.3, colour = "grey60") +
  geom_point(aes(colour = Stage), size = 0.9) +
  facet_wrap(~ Tissue, nrow = 1) +
  scale_colour_manual(values = stage_cols, name = "Stage",
                      breaks = c("Healthy", "No_fibrosis", "Fibrosis"),
                      labels = c("Healthy", "No fibrosis", "Fibrosis")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Clones removed (per donor)", y = "% of donor's T cells in expanded clones",
       title = "Per-donor stability of clonal expansion to removing that donor's largest clone(s)") +
  theme_nature(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right", legend.key.size = unit(0.7, "lines"))
save_nature(pC, file.path(fig_dir, "Fig2_dropTopClones_perDonor_v7"), 18, 7)

# =====================================================================
# Bonus: per-donor normalised Shannon diversity by tissue x stage
# =====================================================================
# normalised Shannon per donor x tissue x stage - base R (split/reduce), masking-proof
SEP <- "|||"
gkey <- paste(all_tcells$Patient_ID, all_tcells$Tissue, all_tcells$Stage, sep = SEP)
div_rows <- lapply(split(all_tcells$CTaa, gkey), function(ctaa) {
  nc <- as.numeric(table(ctaa)); p <- nc / sum(nc); K <- length(nc)
  c(n_cells = sum(nc), nShannon = if (K <= 1) NA_real_ else -sum(p * log(p)) / log(K))
})
div <- data.frame(do.call(rbind, div_rows))
gp  <- do.call(rbind, strsplit(names(div_rows), SEP, fixed = TRUE))
div$Patient_ID <- gp[,1]; div$Tissue <- gp[,2]; div$Stage <- gp[,3]
div <- div[div$n_cells >= MINCELL & div$Stage %in% stage_lev, ]
div$Stage  <- factor(div$Stage,  levels = stage_lev)
div$Tissue <- factor(div$Tissue, levels = tis_lev)
statD <- do.call(rbind, lapply(split(div, droplevels(div$Tissue)), function(d) {
  if (!nrow(d)) return(NULL)
  p <- tryCatch(wilcox.test(nShannon ~ Stage, data = d)$p.value, error = function(e) NA_real_)
  data.frame(Tissue = d$Tissue[1], p = p, y = max(d$nShannon, na.rm = TRUE) * 1.02, lab = star(p))
}))
write.csv(div,   file.path(tab_dir, "Fig2_perdonor_diversity_v7.csv"), row.names = FALSE)
pD <- ggplot(div, aes(Stage, nShannon, fill = Stage)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.8) +
  geom_text(data = statD, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE, size = 2.2) +
  facet_wrap(~ Tissue, nrow = 1) +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_x_discrete(labels = c("No_fibrosis" = "No fibrosis", "Fibrosis" = "Fibrosis")) +
  labs(x = NULL, y = "Normalised Shannon diversity (per donor)",
       title = "Per-donor TCR diversity by stage (lower = more clonal)") +
  theme_nature(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_nature(pD, file.path(fig_dir, "Fig2_perdonor_diversity_v7"), 18, 7)

cat("\nDone Fig 2 robustness. Stars are per-donor Wilcoxon (No-fib vs Fibrosis); with 3 No-fib vs 17",
    "\nFibrosis donors the per-donor tests are underpowered - the robustness message is the spread (A),",
    "\nthe LOO stability (B) and the clone-drop stability (C), not the p-values.\n")
# End of the script