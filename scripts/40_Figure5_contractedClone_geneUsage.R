###############################################################################
# 40 - Figure 5 (supplementary) - are the FIBROSIS-associated TCR receptors the
#      ones that CONTRACT after surgery?  (runs LIVER and PBMC)
#
# Rationale: in MASLD fibrosis the repertoire is skewed toward TRBV20-1, TRBV28
# and TRBJ2-1 (Extended Data Fig. 5a). If those disease-associated receptors mark
# the pathogenic expanded clones, they should be preferentially found among the
# clones that shrink or disappear after bariatric surgery.
#
# Design: baseline-expanded clones tracked to follow-up (baseline_expanded_clones_
# tracking.csv) are split by FATE -
#     Lost                 = not re-detected at follow-up
#     Contracted           = re-detected, lower repertoire frequency
#     Maintained/Expanded  = re-detected, frequency >= baseline
# For each clone the TCR-beta V and J genes are recovered from the v7 contigs
# (combineTCR CTgene) and gene usage is compared between fates. A per-gene Fisher
# test asks whether TRBV20-1 / TRBV28 / TRBJ2-1 are enriched in the clones that
# contract/are lost (Contracted+Lost) versus those Maintained/Expanded.
#
# Runs for LIVER and PBMC and writes a combined figure.
#
# Inputs : baseline_expanded_clones_tracking.csv, v7_tcr_TR contigs
# Outputs: per-tissue Fig5_contracted_receptor_enrichment_<TISSUE>_v7.csv +
#          V/J usage csvs, per-tissue plots, and a combined bothTissues figure.
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
suppressMessages({ library(scRepertoire); library(dplyr); library(tidyr); library(ggplot2); library(readr); library(patchwork) })
source("00_nature_style.R")

fib_genes <- c("TRBV20-1", "TRBV28", "TRBJ2-1")       # fibrosis-enriched (ED Fig 5a)
tcr_dir <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
ann_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/CloneAnnotation"
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure5"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure5"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# baseline patient -> follow-up capture id (same map as scripts 34/35)
fu_map <- c("10738" = "11596", "10291-1" = "11303", "11040" = "11816", "11183" = "11937")
caps_avail <- list.dirs(tcr_dir, recursive = FALSE, full.names = FALSE)

run_tissue <- function(TISSUE) {
  cat("\n=====================", TISSUE, "=====================\n")

  # ---- 1. clone fate from the tracking table -------------------------------
  trk <- readr::read_csv(file.path(ann_dir, "baseline_expanded_clones_tracking.csv"),
                         show_col_types = FALSE) %>%
    filter(Tissue == TISSUE) %>%
    mutate(Baseline_Percent = as.numeric(Baseline_Percent),
           Followup_Percent = as.numeric(Followup_Percent),
           fate = dplyr::case_when(
             Presence != "Present in Both"        ~ "Lost",
             Followup_Percent <  Baseline_Percent ~ "Contracted",
             TRUE                                 ~ "Maintained/Expanded"),
           contracted = fate %in% c("Lost", "Contracted"))
  trk$fate <- factor(trk$fate, levels = c("Lost", "Contracted", "Maintained/Expanded"))
  if (!nrow(trk)) { cat("  no clones tracked for", TISSUE, "- skipping.\n"); return(NULL) }
  cat("clones tracked:", nrow(trk), "\n"); print(table(trk$fate))

  # ---- 2. TCR-beta V/J gene per CTaa (from v7 contigs) ---------------------
  tis_rx <- if (TISSUE == "LIVER") "LIVER" else "PBMC|SVP"
  find_cap <- function(id) { hit <- grep(paste0("^GC-WL-", id, "-.*(", tis_rx, ")"), caps_avail, value = TRUE)
                             if (length(hit)) hit[1] else NA_character_ }
  ids  <- unique(c(names(fu_map), unname(fu_map)))
  caps <- na.omit(unname(vapply(ids, find_cap, character(1))))
  cat("TCR captures used:", paste(caps, collapse = ", "), "\n")
  cl   <- lapply(caps, function(c) read.csv(file.path(tcr_dir, c, "outs", "filtered_contig_annotations.csv"),
                                            stringsAsFactors = FALSE))
  comb <- combineTCR(cl, samples = caps, removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)
  allc <- bind_rows(comb)
  beta <- sub("^[^_]*_", "", allc$CTgene)
  Vb <- sub("\\..*", "", beta); Vb[!grepl("^TRBV", Vb)] <- NA
  Jb <- sub(".*?(TRBJ[0-9-]+).*", "\\1", beta); Jb[!grepl("TRBJ", Jb)] <- NA
  gmap <- data.frame(CTaa = allc$CTaa, Vbeta = Vb, Jbeta = Jb, stringsAsFactors = FALSE) %>%
    filter(!is.na(CTaa)) %>% distinct(CTaa, Vbeta, Jbeta) %>%
    group_by(CTaa) %>% dplyr::slice(1) %>% ungroup()

  dat <- trk %>% left_join(gmap, by = "CTaa")
  cat("clones with a recovered beta V gene:", sum(!is.na(dat$Vbeta)), "/", nrow(dat), "\n")

  # ---- 3. per-gene enrichment: Contracted+Lost vs Maintained/Expanded ------
  test_gene <- function(col, gene) {
    d <- dat[!is.na(dat[[col]]), ]
    a <- sum(d[[col]] == gene &  d$contracted); b <- sum(d[[col]] == gene & !d$contracted)
    cc<- sum(d[[col]] != gene &  d$contracted); dd<- sum(d[[col]] != gene & !d$contracted)
    ft <- fisher.test(matrix(c(a, cc, b, dd), 2))
    data.frame(Tissue = TISSUE, gene = gene, n_contracted_lost = a, n_maintained = b,
               pct_contracted_lost = round(100 * a / max(1, a + cc), 1),
               pct_maintained = round(100 * b / max(1, b + dd), 1),
               odds_ratio = round(unname(ft$estimate), 2), p = signif(ft$p.value, 3))
  }
  res <- bind_rows(lapply(fib_genes, function(g) test_gene(if (grepl("TRBJ", g)) "Jbeta" else "Vbeta", g)))
  cat("\n== fibrosis-associated receptors by clone fate (", TISSUE, ") ==\n"); print(res, row.names = FALSE)
  write.csv(res, file.path(tab_dir, sprintf("Fig5_contracted_receptor_enrichment_%s_v7.csv", TISSUE)), row.names = FALSE)

  usage <- function(col) dat %>% filter(!is.na(.data[[col]])) %>%
    count(fate, gene = .data[[col]]) %>% group_by(fate) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
  write.csv(usage("Vbeta"), file.path(tab_dir, sprintf("Fig5_contracted_Vgene_usage_%s_v7.csv", TISSUE)), row.names = FALSE)
  write.csv(usage("Jbeta"), file.path(tab_dir, sprintf("Fig5_contracted_Jgene_usage_%s_v7.csv", TISSUE)), row.names = FALSE)

  # ---- 4. per-tissue plot --------------------------------------------------
  pl <- res %>% select(gene, pct_contracted_lost, pct_maintained) %>%
    pivot_longer(-gene, names_to = "grp", values_to = "pct") %>%
    mutate(grp = recode(grp, pct_contracted_lost = "Contracted / Lost", pct_maintained = "Maintained / Expanded"),
           gene = factor(gene, levels = fib_genes))
  p <- ggplot(pl, aes(gene, pct, fill = grp)) +
    geom_col(position = position_dodge(0.7), width = 0.62) +
    geom_text(data = res, inherit.aes = FALSE,
              aes(x = gene, y = pmax(pct_contracted_lost, pct_maintained) + 3,
                  label = paste0("p=", signif(p, 2))), size = 2.2) +
    scale_fill_manual(values = c("Contracted / Lost" = "#2166ac", "Maintained / Expanded" = "#b2182b"), name = NULL) +
    labs(x = NULL, y = "% of clones using the gene",
         title = sprintf("%s: fibrosis-associated TCR genes by clone fate", TISSUE)) +
    theme_nature(7) + theme(legend.position = "bottom")
  # per-tissue panel not saved on its own - only the combined LIVER|PBMC figure below is used.

  for (i in seq_len(nrow(res))) with(res[i, ],
    cat(sprintf("  %s: %.1f%% contracted/lost vs %.1f%% maintained (OR=%.2f, p=%s) -> %s\n",
        gene, pct_contracted_lost, pct_maintained, odds_ratio, format(p, scientific = TRUE),
        ifelse(p < 0.05 & odds_ratio > 1, "ENRICHED in contracted",
        ifelse(p < 0.05 & odds_ratio < 1, "enriched in PERSISTING", "no significant difference")))))
  list(res = res, plot = p)
}

out <- lapply(c("LIVER", "PBMC"), run_tissue)
out <- out[!vapply(out, is.null, logical(1))]

# ---- combined figure (liver | pbmc) + combined table -----------------------
if (length(out) >= 2) {
  combined <- out[[1]]$plot | out[[2]]$plot
  save_nature(combined, file.path(fig_dir, "Fig5_contracted_receptor_enrichment_bothTissues_v7"), 18, 8)
}
allres <- bind_rows(lapply(out, `[[`, "res"))
write.csv(allres, file.path(tab_dir, "Fig5_contracted_receptor_enrichment_bothTissues_v7.csv"), row.names = FALSE)
cat("\nDone Fig 5 contracted-clone receptor usage (LIVER + PBMC).\n")
# End of the script
