###############################################################################
# Figure 2d - antigen-annotation donut, ALL 4 TISSUES as concentric rings
#
# Same categorisation / colours / assign_cat as Fig2d_annotation_donut.R
# (TCRdb condition + VDJdb species + IEDB bacterial/viral), but split by tissue
# and drawn as 4 concentric rings in one plot:
#     inner -> outer :  LIVER, PBMC, SAT, VAT
# A clone shared across tissues appears in each tissue's own ring.
#
# Source: v7_all_clones_annotation.csv (per tissue; TCRdb_Condition, BacterialViral)
#         + VDJdb_species joined from v7_clones_annotation_combined.csv by CDR3.
###############################################################################
library(dplyr); library(ggplot2); library(readr)

ann_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/CloneAnnotation"
out_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- keyword -> category + assign_cat (identical to Fig 2d) -----------------
categorize <- function(cond) {
 x <- tolower(cond)
 if (grepl("cancer|tumou?r|carcinoma|melanoma|leuk[a]?emia|lymphoma|sarcoma|myeloma|glioblastoma|nsclc|neoplas|adenocarcinoma|metasta|malignan|adoptive t-cell therapy|hsct|\\baml\\b|lymphoblastic|hepatocellular|urothelial|osteosarcoma|liposarcoma|glioma|solid tumor", x, perl = TRUE)) return("Cancer")
 if (grepl("cmv|ebv|influenza|sars|cov|hiv|hcv|hbv|hepatitis|dengue|denv|htlv|vzv|yfv|yellow fever|mcmv|lcmv|mcpyv|tuberculosis|salmonella|covid|cytomegalovirus|mononucleosis|vaccin|streptomyces|villitis|infectious|virus|viral|bacter|paratyphi|typhi|leukoencephalopathy", x, perl = TRUE)) return("Microbial Antigen")
 if (grepl("autoimmun|arthritis|sclerosis|sjogren|graves|spondylitis|lupus|psoriasis|crohn|colitis|inflammatory bowel|sarcoidosis|t1d|type 1 diab|narcolepsy|susac|aplastic anemia|transplant|graft-versus|gvhd|aire|ipex|immunodeficiency|scid|cvid|granulomatous|fanconi|depleted|wiskott", x, perl = TRUE)) return("Autoimmune & Transplantation")
 if (grepl("t2d|type 2 diab|metabolic|masld|nafld|obesity|steato|nash", x, perl = TRUE)) return("Metabolic diseases")
 if (grepl("healthy|homosapiens|self", x, perl = TRUE)) return("Healthy")
 "Unknown antigen"
}
priority <- c("Metabolic diseases", "Autoimmune & Transplantation", "Microbial Antigen",
              "Cancer", "Healthy", "Unknown antigen")
assign_cat <- function(tcrdb, vdj, bv) {
 conds <- trimws(unlist(strsplit(paste(tcrdb, vdj, sep = ";"), "[;,]")))
 conds <- conds[conds != "" & !conds %in% c("#N/A", "NA", "Unknown")]
 cats  <- unique(vapply(conds, categorize, character(1)))
 if (identical(as.character(bv), "Yes")) cats <- union(cats, "Microbial Antigen")
 if (length(cats) == 0) return("Unannotated")
 hit <- priority[priority %in% cats]; if (length(hit)) hit[1] else "Unknown antigen"
}

# ---- load per-tissue clones + join VDJdb species ----------------------------
ann <- readr::read_csv(file.path(ann_dir, "v7_all_clones_annotation.csv"), show_col_types = FALSE)
cmb <- readr::read_csv(file.path(ann_dir, "v7_clones_annotation_combined.csv"), show_col_types = FALSE)
vdj <- cmb %>% select(CTaaA, CTaaB, VDJdb_species) %>% distinct()

d <- ann %>% left_join(vdj, by = c("CTaaA", "CTaaB")) %>%
 mutate(TCRdb_Condition = ifelse(is.na(TCRdb_Condition), "", TCRdb_Condition),
        VDJdb_species   = ifelse(is.na(VDJdb_species),   "", VDJdb_species),
        BacterialViral  = ifelse(is.na(BacterialViral),  "", BacterialViral),
        Category = mapply(assign_cat, TCRdb_Condition, VDJdb_species, BacterialViral))
write.csv(d[, c("Tissue","CTaaA","CTaaB","TCRdb_Condition","VDJdb_species","BacterialViral","Category")],
          file.path(out_dir, "Fig2d_clone_categories_4tissue_V7.csv"), row.names = FALSE)

# ---- per-tissue fractions ---------------------------------------------------
cat_order <- c("Cancer","Healthy","Microbial Antigen","Autoimmune & Transplantation",
               "Metabolic diseases","Unknown antigen","Unannotated")
tissues   <- c("LIVER", "PBMC", "SAT", "VAT")                    # inner -> outer
comp <- d %>% filter(Tissue %in% tissues) %>%
 dplyr::count(Tissue, Category) %>% group_by(Tissue) %>% mutate(frac = n / sum(n)) %>% ungroup() %>%
 mutate(Category = factor(Category, levels = cat_order)) %>% arrange(Tissue, Category) %>%
 group_by(Tissue) %>% mutate(ymax = cumsum(frac), ymin = dplyr::lag(ymax, default = 0)) %>% ungroup()
ring <- setNames(seq_along(tissues), tissues)
comp$xmin <- ring[comp$Tissue]; comp$xmax <- comp$xmin + 0.9
for (tt in tissues) { s <- d %>% filter(Tissue == tt)
cat(sprintf("%s: n=%d, matched = %.1f%%\n", tt, nrow(s), 100 * mean(s$Category != "Unannotated"))) }
write.csv(comp, file.path(out_dir, "Fig2d_donut_counts_4tissue_V7.csv"), row.names = FALSE)

# ---- 4-ring donut -----------------------------------------------------------
cat_cols <- c("Cancer"="#3C3C3C", "Healthy"="#007E7E", "Microbial Antigen"="#FCD860",
              "Autoimmune & Transplantation"="#CCE4D2", "Metabolic diseases"="#F6EAC6",
              "Unknown antigen"="#424242", "Unannotated"="#1E72A2")
mtab <- d %>% filter(Tissue %in% tissues) %>% group_by(Tissue) %>%
 summarise(m = round(100 * mean(Category != "Unannotated")), .groups = "drop")
lab <- data.frame(Tissue = tissues, x = ring[tissues] + 0.45, y = 0,
                  label = paste0(tissues, " (", mtab$m[match(tissues, mtab$Tissue)], "%)"))

p <- ggplot(comp) +
 geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Category),
           colour = "white", linewidth = 0.3) +
 coord_polar(theta = "y") + xlim(0, max(comp$xmax) + 0.15) +
 scale_fill_manual(values = cat_cols, name = NULL, breaks = cat_order) +
 geom_text(data = lab, aes(x = x, y = y, label = label), size = 2.0, vjust = -0.3) +
 theme_void(base_size = 7) +
 theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"),
       legend.position = "right", plot.title = element_text(size = 8, hjust = 0.5)) +
 ggtitle("Expanded clone antigen annotation by tissue")

grDevices::cairo_pdf(file.path(out_dir, "Fig2d_annotation_donut_4tissue_V7.pdf"), width = 12/2.54, height = 8/2.54)
print(p); dev.off()
tiff(file.path(out_dir, "Fig2d_annotation_donut_4tissue_V7.tiff"), width = 12, height = 8, units = "cm", res = 600, compression = "lzw")
print(p); dev.off()
message("Saved Fig2d_annotation_donut_4tissue_V7 (pdf + tiff) + category/donut CSVs")
# End of the script
