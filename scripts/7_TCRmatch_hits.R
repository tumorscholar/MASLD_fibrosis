###############################################################################
# 7 - Apply the gap-run TCRMatch (IEDB) hits to v7_all_clones_annotation.csv  (R)
#
# The annotation build already added the 129 gap clones (10291-1 PBMC) as rows
# (Annotated = No). TCRMatch runs 100 seqs at a time, so the 111 gap betas gave
# TWO result files. This flips ONLY those 129 gap CTaa whose CDR3-beta got an
# IEDB hit to BacterialViral = Yes / Annotated = Yes. In-place update: no new
# rows, existing Yes annotations and all liver/cross-sectional rows untouched.
# Re-run 36_Figure5_TCRdonut.R afterwards.
###############################################################################
ann_dir  <- "/data/Blizard-AlazawiLab/rk/Tcell/CloneAnnotation"
ann_csv  <- file.path(ann_dir, "v7_all_clones_annotation.csv")
gap_csv  <- file.path(ann_dir, "v7_missing_clones_to_annotate.csv")     # Tissue,Donor,CTaa,CDR3b
# both TCRMatch result files (gzip-compressed .tsv). Add more if you split >200 seqs.
tcrm_files <- file.path(ann_dir, c("tcrmatch_gap_10291.tsv", "tcrmatch_gap_10291_1.tsv"))

# ---- 1. TCRMatch hit set (trimmed_input_sequence; already C.../...F-trimmed) --
# files are gzip-compressed despite the .tsv name; gzfile() reads them (and reads
# a plain uncompressed .tsv too, so it is safe either way). Let read.delim open
# and close the connection - do NOT close it yourself (double-close = "invalid
# connection").
read_tcrmatch <- function(f) {
 tryCatch(read.delim(gzfile(f), stringsAsFactors = FALSE, check.names = FALSE),
          error = function(e) { message("could not read ", f, ": ", conditionMessage(e)); NULL })
}
hits <- unique(unlist(lapply(tcrm_files, function(f) {
 d <- read_tcrmatch(f)
 if (is.null(d) || !nrow(d)) return(character(0))
 col <- if ("trimmed_input_sequence" %in% names(d)) "trimmed_input_sequence" else names(d)[1]
 d[[col]]
})))
hits <- hits[!is.na(hits) & hits != ""]
cat("TCRMatch hit sequences:", length(hits), "\n")

# ---- 2. gap CTaa set (the 129 10291-1 PBMC clones) ---------------------------
gap <- read.csv(gap_csv, stringsAsFactors = FALSE, check.names = FALSE)
gap_ctaa <- unique(gap$CTaa)
cat("gap CTaa:", length(gap_ctaa), "\n")

# TCRMatch strips a leading C and trailing F from the CDR3-beta -> same key here
trim_key <- function(b) { b <- ifelse(is.na(b), "", b)
sub("F$", "", sub("^C", "", b)) }

# ---- 3. flip matched gap rows in place ---------------------------------------
ann <- read.csv(ann_csv, stringsAsFactors = FALSE, check.names = FALSE)
sel <- ann$CTaa %in% gap_ctaa & trim_key(ann$CTaaB) %in% hits & ann$Annotated != "Yes"
cat("gap rows flipped to Yes:", sum(sel), "\n")
print(ann$CTaa[sel])
ann$BacterialViral[sel] <- "Yes"
ann$Annotated[sel]      <- "Yes"

write.csv(ann, ann_csv, row.names = FALSE)   # quotes character cols (matches file)
cat(sprintf("Wrote %s | Annotated=Yes now %d / %d (%.1f%%)\n",
            ann_csv, sum(ann$Annotated == "Yes"), nrow(ann),
            100 * mean(ann$Annotated == "Yes")))
cat("Now re-run 36_Figure5_TCRdonut.R.\n")
# End of the script