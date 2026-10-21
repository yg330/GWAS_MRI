# ===========================================================================
# resolve_GarciaMarin_positions.R
# Manuscript item: Novel-loci comparison
#
# What this does:
#   Resolves positions for the Garcia-Marin et al. SNP list.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Resolve hg19 positions for Garcia_Marin_NG sheet (Unique_SNPs.xlsx)
# Garcia_Marin_NG has rsIDs only, no chromosome/position at all. This
# script looks up the hg19 (GRCh37) chr/pos for every rsID directly
# from dbSNP via the offline Bioconductor SNPlocs package, and writes
# the result as new CHR_hg19 / POS_hg19 columns, WITHOUT touching any
# other sheet in the workbook.
#
# Requires: openxlsx, SNPlocs.Hsapiens.dbSNP144.GRCh37 (Bioconductor)
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
#   install.packages("openxlsx")
#
# Note: this is a ~200MB one-off download the first time it's used.
# If you'd rather not install it, the alternative is an online lookup
# via biomaRt against the GRCh37 Ensembl archive (grch37.ensembl.org)
# -- let me know if you'd prefer that version instead.

source("config.R")

suppressPackageStartupMessages({
  library(openxlsx)
})

# ---- EDIT THIS if the file isn't in your working directory ----
xlsx_in  <- "Unique_SNPs.xlsx"
xlsx_out <- "Unique_SNPs_updated.xlsx"   # writes to a NEW file -- doesn't touch your original
# If you already produced Unique_SNPs_updated.xlsx via the Sun_NComms liftover
# script, point xlsx_in at THAT file so both fixes end up in the same workbook:
if (file.exists("Unique_SNPs_updated.xlsx")) xlsx_in <- "Unique_SNPs_updated.xlsx"

sheet_name <- "Garcia_Marin_NG"

# ---- 1. Load the offline dbSNP position reference for GRCh37/hg19 ----
if (!requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP144.GRCh37))
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37

# ---- 2. Read the sheet ----
wb <- loadWorkbook(xlsx_in)
df <- read.xlsx(xlsx_in, sheet = sheet_name, colNames = TRUE)
names(df)[1] <- "SNP"
df$SNP <- as.character(df$SNP)
cat(sprintf("Read %d rsIDs from %s\n", nrow(df), sheet_name))

# ---- 3. Look up hg19 chr/pos for each rsID ----
res <- tryCatch(
  snpsById(snps, df$SNP, ifnotfound = "drop"),
  error = function(e) { cat("snpsById error:", conditionMessage(e), "\n"); NULL }
)

if (is.null(res) || length(res) == 0) {
  stop("No rsIDs resolved -- check the rsID format and SNPlocs package version.")
}

hits <- as.data.frame(res)

# Diagnostic -- if matching still fails, paste this output back so we can see
# exactly what format RefSNP_id is coming back in on your installed package version.
cat("Preview of resolved hits (first 3 rows):\n")
print(utils::head(hits, 3))
cat("Preview of input rsIDs:", paste(utils::head(df$SNP, 3), collapse = ", "), "\n")

# Normalize BOTH sides to the bare numeric rsID (strip any "rs" prefix) before
# matching. Different SNPlocs package versions store RefSNP_id with or without
# the "rs" prefix already -- assuming one convention (as the previous version of
# this script did) silently fails to match ANYTHING if the version installed
# uses the other convention, which is what happened here.
hits$rs_bare <- gsub("^rs", "", as.character(hits$RefSNP_id))
df$rs_bare   <- gsub("^rs", "", as.character(df$SNP))

df$CHR_hg19 <- gsub("^chr|^ch", "", as.character(hits$seqnames[match(df$rs_bare, hits$rs_bare)]))
df$POS_hg19 <- hits$pos[match(df$rs_bare, hits$rs_bare)]

n_found <- sum(!is.na(df$POS_hg19))
cat(sprintf(
  "Resolved hg19 position for %d / %d rsIDs (%d not found -- likely merged/withdrawn/renamed rs numbers)\n",
  n_found, nrow(df), nrow(df) - n_found
))
if (n_found < nrow(df)) {
  cat("Unresolved rsIDs:\n")
  print(df$SNP[is.na(df$POS_hg19)])
}

# ---- 4. Write back ----
writeData(wb, sheet = sheet_name, x = df, startRow = 1, colNames = TRUE)
saveWorkbook(wb, xlsx_out, overwrite = TRUE)

cat("\nDone. Wrote", xlsx_out, "with CHR_hg19 / POS_hg19 columns added to the", sheet_name, "sheet.\n")
cat("  All other sheets in the workbook are unchanged.\n")
