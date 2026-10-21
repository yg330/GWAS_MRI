# ===========================================================================
# liftover_Sun_NComms.R
# Manuscript item: Novel-loci comparison
#
# What this does:
#   Lifts the SNP list from Sun et al. (Nat Commun) from hg38 to hg19.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Liftover Sun_NComms sheet (Unique_SNPs.xlsx) from hg38 -> hg19
# Sun_NComms is the only sheet in Unique_SNPs.xlsx reported on hg38
# (column "POS (hg38)") -- everything else, including our own ST3
# meta-analysis, is hg19. This script lifts SUN_NCOMMS's reported
# hg38 positions over to hg19 using the standard UCSC chain file, and
# writes the result as a new POS_hg19 column, WITHOUT touching any
# other sheet in the workbook.
#
# Requires: openxlsx, rtracklayer, GenomicRanges (Bioconductor)
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("rtracklayer", "GenomicRanges"))
#   install.packages("openxlsx")

source("config.R")

suppressPackageStartupMessages({
  library(openxlsx)
  library(rtracklayer)
  library(GenomicRanges)
})

# ---- EDIT THIS if the file isn't in your working directory ----
xlsx_in  <- "Unique_SNPs.xlsx"
xlsx_out <- "Unique_SNPs_updated.xlsx"   # writes to a NEW file -- doesn't touch your original
sheet_name <- "Sun_NComms"
header_row <- 3   # actual column header is on Excel row 3 (1-indexed) in this sheet

# ---- 1. Read the sheet (openxlsx keeps the rest of the workbook intact) ----
wb <- loadWorkbook(xlsx_in)
df <- read.xlsx(xlsx_in, sheet = sheet_name, startRow = header_row, colNames = TRUE)

pos_col <- grep("^POS", names(df), value = TRUE)[1]
cat("Detected hg38 position column:", pos_col, "\n")

df$CHR <- as.character(df$CHR)
df[[pos_col]] <- as.numeric(df[[pos_col]])

# ---- 2. Download the hg38 -> hg19 chain file (one-off, cached in cwd) ----
chain_path <- "hg38ToHg19.over.chain"
if (!file.exists(chain_path)) {
  chain_gz <- paste0(chain_path, ".gz")
  cat("Downloading UCSC hg38->hg19 chain file...\n")
  download.file(
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
    chain_gz, mode = "wb"
  )
  con_in  <- gzfile(chain_gz, "rt")
  con_out <- file(chain_path, "wt")
  writeLines(readLines(con_in), con_out)
  close(con_in); close(con_out)
}
chain <- import.chain(chain_path)

# ---- 3. Build GRanges from hg38 coordinates and liftOver ----
valid <- !is.na(df$CHR) & !is.na(df[[pos_col]]) & df$CHR != ""
gr <- GRanges(
  seqnames = paste0("chr", df$CHR[valid]),
  ranges   = IRanges(start = df[[pos_col]][valid], width = 1),
  idx      = which(valid)   # remember original row index
)

lifted <- liftOver(gr, chain)

# liftOver can return 0 or >1 mapped ranges per input -- only keep clean 1:1 mappings;
# ambiguous / unmapped variants are left as NA and should be checked manually if needed
n_hits <- lengths(lifted)
cat(sprintf(
  "Lifted %d / %d variants to a unique hg19 position; %d failed to map, %d mapped ambiguously (>1 region)\n",
  sum(n_hits == 1), length(gr), sum(n_hits == 0), sum(n_hits > 1)
))

df$POS_hg19 <- NA_real_
ok <- which(n_hits == 1)
lifted_ok <- unlist(lifted[ok])
df$POS_hg19[gr$idx[ok]] <- start(lifted_ok)

# ---- 4. Write the updated sheet into a copy of the workbook ----
writeData(wb, sheet = sheet_name, x = df, startRow = header_row, colNames = TRUE)
saveWorkbook(wb, xlsx_out, overwrite = TRUE)

cat("\nDone. Wrote", xlsx_out, "with a new POS_hg19 column added to the", sheet_name, "sheet.\n")
cat(sprintf("  %d / %d rows now have an hg19 position.\n", sum(!is.na(df$POS_hg19)), nrow(df)))
cat("  All other sheets in the workbook are unchanged.\n")
