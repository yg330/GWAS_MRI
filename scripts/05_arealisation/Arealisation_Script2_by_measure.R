# ===========================================================================
# Arealisation_Script2_by_measure.R
# Manuscript item: Arealisation validation, step 2
#
# What this does:
#   Reads each regional GWAS file once per measure and pulls out all index
#   SNPs, resumable by measure. Meant for the cluster; set SUMSTAT_DIR to
#   the folder of regional fastGWA files.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Arealisation validation -- Script 2, BY-MEASURE, RESUMABLE
# Single self-contained script for interactive use in RStudio.
# For each of the ~12 measures (SA, ICVF, MD, ISOVF, GMV, CT, OD,
# IC, FI, GC, FA, MC), reads every regional parcel file for that
# measure exactly ONCE. Each file is genome-wide, so a single awk
# pass per file can pull out ALL of that measure's index SNPs at
# once (0 to all of them), rather than re-opening the same file once
# per SNP.
#
# RESUMABLE AT THE MEASURE LEVEL: if
# Arealisation_betas/by_measure/{measure}.csv already exists, that
# measure is skipped entirely and the script moves on. To resume
# after an interruption, just re-source this script -- completed
# measures are detected and skipped automatically. Once every
# measure's file exists, it also writes the combined
# index_snp_betas_by_parcel.csv (Script 3's input).
#
# INPUT (place in the working directory)
#   index_snps_per_gene_measure.csv   (Script 1's output)
# CLUSTER PATH
#   .../yg330/meta_MRI_GWAS/UKB_all_GWAS/regional/{MEASURE}/  -- one
#   subfolder per measure (CT, FA, FI, GC, GMV, IC, ICVF, ISOVF, MC,
#   MD, OD, SA), each containing plain-text (NOT gzipped) native
#   GCTA .fastGWA files, ~1.32GB each, one per parcel. This is the
#   REAL raw fastGWA output -- the earlier GWAS_sumstat/*.sumstats.gz
#   folder turned out to be LDSC-munged (SNP/A1/A2/N/Z only, no
#   CHR/POS/BETA, HapMap3-restricted) and can't be used for this.

source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(parallel)
  library(purrr)
  library(data.table)
})

ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

sumstat_dir    <- Sys.getenv("SUMSTAT_DIR")
index_snp_path <- "index_snps_per_gene_measure.csv"
out_dir        <- "Arealisation_betas"
by_measure_dir <- file.path(out_dir, "by_measure")
dir.create(by_measure_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PART 1: manifest -- parse MEASURE/PARCEL out of every filename.
# Files live under sumstat_dir/{MEASURE}/combined_GWAS_{TOOL}_regional_
# {MEASURE}_{PARCEL}_avg_scaled_as_pheno{N}.fastGWA -- scan
# recursively across all measure subfolders in one pass. {TOOL} is
# NOT constant: structural measures (CT/FI/GC/GMV/IC/MC/SA) use
# "STR", diffusion measures (FA/ICVF/ISOVF/MD/OD) use "DTI" --
# the regex below accepts any alphanumeric token there rather than
# hardcoding one, so a third convention wouldn't silently fail again.
#
# CACHED: the recursive scan over ~4,700 files is real filesystem
# metadata traffic on a networked mount, not free -- re-running it
# every time you re-source the script just to resume is wasted work.
# Cached to manifest_cache.csv after the first successful scan; every
# later re-source loads that instead. Set RESCAN <- TRUE (or just
# delete manifest_cache.csv) if files on disk have actually changed.
# ============================================================
RESCAN <- FALSE
manifest_cache_path <- file.path(out_dir, "manifest_cache.csv")

if (!RESCAN && file.exists(manifest_cache_path)) {
  ts("Loading cached manifest (set RESCAN <- TRUE to force a rescan)...")
  manifest <- read_csv(manifest_cache_path, show_col_types = FALSE)
  cat(sprintf("  %d regional files loaded from cache across %d measures\n",
              nrow(manifest), n_distinct(manifest$measure)))
} else {
  ts("Scanning sumstat directory (recursive across measure subfolders)...")
  all_paths <- list.files(sumstat_dir, pattern = "\\.fastGWA$", recursive = TRUE, full.names = TRUE)
  fname_re  <- "^combined_GWAS_[A-Za-z]+_regional_([A-Za-z]+)_(.+)_avg_scaled_as_pheno(\\d+)\\.fastGWA$"

  manifest <- tibble(path = all_paths, file = basename(all_paths)) %>%
    mutate(measure = str_match(file, fname_re)[, 2],
           parcel  = str_match(file, fname_re)[, 3]) %>%
    dplyr::filter(!is.na(measure))
  n_unparsed <- length(all_paths) - nrow(manifest)
  if (n_unparsed > 0) {
    cat(sprintf("  WARNING: %d .fastGWA file(s) didn't match the expected naming pattern -- inspect manually\n", n_unparsed))
  }
  cat(sprintf("  %d regional files parsed across %d measures\n", nrow(manifest), n_distinct(manifest$measure)))
  print(table(manifest$measure))
  write_csv(manifest, manifest_cache_path)
  cat(sprintf("  Cached manifest -> %s\n", manifest_cache_path))
}

# ============================================================
# PART 2: index SNPs from Script 1
# ============================================================
ts("Loading index SNPs from Script 1...")
index_snps <- read_csv(index_snp_path, show_col_types = FALSE) %>%
  dplyr::filter(status == "ok") %>%
  dplyr::distinct(symbol, phenotype_type, locus_rank, index_snp, coord_id, chr, pos) %>%
  dplyr::mutate(chr = as.character(chr))
cat(sprintf("  %d resolved index SNPs across %d measures\n",
            nrow(index_snps), n_distinct(index_snps$phenotype_type)))
write_csv(index_snps, file.path(out_dir, "extraction_tasks.csv"))

missing_from_manifest <- setdiff(unique(index_snps$phenotype_type), manifest$measure)
if (length(missing_from_manifest) > 0) {
  cat("  WARNING: no manifest files at all for these measures:\n")
  print(missing_from_manifest)
}

# Order measures by ASCENDING file count (not by target SNP count --
# runtime is driven by how many parcel files must be read, which is
# independent of how many target SNPs that measure has). This means
# whichever measure has the fewest parcel files runs first, giving
# you a quick correctness check and an early per-file timing signal
# before committing to the slower/larger measures.
files_per_measure <- manifest %>%
  dplyr::filter(measure %in% unique(index_snps$phenotype_type)) %>%
  dplyr::count(measure, name = "n_files") %>%
  dplyr::arrange(n_files)
cat("\n  Processing order (smallest file count first):\n")
print(files_per_measure)

measures_to_do <- files_per_measure$measure
# Any measure entirely missing from the manifest (0 files) won't
# appear in files_per_measure at all -- append it at the end so it's
# still handled (and flagged) rather than silently skipped.
measures_to_do <- c(measures_to_do, setdiff(missing_from_manifest, measures_to_do))

# ============================================================
# PART 3: per-file extraction helper -- ONE awk pass per file,
# matching against every target (chr,pos) for that file's measure
# at once. Target set is passed inline via -v (small: <=61 entries
# currently), no scratch files needed.
#
# These .fastGWA files are PLAIN TEXT (not gzipped, unlike the old
# GWAS_sumstat/*.sumstats.gz folder) -- awk reads the path directly,
# no zcat needed. Column layout is native GCTA fastGWA, tab-
# delimited, 1 header row: CHR SNP POS A1 A2 N AF1 BETA SE P.
# NOTE: this column layout has NOT been independently verified for
# THIS specific file set yet (the earlier "confirmed plain fastGWA
# columns" turned out to describe the wrong folder entirely) -- if
# the first measure processed below comes back with an unexpectedly
# low hit rate again, check `head -3` of one real file here before
# assuming it's a data-availability issue.
# ============================================================
# N_CORES: how many files to read CONCURRENTLY. Do NOT trust
# parallel::detectCores() on a shared HPC node -- it often reports
# the whole physical node's core count, not what your interactive
# RStudio session was actually allocated, and oversubscribing will
# make things slower, not faster (context-switching + I/O
# contention). Set this to whatever core count you requested when
# launching the RStudio interactive session.
#
# Set to 3 (matches a 3-core RStudio allocation). Note: requesting
# multiple NODES does not add to this -- mclapply forks within a
# single machine only, so an interactive RStudio session is capped
# at whatever cores it has on the ONE node it's running on, regardless
# of how many nodes were requested for the session. The parent R
# process is essentially idle while its forked children wait on the
# awk subprocesses, so using all available cores for children is
# fine here -- expect roughly ~3x over the fully serial version.
N_CORES <- 3

extract_one_file <- function(path, target_keys) {
  # LC_ALL=C disables locale-aware string handling in awk, which is
  # unnecessary here (plain ASCII numeric/rsID fields) and measurably
  # slows down byte-for-byte comparison/hashing in some locales.
  cmd <- sprintf(
    "LC_ALL=C awk -F'\\t' -v targets=%s 'BEGIN{n=split(targets,arr,\",\"); for(i=1;i<=n;i++) want[arr[i]]=1} NR>1{key=$1\"_\"$3; if(key in want) print}' %s",
    shQuote(paste(target_keys, collapse = ",")), shQuote(path))
  res <- tryCatch(fread(cmd = cmd, header = FALSE, sep = "\t"), error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(tibble())
  tibble(chr = as.character(res[[1]]), rs_snp = as.character(res[[2]]),
         pos = as.numeric(res[[3]]), a1 = as.character(res[[4]]), a2 = as.character(res[[5]]),
         n = as.numeric(res[[6]]), af1 = as.numeric(res[[7]]),
         beta = as.numeric(res[[8]]), se = as.numeric(res[[9]]), p = as.numeric(res[[10]]))
}

process_measure <- function(m) {
  out_path <- file.path(by_measure_dir, sprintf("%s.csv", m))
  targets <- index_snps %>% dplyr::filter(phenotype_type == m) %>% dplyr::distinct(chr, pos)
  target_keys <- paste(targets$chr, targets$pos, sep = "_")

  files_m <- manifest %>% dplyr::filter(measure == m)
  cat(sprintf("  %s: %d target SNPs, %d parcel files to scan\n", m, nrow(targets), nrow(files_m)))

  if (nrow(files_m) == 0) {
    write_csv(tibble(), out_path)
    return(invisible(NULL))
  }

  # Read files CONCURRENTLY (mc.cores) -- each file's awk pass is an
  # independent OS subprocess, so this is safe: order doesn't matter,
  # results are collected and combined after all finish. mclapply
  # uses fork() (Linux/Mac only, fine on the cluster). If any single
  # file's extraction errors, mclapply returns a try-error object for
  # that element rather than killing the whole batch -- filtered out
  # below rather than crashing the measure.
  ts(sprintf("  %s: reading %d files across %d cores...", m, nrow(files_m), N_CORES))
  hit_list <- parallel::mclapply(seq_len(nrow(files_m)), function(i) {
    h <- extract_one_file(files_m$path[i], target_keys)
    if (nrow(h) > 0) h$parcel <- files_m$parcel[i]
    h
  }, mc.cores = N_CORES)
  hit_list <- hit_list[!vapply(hit_list, function(x) inherits(x, "try-error"), logical(1))]
  hits <- dplyr::bind_rows(hit_list)

  # Expected grid for this measure: every index SNP x every parcel
  # file that exists for it -- left-joined to whatever was actually
  # found, so genuinely-missing (chr,pos,parcel) combos are recorded
  # as "not_found_in_file" rather than silently absent.
  expected <- index_snps %>%
    dplyr::filter(phenotype_type == m) %>%
    tidyr::crossing(parcel = files_m$parcel)

  if (nrow(hits) == 0) {
    merged <- expected %>% dplyr::mutate(status = "not_found_in_file",
                                         rs_snp = NA_character_, a1 = NA_character_, a2 = NA_character_,
                                         n = NA_real_, af1 = NA_real_, beta = NA_real_, se = NA_real_, p = NA_real_)
  } else {
    merged <- expected %>%
      dplyr::left_join(hits, by = c("chr", "pos", "parcel")) %>%
      dplyr::mutate(status = ifelse(is.na(beta), "not_found_in_file", "ok"))
  }

  write_csv(merged, out_path)
  cat(sprintf("  %s: wrote %d rows (%d ok, %d not_found) -> %s\n",
              m, nrow(merged), sum(merged$status == "ok"), sum(merged$status == "not_found_in_file"), out_path))
}

# ============================================================
# PART 4: loop over measures, skipping any already done
# ============================================================
for (m in measures_to_do) {
  out_path <- file.path(by_measure_dir, sprintf("%s.csv", m))
  if (file.exists(out_path)) {
    cat(sprintf("SKIP %s -- already done (%s)\n", m, out_path))
    next
  }
  ts(sprintf("=== Processing measure: %s ===", m))
  process_measure(m)
}

# ============================================================
# PART 5: if every measure is done, write the combined output
# ============================================================
done_measures <- measures_to_do[file.exists(file.path(by_measure_dir, sprintf("%s.csv", measures_to_do)))]
if (setequal(done_measures, measures_to_do)) {
  ts("All measures done -- writing combined output...")
  combined <- purrr::map_dfr(measures_to_do, ~ read_csv(file.path(by_measure_dir, sprintf("%s.csv", .x)),
                                                         show_col_types = FALSE))
  write_csv(combined, file.path(out_dir, "index_snp_betas_by_parcel.csv"))
  cat(sprintf("ALL DONE. Wrote %d rows -> %s  <- Script 3 input\n",
              nrow(combined), file.path(out_dir, "index_snp_betas_by_parcel.csv")))
} else {
  cat(sprintf("\n%d / %d measures done. Re-source this script to continue with: %s\n",
              length(done_measures), length(measures_to_do),
              paste(setdiff(measures_to_do, done_measures), collapse = ", ")))
}
ts("=== script finished (this pass) ===")
