# ===========================================================================
# emVAR_FLAMES_geneprioritisationcomparison.R
# Manuscript item: Results: emVar / ABC concordance
#
# What this does:
#   Compares FLAMES-prioritised genes with ABC-linked genes at loci that
#   contain emVars, using a locus-stratified permutation test (10,000
#   permutations).
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# emVar (ABC-prioritised) x FLAMES gene concordance test
# PHENOTYPE-AWARE VERSION
#
# QUESTION: Within loci where a GWAS signal and an emVar/MPRA
# signal overlap, do ABC and FLAMES converge on the SAME gene
# FOR THE SAME PHENOTYPE more than expected by chance?
#
# WHY PHENOTYPE-MATCHING IS NECESSARY (not just position):
# This project runs many GWAS (one per imaging phenotype). A locus
# is only meaningfully "the same locus" for ABC and FLAMES if
# they're both being credited for the SAME trait -- a FLAMES hit
# for cortical thickness sitting near an ABC hit that was actually
# linked (via the overlap file) to an orientation-dispersion GWAS
# signal is a coincidence of geography, not evidence of
# concordance. Positional overlap alone (the previous version of
# this script) does not enforce this and would silently inflate
# the overlap count with phenotype-mismatched pairs.
#
# THREE-WAY MATCHING:
#   1. LOCUS -> phenotype: each lead_snp in the overlap file carries
#      a fixed gwas_phenotype (confirmed constant across its
#      repeated rows -- multiple emVars can share one lead_snp, but
#      the phenotype attached to that lead_snp doesn't change).
#   2. ABC -> phenotype: the ABC file itself has NO phenotype column.
#      The only way to attach one is to join each ABC row back to
#      the overlap file by EXACT variant identity (chr:pos + alleles,
#      checked in both possible orders since ref/alt vs A1/A2 could
#      be flipped), inheriting whichever gwas_phenotype that SPECIFIC
#      emVar was matched to -- not just "falls somewhere in the
#      merged window."
#   3. FLAMES -> phenotype: FLAMES already carries its own phenotype
#      per row. A FLAMES gene only counts as a hit at a locus if its
#      phenotype is actually present among that locus's
#      gwas_phenotype value(s).
#
# "Overlap" = same gene, same locus, AND matching phenotype on both
# sides -- not just positional co-occurrence.
#
# WINDOW MERGING: loci are still built as +-WINDOW_BP windows around
# each lead_snp, merged wherever they overlap (e.g. the MAPT/
# 17q21.31 region, where dozens of lead SNPs sit within ~800kb of
# each other) -- this remains purely geometric and is used only to
# define the candidate gene BACKGROUND (what could plausibly have
# been prioritised near these loci at all), not to launder phenotype
# matching across genuinely different GWAS signals that happen to
# be physically close.

source("config.R")

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(purrr)

ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
set.seed(42)

# ---- Paths ----
base_dir  <- proj_dir
out_dir   <- file.path(base_dir, "Outputs", "emVar_FLAMES_concordance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

overlap_path  <- file.path(base_dir, "Data", "emvar_gwas_overlap_r2_0.8.xlsx")  # sheet 1
flames_path   <- file.path(proj_dir, "Data", "FLAMES_gene_prioritisation_results.txt")
gene_loc_path <- file.path(base_dir, "Data", "NCBI37.3.gene.loc")

# Confirmed by the user -- .tsv, not the earlier .xlsx guess.
# abc_path <- "~/Downloads/NPC_ABC_clean_variant_gene_locus_z_mapping_ABC_proxy_z_gt1.5.tsv"
abc_path <- file.path(proj_dir, "Data", "NPC_ABC_all_tested_GWAS_LD_filtered_clean_variant_gene_locus_z_mapping_ABC_proxy_z_gt1.5.tsv")


WINDOW_BP <- 500000
N_PERM    <- 10000

# ---- FLAMES phenotype -> short code, matching the convention used ----
# throughout the rest of this project (e.g. "regional_SA_47l" -> "SA")
extract_phenotype_type <- function(phenotype_col) {
  cleaned <- phenotype_col %>% str_remove("^global_") %>% str_remove("^regional_")
  case_when(
    str_detect(cleaned, "^tract_") ~ str_extract(cleaned, "(?<=^tract_)[A-Za-z0-9]+"),
    TRUE ~ str_split(cleaned, "_") %>% map_chr(1)
  )
}

# ================================================================
# SECTION 1: Load gene positions (hg19, genome-wide)
# ================================================================
ts("Loading gene position table...")
gene_loc <- read.table(gene_loc_path, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("entrez_id", "chr", "start", "end", "strand", "symbol")) %>%
  dplyr::filter(chr %in% c(1:22, "X", "Y")) %>%
  mutate(chr = as.character(chr)) %>%
  distinct(symbol, .keep_all = TRUE)
cat(sprintf("  %d genes with hg19 positions\n", nrow(gene_loc)))

# ================================================================
# SECTION 2: Load the emVar-GWAS overlap file, extract per-lead_snp
#            phenotype, build loci
# ================================================================
ts("Loading emVar-GWAS overlap file...")
overlap_raw <- read_excel(overlap_path, sheet = 1)
cat(sprintf("  %d overlap rows (direct + LD_proxy)\n", nrow(overlap_raw)))

# One row per lead_snp with its (fixed) phenotype list, split into a
# character vector per lead_snp for later matching.
lead_snp_phenotypes <- overlap_raw %>%
  distinct(lead_snp, gwas_phenotype) %>%
  distinct(lead_snp, .keep_all = TRUE) %>%   # gwas_phenotype is constant per lead_snp; guard anyway
  mutate(phenotype_list = str_split(gwas_phenotype, ",\\s*")) %>%
  mutate(phenotype_types = map(phenotype_list, extract_phenotype_type))

lead_loci <- lead_snp_phenotypes %>%
  mutate(chr = str_remove(str_extract(lead_snp, "^[^:]+"), "^chr"),
         pos = as.integer(str_extract(lead_snp, "(?<=:)[0-9]+")),
         start = pmax(1L, pos - WINDOW_BP),
         end   = pos + WINDOW_BP) %>%
  dplyr::filter(!is.na(chr), !is.na(pos))
cat(sprintf("  %d unique lead SNPs -> %d unmerged locus windows\n",
            n_distinct(overlap_raw$lead_snp), nrow(lead_loci)))

# ---- Merge overlapping windows (geometry only -- background universe) ----
merge_intervals <- function(df) {
  df <- df %>% arrange(start)
  merged <- list()
  cur_start <- df$start[1]; cur_end <- df$end[1]
  cur_lead_snps <- df$lead_snp[1]
  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      if (df$start[i] <= cur_end) {
        cur_end <- max(cur_end, df$end[i])
        cur_lead_snps <- paste(cur_lead_snps, df$lead_snp[i], sep = ";")
      } else {
        merged[[length(merged) + 1]] <- tibble(start = cur_start, end = cur_end,
                                               lead_snps = cur_lead_snps)
        cur_start <- df$start[i]; cur_end <- df$end[i]; cur_lead_snps <- df$lead_snp[i]
      }
    }
  }
  merged[[length(merged) + 1]] <- tibble(start = cur_start, end = cur_end,
                                         lead_snps = cur_lead_snps)
  bind_rows(merged)
}

merged_loci <- lead_loci %>%
  group_by(chr) %>%
  group_modify(~ merge_intervals(.x)) %>%
  ungroup() %>%
  mutate(locus_id = paste0("locus_", row_number()),
         n_lead_snps = str_count(lead_snps, ";") + 1)

cat(sprintf("  %d unmerged windows -> %d merged loci\n", nrow(lead_loci), nrow(merged_loci)))
if (any(merged_loci$n_lead_snps > 5)) {
  cat("  NOTE: some merged loci absorbed >5 lead SNPs (likely MAPT/17q21.31 or similar):\n")
  print(merged_loci %>% dplyr::filter(n_lead_snps > 5) %>% dplyr::select(chr, start, end, n_lead_snps))
}

write.csv(merged_loci %>% dplyr::select(-lead_snps),
          file.path(out_dir, "merged_loci.csv"), row.names = FALSE)

# Map: which lead_snp(s) (and hence which phenotypes) fall inside
# each merged locus -- needed downstream for phenotype matching.
locus_lead_snp_map <- merged_loci %>%
  dplyr::select(locus_id, lead_snps) %>%
  mutate(lead_snp = str_split(lead_snps, ";")) %>%
  dplyr::select(-lead_snps) %>%
  unnest(lead_snp) %>%
  left_join(lead_snp_phenotypes %>% dplyr::select(lead_snp, phenotype_types), by = "lead_snp")

# ================================================================
# SECTION 3: Candidate gene universe -- positional, per locus
#            (background only -- NOT phenotype-restricted)
# ================================================================
ts("Building candidate gene universe (positional, per merged locus)...")

assign_genes_to_loci <- function(genes_df, loci_df) {
  map_dfr(seq_len(nrow(loci_df)), function(i) {
    loc <- loci_df[i, ]
    hits <- genes_df %>% dplyr::filter(chr == loc$chr, start <= loc$end, end >= loc$start)
    if (nrow(hits) == 0) return(NULL)
    tibble(locus_id = loc$locus_id, symbol = hits$symbol)
  })
}

candidate_universe <- assign_genes_to_loci(gene_loc, merged_loci)
cat(sprintf("  %d unique candidate genes across %d loci (background for the Fisher test)\n",
            n_distinct(candidate_universe$symbol), n_distinct(candidate_universe$locus_id)))

write.csv(candidate_universe, file.path(out_dir, "candidate_universe_by_locus.csv"), row.names = FALSE)

# ================================================================
# SECTION 4: ABC genes, phenotype inherited via EXACT variant match
#            back to the overlap file (not positional containment)
# ================================================================
ts("Loading ABC gene-variant mapping and joining to exact emVar identity...")
abc_raw <- read.table(abc_path, header = TRUE)
cat(sprintf("  %d ABC variant-gene rows (pre-filtered to z > 1.5)\n", nrow(abc_raw)))

# Exact-variant lookup from the overlap file, checked in BOTH allele
# orders since ref/alt (ABC) vs A1/A2 (overlap file) may be flipped.
emvar_lookup <- overlap_raw %>%
  distinct(emvar_chr_hg38, emvar_pos_hg38, emvar_A1, emvar_A2, lead_snp) %>%
  left_join(lead_snp_phenotypes %>% dplyr::select(lead_snp, phenotype_types), by = "lead_snp")

match_key_fwd <- emvar_lookup %>%
  transmute(key = paste(emvar_chr_hg38, emvar_pos_hg38, emvar_A1, emvar_A2, sep = "_"),
            lead_snp, phenotype_types)
match_key_rev <- emvar_lookup %>%
  transmute(key = paste(emvar_chr_hg38, emvar_pos_hg38, emvar_A2, emvar_A1, sep = "_"),
            lead_snp, phenotype_types)
emvar_key_table <- bind_rows(match_key_fwd, match_key_rev) %>% distinct(key, .keep_all = TRUE)

abc_matched <- abc_raw %>%
  mutate(variant_chr = str_remove(as.character(variant_chr), "^chr"),
         # ABC's variant_pos is BED-style 0-based; emvar_pos_hg38 (and
         # the rest of this pipeline) is 1-based -- confirmed via
         # multiple examples where ABC's position was consistently
         # exactly 1bp less than the matching emVar's hg38 position.
         variant_pos_1based = variant_pos + 1L,
         key = paste(variant_chr, variant_pos_1based, variant_ref, variant_alt, sep = "_")) %>%
  left_join(emvar_key_table, by = "key")

n_abc_matched <- sum(!is.na(abc_matched$lead_snp))
cat(sprintf("  %d / %d ABC rows matched an exact emVar in the overlap file\n",
            n_abc_matched, nrow(abc_raw)))
if (n_abc_matched < nrow(abc_raw)) {
  cat(sprintf("  NOTE: %d ABC row(s) did NOT match any emVar exactly -- these are dropped\n",
              nrow(abc_raw) - n_abc_matched))
  cat("  from the phenotype-aware analysis (no phenotype can be attached to them).\n")
  cat("  If this number is large, check for a build mismatch or allele-encoding difference.\n")
}

abc_by_locus_phenotype <- abc_matched %>%
  dplyr::filter(!is.na(lead_snp)) %>%
  left_join(locus_lead_snp_map %>% dplyr::select(locus_id, lead_snp), by = "lead_snp",
            relationship = "many-to-many") %>%
  dplyr::filter(!is.na(locus_id)) %>%
  dplyr::select(locus_id, lead_snp, symbol = gene_symbol, phenotype_types) %>%
  # Same universe restriction FLAMES already gets in Section 5 --
  # without this, an ABC gene outside the +-500kb window (ABC's
  # activity-by-contact scores aren't strictly distance-capped)
  # would count toward the Fisher table without being part of the
  # background, silently breaking the a+b+c+d = universe arithmetic.
  semi_join(candidate_universe, by = c("locus_id", "symbol")) %>%
  distinct()

n_abc_outside_universe <- abc_matched %>%
  dplyr::filter(!is.na(lead_snp)) %>%
  left_join(locus_lead_snp_map %>% dplyr::select(locus_id, lead_snp), by = "lead_snp",
            relationship = "many-to-many") %>%
  dplyr::filter(!is.na(locus_id)) %>%
  distinct(locus_id, symbol = gene_symbol) %>%
  anti_join(candidate_universe, by = c("locus_id", "symbol")) %>%
  nrow()
if (n_abc_outside_universe > 0) {
  cat(sprintf("  NOTE: %d ABC (gene, locus) pair(s) fell OUTSIDE the +-%dkb candidate\n",
              n_abc_outside_universe, WINDOW_BP / 1000))
  cat("  universe and were dropped from the Fisher test background for consistency\n")
  cat("  with how FLAMES genes are treated. If this number is large, the window may\n")
  cat("  be too narrow for how far ABC's activity-by-contact links can reach.\n")
}

cat(sprintf("  %d unique ABC (gene, locus, phenotype-linked) hits\n", nrow(abc_by_locus_phenotype)))
write.csv(abc_by_locus_phenotype %>% mutate(phenotype_types = map_chr(phenotype_types, paste, collapse = ";")),
          file.path(out_dir, "abc_gene_locus_phenotype.csv"), row.names = FALSE)

# ================================================================
# SECTION 5: FLAMES genes, phenotype-matched to each locus
# ================================================================
ts("Loading FLAMES genes and phenotype-matching to loci...")
flames_raw <- read.table(flames_path, header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(phenotype_type = extract_phenotype_type(phenotype))
cat(sprintf("  FLAMES genome-wide: %d unique genes across %d phenotype rows\n",
            n_distinct(flames_raw$symbol), nrow(flames_raw)))

# Position-restrict FLAMES genes to the candidate universe first
flames_in_universe <- flames_raw %>%
  dplyr::select(symbol, phenotype_type) %>%
  distinct() %>%
  inner_join(candidate_universe, by = "symbol", relationship = "many-to-many")

# Now require the FLAMES row's phenotype to actually be among the
# phenotype(s) attached to at least one lead_snp within that SAME
# locus -- this is the phenotype-matching step positional overlap
# alone does not enforce.
flames_by_locus_phenotype <- flames_in_universe %>%
  left_join(locus_lead_snp_map %>% dplyr::select(locus_id, phenotype_types), by = "locus_id",
            relationship = "many-to-many") %>%
  rowwise() %>%
  mutate(phenotype_match = phenotype_type %in% phenotype_types) %>%
  ungroup() %>%
  dplyr::filter(phenotype_match) %>%
  distinct(locus_id, symbol, phenotype_type)

cat(sprintf("  FLAMES hits after BOTH position AND phenotype restriction: %d unique genes\n",
            n_distinct(flames_by_locus_phenotype$symbol)))
cat("  (compare this to the naive genome-wide 848, and to the position-only-restricted\n")
cat("   count -- the gap between the two restricted numbers is exactly the effect of\n")
cat("   phenotype matching)\n")

write.csv(flames_by_locus_phenotype, file.path(out_dir, "flames_gene_locus_phenotype.csv"),
          row.names = FALSE)

# ================================================================
# SECTION 6: Fisher's exact test -- gene+locus+phenotype concordance
# ================================================================
ts("Running Fisher's exact test (phenotype-aware overlap definition)...")

# A true concordant hit now requires: same gene, same locus, AND at
# least one shared phenotype between the ABC row's phenotype_types
# and the FLAMES row's phenotype_type.
overlap_pairs <- abc_by_locus_phenotype %>%
  inner_join(flames_by_locus_phenotype, by = c("locus_id", "symbol"),
             relationship = "many-to-many") %>%
  rowwise() %>%
  mutate(shared_phenotype = phenotype_type %in% phenotype_types) %>%
  ungroup() %>%
  dplyr::filter(shared_phenotype)

concordant_genes <- unique(overlap_pairs$symbol)
universe_genes   <- unique(candidate_universe$symbol)
abc_genes        <- unique(abc_by_locus_phenotype$symbol)
flames_genes     <- unique(flames_by_locus_phenotype$symbol)

a <- length(concordant_genes)
b <- length(setdiff(abc_genes, concordant_genes))
c <- length(setdiff(flames_genes, concordant_genes))
d <- length(universe_genes) - a - b - c

ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

cat(sprintf("\n  2x2 table (locus-restricted background, n=%d candidate genes):\n", length(universe_genes)))
cat(sprintf("    ABC+FLAMES, SAME phenotype:      %d\n", a))
cat(sprintf("    ABC only (or mismatched phenotype): %d\n", b))
cat(sprintf("    FLAMES only (or mismatched phenotype): %d\n", c))
cat(sprintf("    Neither:                          %d\n", d))
cat(sprintf("  OR = %.2f, p = %.3e (one-sided, greater)\n", ft$estimate, ft$p.value))

fisher_summary <- tibble(
  n_candidate_universe = length(universe_genes),
  n_abc = length(abc_genes), n_flames_restricted = length(flames_genes),
  n_concordant = a, OR = ft$estimate, p_value = ft$p.value,
  concordant_genes = paste(sort(concordant_genes), collapse = ";")
)
write.csv(fisher_summary, file.path(out_dir, "fisher_test_result.csv"), row.names = FALSE)
write.csv(overlap_pairs %>% mutate(phenotype_types = map_chr(phenotype_types, paste, collapse = ";")),
          file.path(out_dir, "concordant_gene_locus_phenotype_pairs.csv"), row.names = FALSE)

# ================================================================
# SECTION 7: Locus-stratified permutation test (phenotype-aware)
# ================================================================
# For each locus, the phenotype-matched FLAMES set is held FIXED.
# The candidate pool for permutation is restricted to genes that
# are phenotype-eligible at that locus (i.e. genes physically in
# the locus AND that appear anywhere in the FLAMES phenotype-eligible
# set for that locus's phenotype(s)) is NOT required here -- instead
# we permute across the full positional candidate pool per locus
# (same as the geometric background) and ask how often a random
# ABC-sized draw would hit the FIXED phenotype-matched FLAMES set.
# This still tests "given ABC found n hits at this locus, how often
# would n random candidate genes match FLAMES's actual phenotype-
# matched genes by chance" -- i.e. permutes ABC's specificity while
# holding the real (phenotype-filtered) FLAMES evidence fixed.
ts("Running locus-stratified permutation test (phenotype-aware)...")

locus_candidates <- candidate_universe %>%
  group_by(locus_id) %>% summarise(candidates = list(unique(symbol)), .groups = "drop")

locus_abc_n <- abc_by_locus_phenotype %>% distinct(locus_id, symbol) %>%
  count(locus_id, name = "n_abc")

locus_flames_set <- flames_by_locus_phenotype %>%
  group_by(locus_id) %>% summarise(flames_set = list(unique(symbol)), .groups = "drop")

perm_input <- locus_candidates %>%
  left_join(locus_abc_n, by = "locus_id") %>%
  left_join(locus_flames_set, by = "locus_id") %>%
  mutate(n_abc = replace_na(n_abc, 0L)) %>%
  dplyr::filter(n_abc > 0)

observed_overlap <- a  # from Section 6

perm_overlaps <- replicate(N_PERM, {
  sum(map_dbl(seq_len(nrow(perm_input)), function(i) {
    cands <- perm_input$candidates[[i]]
    n_abc <- perm_input$n_abc[i]
    fl    <- perm_input$flames_set[[i]]
    if (is.null(fl)) fl <- character(0)
    if (length(cands) < n_abc) return(0)
    sampled_abc <- sample(cands, n_abc, replace = FALSE)
    length(intersect(sampled_abc, fl))
  }))
})

perm_p <- mean(perm_overlaps >= observed_overlap)
cat(sprintf("\n  Permutation test (N=%d): observed concordant overlap = %d\n", N_PERM, observed_overlap))
cat(sprintf("  Null distribution: mean = %.2f, sd = %.2f\n", mean(perm_overlaps), sd(perm_overlaps)))
cat(sprintf("  Empirical p-value = %.4f\n", perm_p))

write.csv(tibble(observed_overlap = observed_overlap, perm_mean = mean(perm_overlaps),
                 perm_sd = sd(perm_overlaps), perm_p = perm_p, n_perm = N_PERM),
          file.path(out_dir, "permutation_test_result.csv"), row.names = FALSE)

cat("\n=== DONE ===\n")
cat("Outputs in:", out_dir, "\n")
cat("  merged_loci.csv\n")
cat("  candidate_universe_by_locus.csv\n")
cat("  abc_gene_locus_phenotype.csv\n")
cat("  flames_gene_locus_phenotype.csv\n")
cat("  concordant_gene_locus_phenotype_pairs.csv\n")
cat("  fisher_test_result.csv\n")
cat("  permutation_test_result.csv\n")