# ===========================================================================
# Arealisation_leadSNP_LD_check.R
# Manuscript item: Arealisation validation, step 1
#
# What this does:
#   Finds the lead SNP of each arealisation gene in each hit parcel (from
#   Supplementary Table 1) and checks, via LDlink, that the parcels share
#   one signal rather than several independent ones. Needs a personal LDlink
#   token in the environment variable LDLINK_TOKEN.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Arealisation validation -- Script 1 of 3
# Lead SNP identification + cross-parcel LD consistency check
#
# PURPOSE
# For each arealization gene x phenotype_type (from Figure 6 Part 1's
# arealisation_lookup.csv), and for EVERY hit_parcel that gene was
# FLAMES-prioritised in for that phenotype, identify the lead
# (lowest p-value) SNP within the gene's +/-500kb window in that
# specific parcel's GWAS summary statistics. Then check whether the
# lead SNP is the same (or in high LD) across all of that gene's hit
# parcels -- i.e. is this really ONE genetic signal producing a
# spatial gradient, or several coincidentally-overlapping independent
# signals that each separately cleared FLAMES significance in
# different parcels.
#
# STAGE A: extract lead SNP per gene x phenotype_type x hit_parcel.
#          Pulled directly from Supplementary Table 1 (genome-wide
#          significant SNPs, UKB-only discovery GWAS) rather than raw
#          per-parcel sumstats -- ST1 already contains, per regional
#          phenotype (MRI_scale/MRI_measure/region), every genome-wide
#          significant SNP with CHR/POS/UKB_BETA/UKB_SE/UKB_P. No
#          cluster access needed for this stage at all.
#
# STAGE B: for every gene x phenotype_type with >=2 hit_parcels and
#          >1 DISTINCT lead SNP, query LDlink (LDmatrix) for pairwise
#          r^2 / D' between the distinct lead SNPs, and flag any pair
#          below the LD threshold. This stage answers the immediate
#          question: "where are we getting two SNPs NOT in LD across
#          parcels for the same gene?" -- nothing here yet decides
#          which single SNP to carry forward to Script 2; that's a
#          separate decision once we've seen how big this problem is.
#
# REQUIRES: LDlinkR package + a personal LDlink API token
#   (register free at https://ldlink.nih.gov/?tab=apiaccess)
#   install.packages("LDlinkR")

source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(readxl)
  library(LDlinkR)
})

ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# ---- Paths -- ADJUST FOR CLUSTER ----
base_dir     <- proj_dir
tables_dir   <- file.path(base_dir, "Outputs", "Figure6_arealisation", "tables")
# arealisation_lookup.csv = ALL FDR-significant spatial clusters from
# Figure 6 Part 1, BEFORE the Part 2 quality filter (n_hits<=60,
# LCC_fraction>=0.30, n_components<=8) -- this is NOT the bonafide
# arealization gene set and includes many more gene x phenotype_type
# rows than actually get reported as arealization genes.
lookup_path  <- file.path(tables_dir, "arealisation_lookup.csv")
# arealisation_filtered.csv = the 81 bonafide arealization genes
# (126 gene x phenotype_type rows) that survive the quality filter --
# this is Supplementary Table 8, and the set this script should
# actually be restricted to. It has no chr/window columns of its own,
# so we use it to filter arealisation_lookup.csv down to just these
# rows rather than as a standalone source.
filtered_path <- file.path(tables_dir, "arealisation_filtered.csv")
out_dir      <- file.path(base_dir, "Outputs", "Arealisation_leadSNP")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ST1 = Supplementary Table 1 (genome-wide significant SNPs, UKB
# discovery + ABCD + meta-analysis). We only use the UKB-only columns
# here (UKB_BETA/UKB_SE/UKB_P), per instruction -- this is the UKB-only
# GWAS, so the lead SNP should be defined on UKB significance, not the
# meta-analysed result.
supp_table_path <- file.path(proj_dir, "Data", "Manuscript_Supplementary_Tables_MRI_GWAS.xlsx")

# LDlink personal access token -- required, get one free at the URL
# in the header comment above.
LDLINK_TOKEN <- Sys.getenv("LDLINK_TOKEN")

# Population code for LDlink queries. UK Biobank discovery sample is
# predominantly European ancestry -- GBR is the closest 1000G
# population; swap for CEU or a custom multi-population code if you'd
# rather match a broader EUR reference.
LD_POP <- "GBR"

# LD threshold below which a pair of lead SNPs is treated as a
# genuinely distinct signal rather than the same underlying variant.
# Also used in Stage C to decide when two lead SNPs collapse into one
# index SNP vs. get carried forward as separate index SNPs.
R2_THRESHOLD <- 0.4

# Per-gene window overrides. The default arealisation window is
# +/-500kb (baked into arealisation_lookup.csv's start_window /
# end_window from Figure 6 Part 1) -- widen it here for specific genes
# where a manual check of ST1 showed real genome-wide significant SNPs
# just outside that default window rather than genuinely absent
# signal. MEIS1/SA: nearest ST1 SNPs in its 9 hit_parcels sit
# 582-696kb from the gene body (p down to 1e-19) -- all missed by the
# default 500kb window by 82-196kb. Add more genes here as needed;
# units are +/- bp from start_hg19/end_hg19 (the gene body itself, NOT
# from the default window).
WINDOW_OVERRIDES <- c(MEIS1 = 1e6)

# ---- Load ST1 (genome-wide significant SNPs, UKB-only columns) ----
# Header layout matches the rest of the supplementary tables: row 1
# title, row 2 legend, row 3 blank, row 4 column headers, data from
# row 5 -- i.e. skip = 3 for read_excel.
ts("Loading ST1...")
st1_regional <- readxl::read_excel(supp_table_path, sheet = "ST1", skip = 3) %>%
  dplyr::filter(MRI_scale == "regional") %>%
  dplyr::transmute(
    phenotype_type = MRI_measure,
    parcel          = region,
    CHR             = as.character(CHR),
    SNP             = SNP,
    BP              = as.numeric(POS),
    BETA            = as.numeric(UKB_BETA),
    SE              = as.numeric(UKB_SE),
    P               = as.numeric(UKB_P))
cat(sprintf("  ST1 regional rows (genome-wide significant SNPs): %d\n", nrow(st1_regional)))

# ============================================================
# STAGE A: lead SNP per gene x phenotype_type x hit_parcel
# ============================================================
ts("=== STAGE A: extracting lead SNPs per gene x parcel ===")

lookup_raw <- readr::read_csv(lookup_path, show_col_types = FALSE)

# Restrict to the bonafide 81-gene arealization set (arealisation_filtered.csv
# / Supplementary Table 8) -- semi-join on symbol + phenotype_type so we
# keep lookup's chr/window columns but only for rows that survived the
# Part 2 quality filter. hit_parcels is identical between the two files
# for a given symbol x phenotype_type (filtering drops whole rows, not
# individual parcels), so this join doesn't need to touch hit_parcels.
bonafide <- readr::read_csv(filtered_path, show_col_types = FALSE) %>%
  dplyr::distinct(symbol, phenotype_type)

cat(sprintf("  arealisation_lookup.csv: %d rows (pre-quality-filter)\n", nrow(lookup_raw)))
cat(sprintf("  arealisation_filtered.csv (bonafide): %d gene x phenotype_type rows, %d unique genes\n",
            nrow(bonafide), n_distinct(bonafide$symbol)))

lookup <- lookup_raw %>%
  dplyr::semi_join(bonafide, by = c("symbol", "phenotype_type")) %>%
  dplyr::mutate(hit_parcel_list = str_split(hit_parcels, ";"))

cat(sprintf("  After restricting to bonafide arealization genes: %d rows (%d unique genes)\n",
            nrow(lookup), n_distinct(lookup$symbol)))
if (n_distinct(lookup$symbol) != n_distinct(bonafide$symbol)) {
  cat(sprintf("  WARNING: expected %d unique genes after the join but got %d -- check for symbol/phenotype_type mismatches between the two files\n",
              n_distinct(bonafide$symbol), n_distinct(lookup$symbol)))
}

# Explicitly name any bonafide gene x phenotype_type rows that aren't
# in arealisation_lookup.csv AT ALL (not a join-key mismatch in this
# script -- the row is simply absent from the source table) -- these
# have no chr/window info to search against and so can never get a
# lead SNP, regardless of anything downstream.
missing_bonafide <- bonafide %>% dplyr::anti_join(lookup_raw, by = c("symbol", "phenotype_type"))
if (nrow(missing_bonafide) > 0) {
  cat(sprintf("  WARNING: %d bonafide gene x phenotype_type row(s) missing from arealisation_lookup.csv entirely (check Figure 6 Part 1 output for these):\n",
              nrow(missing_bonafide)))
  print(missing_bonafide)
}

# Gene body coordinates (hg19), one row per gene x phenotype_type --
# used in Stage C to rank multiple independent index SNPs by distance
# to the gene itself.
gene_windows <- lookup %>%
  dplyr::distinct(symbol, phenotype_type, chr, start_hg19, end_hg19) %>%
  dplyr::mutate(chr = as.character(chr),
               start_hg19 = as.numeric(start_hg19), end_hg19 = as.numeric(end_hg19))

# One row per gene x phenotype_type x individual hit_parcel
gene_parcel_rows <- lookup %>%
  dplyr::select(symbol, phenotype_type, chr, start_window, end_window,
                hit_parcel_list) %>%
  tidyr::unnest(hit_parcel_list) %>%
  dplyr::rename(parcel = hit_parcel_list) %>%
  dplyr::distinct()

cat(sprintf("  %d gene x phenotype_type x parcel combinations to look up\n",
            nrow(gene_parcel_rows)))

# Apply WINDOW_OVERRIDES: for the listed genes, replace the default
# +/-500kb start_window/end_window with start_hg19/end_hg19 (gene
# body) +/- the override's bp value, joined in from gene_windows.
if (length(WINDOW_OVERRIDES) > 0) {
  override_tbl <- tibble(symbol = names(WINDOW_OVERRIDES),
                         window_bp = unname(WINDOW_OVERRIDES))
  gene_parcel_rows <- gene_parcel_rows %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    dplyr::left_join(gene_windows, by = c("symbol", "phenotype_type", "chr")) %>%
    dplyr::left_join(override_tbl, by = "symbol") %>%
    dplyr::mutate(
      start_window = dplyr::if_else(!is.na(window_bp), start_hg19 - window_bp, start_window),
      end_window   = dplyr::if_else(!is.na(window_bp), end_hg19 + window_bp, end_window)) %>%
    dplyr::select(-start_hg19, -end_hg19, -window_bp)
  applied <- override_tbl %>% dplyr::filter(symbol %in% gene_parcel_rows$symbol)
  if (nrow(applied) > 0) {
    cat("  Window override(s) applied:\n")
    for (i in seq_len(nrow(applied))) {
      cat(sprintf("    %s: +/-%.0fkb (was +/-500kb default)\n", applied$symbol[i], applied$window_bp[i] / 1000))
    }
  }
}

# Vectorized join against ST1 instead of a per-row file lookup: match
# on phenotype_type + parcel + chr, keep SNPs whose POS falls inside
# the gene's window, then take the lowest-P SNP per gene x
# phenotype_type x parcel. Rows with NO genome-wide significant SNP
# in ST1 for that specific parcel/window get flagged explicitly
# rather than silently dropped -- this can legitimately happen since
# FLAMES prioritisation isn't necessarily gated on every single hit
# parcel individually clearing genome-wide significance on its own.
gene_parcel_rows_id <- gene_parcel_rows %>%
  dplyr::mutate(chr = as.character(chr), row_id = dplyr::row_number())

matched <- gene_parcel_rows_id %>%
  dplyr::inner_join(st1_regional, by = c("phenotype_type", "parcel", "chr" = "CHR"),
                    relationship = "many-to-many") %>%
  dplyr::filter(BP >= start_window, BP <= end_window)

lead_matched <- matched %>%
  dplyr::group_by(row_id) %>%
  dplyr::slice_min(P, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(row_id, symbol, phenotype_type, parcel, lead_snp = SNP,
                   chr, pos = BP, p = P, beta = BETA, se = SE, status = "ok",
                   # ST1's SNP column mixes rsIDs with chr:pos_ref_alt
                   # fallback notation (used whenever dbSNP has no rsID
                   # for that variant, mostly indels) -- LDlink can't
                   # reliably resolve the latter. CHR/pos are always
                   # clean regardless of which notation SNP uses, so
                   # build a coordinate-format ID here for LDlink to
                   # query instead of relying on lead_snp's format.
                   coord_id = paste0("chr", chr, ":", format(pos, scientific = FALSE, trim = TRUE)))

unmatched <- gene_parcel_rows_id %>%
  dplyr::filter(!row_id %in% lead_matched$row_id) %>%
  dplyr::transmute(row_id, symbol, phenotype_type, parcel,
                   lead_snp = NA_character_, chr, pos = NA_real_,
                   p = NA_real_, beta = NA_real_, se = NA_real_,
                   status = "no_significant_snp_in_window",
                   coord_id = NA_character_)

lead_snps <- dplyr::bind_rows(lead_matched, unmatched) %>%
  dplyr::select(-row_id) %>%
  dplyr::arrange(symbol, phenotype_type, parcel)

write_csv(lead_snps, file.path(out_dir, "lead_snps_per_gene_parcel.csv"))
cat(sprintf("  Done. %d rows written (%d ok, %d with no genome-wide significant SNP in window)\n",
            nrow(lead_snps), sum(lead_snps$status == "ok"),
            sum(lead_snps$status == "no_significant_snp_in_window")))

# Does every bonafide gene x phenotype_type combination have AT LEAST
# ONE genome-wide significant SNP underlying it (in ANY of its hit
# parcels)? Two distinct ways to fail this: (1) missing_bonafide above
# (no window to search at all), or (2) present, searched, but every
# single hit_parcel came back "no_significant_snp_in_window" -- FLAMES
# prioritised the gene without any of its hit parcels individually
# clearing genome-wide significance in ST1's UKB-only results (worth a
# manual look: could be a window-size issue, or FLAMES using a
# different significance basis than single-SNP genome-wide P).
zero_ok_genes <- bonafide %>%
  dplyr::anti_join(missing_bonafide, by = c("symbol", "phenotype_type")) %>%
  dplyr::left_join(
    lead_snps %>% dplyr::group_by(symbol, phenotype_type) %>%
      dplyr::summarise(n_ok = sum(status == "ok"), .groups = "drop"),
    by = c("symbol", "phenotype_type")) %>%
  dplyr::filter(is.na(n_ok) | n_ok == 0)
cat(sprintf("\n  %d bonafide gene x phenotype_type row(s) missing from arealisation_lookup.csv (no window to search)\n",
            nrow(missing_bonafide)))
cat(sprintf("  %d bonafide gene x phenotype_type row(s) searched but with NO genome-wide significant SNP in ANY hit_parcel:\n",
            nrow(zero_ok_genes)))
if (nrow(zero_ok_genes) > 0) print(zero_ok_genes %>% dplyr::select(symbol, phenotype_type))
cat(sprintf("  => %d / %d bonafide gene x phenotype_type rows have at least one underlying significant SNP.\n",
            nrow(bonafide) - nrow(missing_bonafide) - nrow(zero_ok_genes), nrow(bonafide)))
ts("=== STAGE A done ===")

# ============================================================
# STAGE B: cross-parcel LD consistency check (LDlink)
# ============================================================
ts("=== STAGE B: LD consistency check across hit_parcels ===")

# Uniqueness is assessed on coord_id (chr:pos), not the raw lead_snp
# label -- same genomic position is the same SNP regardless of which
# ID notation ST1 happened to record for it.
gene_lead_sets <- lead_snps %>%
  dplyr::filter(status == "ok", !is.na(coord_id)) %>%
  dplyr::group_by(symbol, phenotype_type) %>%
  dplyr::summarise(n_parcels = n_distinct(parcel),
                   distinct_snps = list(sort(unique(coord_id))),
                   n_distinct_snps = n_distinct(coord_id),
                   .groups = "drop") %>%
  dplyr::filter(n_parcels >= 2)

cat(sprintf("  %d gene x phenotype_type combos with >=2 hit_parcels\n", nrow(gene_lead_sets)))
cat(sprintf("  %d of those have >1 distinct lead SNP -- these need an LD check\n",
            sum(gene_lead_sets$n_distinct_snps > 1)))

needs_ld_check <- gene_lead_sets %>% dplyr::filter(n_distinct_snps > 1)

# Visibility into Stage A's status breakdown -- if needs_ld_check ends
# up empty, this tells you WHY: either everything genuinely came back
# "ok" with concordant lead SNPs (good), or a lot of parcels failed to
# resolve for a different reason (missing files / no SNPs in window),
# which would silently shrink gene_lead_sets instead.
cat("\n  Stage A status breakdown (per gene x phenotype_type x parcel row):\n")
print(table(lead_snps$status))
cat(sprintf("\n  needs_ld_check has %d rows.\n", nrow(needs_ld_check)))
if (nrow(needs_ld_check) == 0) {
  cat("  Nothing required an LD check -- every gene with >=2 resolved hit_parcels\n")
  cat("  already had a single, identical lead SNP across all of them.\n")
  cat("  If that's unexpected, check the status breakdown above for the\n")
  cat("  'no_significant_snp_in_window' count first.\n")
}

# Always returns a tibble with this exact column set, even on total
# failure -- an empty-but-correctly-typed 0-row tibble rather than
# NULL, so pmap_dfr() below can never collapse to a columnless result
# and crash the downstream group_by().
empty_ld_result <- function() {
  tibble(symbol = character(), phenotype_type = character(),
         snp_a = character(), snp_b = character(),
         r2 = double(), in_ld = logical())
}

check_ld_matrix <- function(symbol, phenotype_type, snps) {
  # `snps` here are coord_id strings ("chr5:128734008"), not ST1's raw
  # SNP labels -- LDmatrix() accepts chromosome-coordinate queries
  # directly (per-variant, resolved against its own reference panel),
  # which sidesteps ST1's rsID/chr:pos_ref_alt notation split entirely.
  # Row/col names in the returned matrix will typically be the
  # resolved rsID where dbSNP has one for that position, or the
  # coordinate string back if not -- either way they're valid LD-panel
  # identifiers, unlike ST1's fallback chr:pos_ref_alt notation which
  # LDlink can't parse at all for indels.
  snps <- unlist(snps)
  if (length(snps) < 2) return(empty_ld_result())
  res <- tryCatch(
    LDlinkR::LDmatrix(snps = snps, pop = LD_POP, r2d = "r2",
                      token = LDLINK_TOKEN, genome_build = "grch37"),
    error = function(e) {
      cat(sprintf("  LDlink query FAILED for %s (%s): %s\n",
                  symbol, phenotype_type, conditionMessage(e)))
      NULL
    })
  if (is.null(res) || nrow(res) == 0 || ncol(res) < 2) {
    cat(sprintf("  LDlink query returned no usable result for %s (%s)\n", symbol, phenotype_type))
    return(empty_ld_result())
  }
  mat <- tryCatch({
    m <- as.matrix(res[, -1]); rownames(m) <- res[[1]]; m
  }, error = function(e) NULL)
  # LDlink can still drop a coordinate if there's no variant recorded
  # at that exact position in the reference panel (rare once we're
  # querying by position rather than by a possibly-malformed ID
  # string) -- report exactly which one rather than crashing.
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) {
    resolved <- if (!is.null(mat)) rownames(mat) else character(0)
    missing_snps <- setdiff(snps, resolved)
    cat(sprintf("  LDlink resolved only %d/%d SNPs for %s (%s) -- unresolved: %s\n",
                length(resolved), length(snps), symbol, phenotype_type,
                paste(missing_snps, collapse = ", ")))
    return(empty_ld_result())
  }
  pairs <- combn(rownames(mat), 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pr) {
    tibble(symbol = symbol, phenotype_type = phenotype_type,
           snp_a = pr[1], snp_b = pr[2],
           r2 = mat[pr[1], pr[2]],
           in_ld = !is.na(mat[pr[1], pr[2]]) && mat[pr[1], pr[2]] >= R2_THRESHOLD)
  })
}

# Guard against the zero-row case explicitly -- pmap_dfr() over an
# empty input never calls the function, so it can't inherit
# empty_ld_result()'s column types from check_ld_matrix() and would
# otherwise return a bare 0x0 tibble, breaking the group_by() below.
if (nrow(needs_ld_check) == 0) {
  ld_results <- empty_ld_result()
} else {
  ld_results <- purrr::pmap_dfr(
    needs_ld_check %>% dplyr::select(symbol, phenotype_type, distinct_snps),
    function(symbol, phenotype_type, distinct_snps) {
      cat(sprintf("  Checking LD: %s / %s (%d SNPs)\n", symbol, phenotype_type, length(distinct_snps)))
      check_ld_matrix(symbol, phenotype_type, distinct_snps)
    })
}

write_csv(ld_results, file.path(out_dir, "lead_snp_LD_pairs.csv"))

# Summary: genes/measures where AT LEAST ONE pair of lead SNPs across
# hit_parcels is NOT in high LD -- these are the ones worth looking
# at manually before deciding how Script 2 should pick a
# representative SNP.
discordant_genes <- ld_results %>%
  dplyr::group_by(symbol, phenotype_type) %>%
  dplyr::summarise(n_pairs = n(), n_discordant = sum(!in_ld, na.rm = TRUE),
                   min_r2 = min(r2, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(n_discordant > 0) %>%
  dplyr::arrange(min_r2)

write_csv(discordant_genes, file.path(out_dir, "discordant_gene_measure_pairs.csv"))

cat(sprintf("\n  %d gene x phenotype_type combos have at least one lead-SNP pair NOT in LD (r2 < %.1f):\n",
            nrow(discordant_genes), R2_THRESHOLD))
print(discordant_genes)

ts("=== STAGE B done ===")

# ============================================================
# STAGE C: index SNP(s) per gene x phenotype_type
# ============================================================
# The deliverable Script 2 actually needs: ideally ONE index SNP per
# gene x phenotype_type combination. Where a gene's lead SNPs across
# hit_parcels genuinely split into 2+ signals not in high LD (r2 <
# R2_THRESHOLD), we carry ALL of those independent signals forward
# rather than arbitrarily picking one -- each gets interrogated
# separately downstream.
#
# Approach: build a graph per gene x phenotype_type where nodes are
# that gene's distinct lead SNPs (by coord_id) and edges connect pairs
# LDlink confirmed are in LD (r2 >= R2_THRESHOLD). Connected components
# = LD blocks; one representative index SNP is picked per block (the
# member with the lowest p-value across its supporting parcels).
#
# Pairs LDlink never resolved contribute NO edge -- an unresolved SNP
# is therefore left as its own singleton unless a resolved edge
# elsewhere links it to another node. This is deliberately
# conservative: absence of LD evidence is not treated as evidence of
# LD, so an unresolved SNP still gets interrogated as its own index
# SNP downstream rather than being silently folded into another
# cluster. Gene x phenotype_type groups with any unresolved pair are
# flagged (ld_complete = FALSE) so this can be reviewed later.
ts("=== STAGE C: selecting index SNP(s) per gene x phenotype_type ===")

# One row per distinct lead SNP within a gene x phenotype_type, with
# the parcels it was the lead SNP for and its best p-value across them.
snp_summary <- lead_snps %>%
  dplyr::filter(status == "ok", !is.na(coord_id)) %>%
  dplyr::group_by(symbol, phenotype_type, coord_id, lead_snp, chr, pos) %>%
  dplyr::summarise(n_parcels = n_distinct(parcel),
                   parcels = paste(sort(unique(parcel)), collapse = ";"),
                   min_p = min(p, na.rm = TRUE),
                   .groups = "drop")

# Map LDlink's returned SNP label back to coord_id. Every non-indel
# lead SNP in ST1 already carries a real rsID, and LDlink's coordinate
# queries resolve to that same rsID, so a direct label match is
# reliable; indels (which never resolve, see Stage B) simply won't
# match here and their pairs are absent from ld_edges below, which is
# exactly the "unresolved" case this stage is designed to handle
# conservatively.
label_lookup <- snp_summary %>%
  dplyr::distinct(symbol, phenotype_type, coord_id, lead_snp)

ld_edges <- ld_results %>%
  dplyr::inner_join(label_lookup, by = c("symbol", "phenotype_type", "snp_a" = "lead_snp")) %>%
  dplyr::rename(coord_a = coord_id) %>%
  dplyr::inner_join(label_lookup, by = c("symbol", "phenotype_type", "snp_b" = "lead_snp")) %>%
  dplyr::rename(coord_b = coord_id) %>%
  dplyr::mutate(linked = !is.na(r2) & r2 >= R2_THRESHOLD) %>%
  dplyr::select(symbol, phenotype_type, coord_a, coord_b, r2, linked)

# Simple union-find (disjoint set) over a gene x phenotype_type's
# distinct coord_ids, using only the "linked" (r2 >= R2_THRESHOLD)
# edges available for that group. No igraph dependency needed for
# something this small.
cluster_coords <- function(coords, edges) {
  n <- length(coords)
  parent <- seq_len(n)
  find <- function(x) { while (parent[x] != x) x <- parent[x]; x }
  if (nrow(edges) > 0) {
    for (i in seq_len(nrow(edges))) {
      if (isTRUE(edges$linked[i])) {
        ia <- match(edges$coord_a[i], coords)
        ib <- match(edges$coord_b[i], coords)
        if (!is.na(ia) && !is.na(ib)) {
          ra <- find(ia); rb <- find(ib)
          if (ra != rb) parent[ra] <- rb
        }
      }
    }
  }
  roots <- vapply(seq_len(n), find, integer(1))
  tibble(coord_id = coords, cluster_id = match(roots, unique(roots)))
}

gene_measure_keys <- snp_summary %>% dplyr::distinct(symbol, phenotype_type)

clusters <- purrr::pmap_dfr(gene_measure_keys, function(symbol, phenotype_type) {
  coords <- snp_summary %>%
    dplyr::filter(symbol == !!symbol, phenotype_type == !!phenotype_type) %>%
    dplyr::pull(coord_id)
  edges <- ld_edges %>%
    dplyr::filter(symbol == !!symbol, phenotype_type == !!phenotype_type)
  n_pairs_expected <- length(coords) * (length(coords) - 1) / 2
  n_pairs_resolved <- nrow(edges)
  cluster_coords(coords, edges) %>%
    dplyr::mutate(symbol = symbol, phenotype_type = phenotype_type,
                  ld_complete = n_pairs_resolved >= n_pairs_expected)
})

index_snps <- clusters %>%
  dplyr::left_join(snp_summary, by = c("symbol", "phenotype_type", "coord_id")) %>%
  dplyr::group_by(symbol, phenotype_type, cluster_id) %>%
  dplyr::mutate(cluster_size = dplyr::n(),
               cluster_members = paste(sort(unique(lead_snp)), collapse = ";")) %>%
  dplyr::slice_min(min_p, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  # Distance from each cluster's representative SNP to the gene body
  # itself (0 if the SNP falls inside start_hg19-end_hg19, else the
  # distance to the nearer boundary) -- used below to flag, when a
  # gene has 2+ genuinely independent loci, which one sits closest to
  # the gene (most plausibly the signal actually driving the FLAMES
  # call) vs. the other(s), which are farther away and more likely to
  # implicate a different nearby gene.
  dplyr::left_join(gene_windows, by = c("symbol", "phenotype_type", "chr")) %>%
  dplyr::mutate(distance_to_gene = dplyr::case_when(
    is.na(start_hg19) | is.na(end_hg19) ~ NA_real_,
    pos >= start_hg19 & pos <= end_hg19 ~ 0,
    pos < start_hg19 ~ start_hg19 - pos,
    pos > end_hg19 ~ pos - end_hg19,
    TRUE ~ NA_real_)) %>%
  dplyr::group_by(symbol, phenotype_type) %>%
  dplyr::mutate(n_index_snps = dplyr::n()) %>%
  # locus_rank 1 = closest to the gene body (ties broken by p-value);
  # is_primary_locus flags that one. For single-index-SNP genes this
  # is trivially rank 1 / primary.
  dplyr::arrange(distance_to_gene, min_p, .by_group = TRUE) %>%
  dplyr::mutate(locus_rank = dplyr::row_number(),
               is_primary_locus = locus_rank == 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(symbol, phenotype_type, status = "ok",
                   index_snp = lead_snp, coord_id, chr, pos,
                   distance_to_gene, locus_rank, is_primary_locus,
                   min_p, n_supporting_parcels = n_parcels, supporting_parcels = parcels,
                   cluster_size, cluster_members, n_index_snps, ld_complete)

# Explicitly flag the bonafide gene x phenotype_type rows that have NO
# index SNP at all, rather than letting them just be absent from the
# file -- one row each, same column layout, status names them so
# they're impossible to miss and match Stage A's own status vocabulary.
flagged_missing <- missing_bonafide %>%
  dplyr::transmute(symbol, phenotype_type, status = "missing_from_lookup",
                   index_snp = NA_character_, coord_id = NA_character_,
                   chr = NA_character_, pos = NA_real_,
                   distance_to_gene = NA_real_, locus_rank = NA_integer_,
                   is_primary_locus = NA, min_p = NA_real_,
                   n_supporting_parcels = NA_integer_, supporting_parcels = NA_character_,
                   cluster_size = NA_integer_, cluster_members = NA_character_,
                   n_index_snps = 0L, ld_complete = NA)

flagged_no_sig <- zero_ok_genes %>%
  dplyr::transmute(symbol, phenotype_type, status = "no_significant_snp_in_window",
                   index_snp = NA_character_, coord_id = NA_character_,
                   chr = NA_character_, pos = NA_real_,
                   distance_to_gene = NA_real_, locus_rank = NA_integer_,
                   is_primary_locus = NA, min_p = NA_real_,
                   n_supporting_parcels = NA_integer_, supporting_parcels = NA_character_,
                   cluster_size = NA_integer_, cluster_members = NA_character_,
                   n_index_snps = 0L, ld_complete = NA)

index_snps <- dplyr::bind_rows(index_snps, flagged_missing, flagged_no_sig) %>%
  dplyr::arrange(symbol, phenotype_type, locus_rank)

write_csv(index_snps, file.path(out_dir, "index_snps_per_gene_measure.csv"))
cat(sprintf("  Flagged %d bonafide gene x phenotype_type row(s) with no index SNP (%d missing_from_lookup, %d no_significant_snp_in_window)\n",
            nrow(flagged_missing) + nrow(flagged_no_sig), nrow(flagged_missing), nrow(flagged_no_sig)))

resolved_index_snps <- index_snps %>% dplyr::filter(status == "ok")
n_gene_measure_combos <- resolved_index_snps %>% dplyr::distinct(symbol, phenotype_type) %>% nrow()
cat(sprintf("  %d gene x phenotype_type combos -> %d index SNP(s) total\n",
            n_gene_measure_combos, nrow(resolved_index_snps)))
n_multi_index <- resolved_index_snps %>%
  dplyr::distinct(symbol, phenotype_type, n_index_snps) %>%
  dplyr::filter(n_index_snps > 1) %>%
  nrow()
cat(sprintf("  %d gene x phenotype_type combos need >1 index SNP (genuinely distinct signals, r2 < %.1f)\n",
            n_multi_index, R2_THRESHOLD))
cat(sprintf("  of those, %d index SNPs are flagged locus_rank==1/is_primary_locus (closest to the gene body)\n",
            sum(resolved_index_snps$is_primary_locus & resolved_index_snps$n_index_snps > 1)))
cat(sprintf("  and %d are flagged as farther/secondary loci (locus_rank > 1) -- see distance_to_gene\n",
            sum(!resolved_index_snps$is_primary_locus)))
n_incomplete <- resolved_index_snps %>% dplyr::distinct(symbol, phenotype_type, ld_complete) %>%
  dplyr::filter(!ld_complete) %>% nrow()
cat(sprintf("  %d gene x phenotype_type combos have at least one unresolved LD pair (ld_complete = FALSE) --\n",
            n_incomplete))
cat("  review these manually; their index SNP count may be an overestimate since some\n")
cat("  of their 'independent' clusters could turn out to be in LD if LDlink had resolved them.\n")

ts("=== STAGE C done ===")
cat("\nOutputs:\n")
cat("  ", file.path(out_dir, "lead_snps_per_gene_parcel.csv"), "\n")
cat("  ", file.path(out_dir, "lead_snp_LD_pairs.csv"), "\n")
cat("  ", file.path(out_dir, "discordant_gene_measure_pairs.csv"), "\n")
cat("  ", file.path(out_dir, "index_snps_per_gene_measure.csv"), "  <- Script 2 input\n")
