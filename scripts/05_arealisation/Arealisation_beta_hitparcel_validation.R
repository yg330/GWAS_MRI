# ===========================================================================
# Arealisation_beta_hitparcel_validation.R
# Manuscript item: Arealisation validation (beta QC)
#
# What this does:
#   Two checks on each gene x measure: (A) Wilcoxon test of |beta| inside vs
#   outside the hit parcels, (B) spin-test Spearman correlation of |beta|
#   with distance to the nearest hit parcel (expected negative).
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Do arealization genes' lead-SNP betas actually line up with the
# hit_parcels they were originally identified in (Figure 6 Part 1 /
# arealisation_filtered.csv), or is that just an artefact of the
# spatial-clustering test used to call them "arealization genes" in
# the first place?
#
# Two complementary checks, per gene x phenotype_type (primary
# locus only, same index SNPs as PC1_vs_arealisation_beta_correlation.R):
#
#   TEST A -- Within vs outside hit_parcels
#     Is |beta| bigger inside the gene's own hit_parcels than in the
#     rest of cortex? Wilcoxon rank-sum test, |beta| in hit_parcels
#     vs the other (180 - n_hits) parcels.
#
#   TEST B -- Distance-based decay
#     Does |beta| fall off with geodesic distance from the nearest
#     hit_parcel? Spearman correlation of |beta| vs distance-to-
#     nearest-hit-parcel (0 for parcels inside the hit set), with a
#     spin-test null (same spherical rotation method as
#     PC1_spatial_correlation.R / PC1_vs_arealisation_beta_correlation.R,
#     for consistency) -- expect a NEGATIVE rho (farther = weaker).
#
# Inputs:
#   Outputs/Figure6_arealisation/tables/arealisation_filtered.csv  -- hit_parcels per gene x measure
#   Outputs/Arealisation_leadSNP/index_snps_per_gene_measure.csv   -- primary-locus flags
#   Outputs/Arealisation_leadSNP/Arealisation_beta_outputs/{measure}.csv -- per-parcel betas
#   Data/Sydnor_maps_scripts/{glasser_regions.csv, glasser_spherical_coordinates.txt,
#                             rotate.parcellation.R, perm.sphere.p.R}

source("config.R")

if (!requireNamespace("clue", quietly = TRUE)) install.packages("clue")

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---------------- USER INPUTS ----------------
filtered_path   <- "Outputs/Figure6_arealisation/tables/arealisation_filtered.csv"
index_snp_path  <- "Outputs/Arealisation_leadSNP/index_snps_per_gene_measure.csv"
beta_dir        <- "Outputs/Arealisation_leadSNP/Arealisation_beta_outputs"
sydnor_dir      <- "Data/Sydnor_maps_scripts"
out_dir         <- "Outputs/Arealisation_leadSNP/hitparcel_validation"

n_perm <- 10000
seed   <- 42

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

glasser_regions     <- read_csv(file.path(sydnor_dir, "glasser_regions.csv"), show_col_types = FALSE)
glasser_regions_180 <- glasser_regions$region

# ============================================================
# STEP 1 -- hit_parcels per gene x phenotype_type
# ============================================================
areal <- read_csv(filtered_path, show_col_types = FALSE) %>%
  transmute(symbol, phenotype_type, hit_parcels = str_split(hit_parcels, ";"), n_hits)

# ============================================================
# STEP 2 -- primary-locus beta profile per gene x measure
# (identical loading logic to PC1_vs_arealisation_beta_correlation.R)
# ============================================================
index_snps <- read_csv(index_snp_path, show_col_types = FALSE)
primary_loci <- index_snps %>%
  filter(status == "ok", is_primary_locus == TRUE) %>%
  distinct(symbol, phenotype_type, locus_rank)

beta_cache <- list()
get_beta_profile <- function(pheno) {
  if (!is.null(beta_cache[[pheno]])) return(beta_cache[[pheno]])
  path <- file.path(beta_dir, paste0(pheno, ".csv"))
  if (!file.exists(path)) return(NULL)
  df <- read_csv(path, show_col_types = FALSE)
  beta_cache[[pheno]] <<- df
  df
}

gene_profiles <- primary_loci %>%
  pmap_dfr(function(symbol, phenotype_type, locus_rank) {
    df <- get_beta_profile(phenotype_type)
    if (is.null(df)) return(NULL)
    sub <- df %>% filter(symbol == !!symbol, locus_rank == !!locus_rank, status == "ok")
    if (nrow(sub) == 0) return(NULL)
    sub %>% transmute(symbol, phenotype_type, parcel, beta)
  })

cat(sprintf("Beta profiles loaded for %d gene x measure pairs\n",
            n_distinct(gene_profiles[c("symbol","phenotype_type")])))

# ============================================================
# STEP 3 -- spherical coordinates + pairwise geodesic (angular) distance
# ============================================================
sphere_path <- file.path(sydnor_dir, "glasser_spherical_coordinates.txt")
stopifnot(file.exists(sphere_path))

sphere_full <- as.matrix(read.table(sphere_path, header = FALSE))
sphere_lh <- sphere_full[1:180, , drop = FALSE]
rownames(sphere_lh) <- glasser_regions_180
sphere_lh_ordered <- sphere_lh[glasser_regions_180, , drop = FALSE]

# angular (great-circle) distance between all parcel pairs, in radians
unit <- sphere_lh_ordered / sqrt(rowSums(sphere_lh_ordered^2))
cos_sim <- unit %*% t(unit)
cos_sim <- pmin(pmax(cos_sim, -1), 1)   # numerical safety
geo_dist <- acos(cos_sim)
rownames(geo_dist) <- colnames(geo_dist) <- glasser_regions_180
diag(geo_dist) <- 0

dist_to_nearest_hit <- function(hit_parcels) {
  hit_parcels <- intersect(hit_parcels, glasser_regions_180)
  if (length(hit_parcels) == 0) return(rep(NA_real_, 180))
  sub <- geo_dist[, hit_parcels, drop = FALSE]
  apply(sub, 1, min)
}

# ============================================================
# STEP 4 -- spin test setup (true spherical, same infra as before)
# ============================================================
rotate_parcellation_path <- file.path(sydnor_dir, "rotate.parcellation.R")
perm_sphere_p_path       <- file.path(sydnor_dir, "perm.sphere.p.R")
source(rotate_parcellation_path)
source(perm_sphere_p_path)

cat(sprintf("Generating %d spherical rotations...\n", n_perm))
set.seed(seed)
perm_full <- rotate.parcellation(coord.l = sphere_lh_ordered, coord.r = sphere_lh_ordered,
                                  nrot = n_perm, method = "vasa")
spins <- perm_full[1:180, , drop = FALSE]

spin_test <- function(map_vec, beta_vec, spins, method = "spearman") {
  obs <- suppressWarnings(cor(map_vec, beta_vec, method = method, use = "complete.obs"))
  null <- apply(spins, 2, function(idx) {
    suppressWarnings(cor(map_vec, beta_vec[idx], method = method, use = "complete.obs"))
  })
  p_spin <- mean(abs(null) >= abs(obs), na.rm = TRUE)
  list(rho = obs, p_spin = p_spin)
}

# ============================================================
# STEP 5 -- run both tests per gene x measure
# ============================================================
gene_measure_pairs <- gene_profiles %>% distinct(symbol, phenotype_type)
cat(sprintf("Running Test A + Test B for %d gene x measure pairs...\n", nrow(gene_measure_pairs)))

results <- gene_measure_pairs %>%
  pmap_dfr(function(symbol, phenotype_type) {
    beta_sub <- gene_profiles %>% filter(symbol == !!symbol, phenotype_type == !!phenotype_type)
    beta_vec <- beta_sub$beta[match(glasser_regions_180, beta_sub$parcel)]

    hit_row <- areal %>% filter(symbol == !!symbol, phenotype_type == !!phenotype_type)
    if (nrow(hit_row) == 0) return(NULL)
    hit_parcels <- hit_row$hit_parcels[[1]]
    is_hit <- glasser_regions_180 %in% hit_parcels

    n_ok <- sum(!is.na(beta_vec))
    if (n_ok < 30 || sum(is_hit, na.rm = TRUE) < 2 || sum(!is_hit, na.rm = TRUE) < 2) {
      return(tibble(symbol, measure = phenotype_type, n_hits = sum(is_hit), n_parcels = n_ok,
                     medabs_beta_hit = NA_real_, medabs_beta_other = NA_real_,
                     wilcox_p = NA_real_, distdecay_rho = NA_real_, distdecay_p_spin = NA_real_))
    }

    abs_beta <- abs(beta_vec)

    # TEST A: within vs outside hit_parcels
    wt <- suppressWarnings(wilcox.test(abs_beta[is_hit], abs_beta[!is_hit]))

    # TEST B: distance-to-nearest-hit vs |beta|, spin-test null
    dist_vec <- dist_to_nearest_hit(hit_parcels)
    r <- spin_test(dist_vec, abs_beta, spins, method = "spearman")

    tibble(symbol, measure = phenotype_type,
           n_hits = sum(is_hit), n_parcels = n_ok,
           medabs_beta_hit   = median(abs_beta[is_hit], na.rm = TRUE),
           medabs_beta_other = median(abs_beta[!is_hit], na.rm = TRUE),
           wilcox_p = wt$p.value,
           distdecay_rho = r$rho, distdecay_p_spin = r$p_spin)
  }) %>%
  filter(!is.na(wilcox_p) | !is.na(distdecay_p_spin)) %>%
  mutate(wilcox_p_fdr = p.adjust(wilcox_p, method = "BH"),
         distdecay_p_spin_fdr = p.adjust(distdecay_p_spin, method = "BH")) %>%
  arrange(wilcox_p)

write_csv(results, file.path(out_dir, "hitparcel_validation_results.csv"))

cat(sprintf("\nDone: %d gene x measure pairs tested.\n", nrow(results)))
cat(sprintf("Test A (within > outside hit_parcels), FDR<0.05: %d / %d\n",
            sum(results$wilcox_p_fdr < 0.05, na.rm = TRUE), sum(!is.na(results$wilcox_p_fdr))))
cat(sprintf("Test B (distance decay, negative rho, FDR<0.05): %d / %d\n",
            sum(results$distdecay_p_spin_fdr < 0.05 & results$distdecay_rho < 0, na.rm = TRUE),
            sum(!is.na(results$distdecay_p_spin_fdr))))
cat(sprintf("Results written to %s/hitparcel_validation_results.csv\n", out_dir))
print(results, n = Inf)

# ============================================================
# STEP 6 -- merge the within/outside-parcel beta comparison and the
# distance-decay result (both computed from each gene's closest/
# primary-locus lead SNP, per STEP 2 above) back into the master
# arealisation_filtered.csv spreadsheet, so hit_parcels, spatial-
# clustering stats, and this beta-level validation live in one table.
# A backup of the pre-merge file is written first since
# arealisation_filtered.csv is regenerated by
# Final_Figure6_arealizationscript.R and shouldn't silently drift.
# ============================================================
areal_full <- read_csv(filtered_path, show_col_types = FALSE)

file.copy(filtered_path,
          str_replace(filtered_path, "\\.csv$", "_pre_hitparcelmerge_backup.csv"),
          overwrite = TRUE)

areal_merged <- areal_full %>%
  left_join(
    results %>%
      dplyr::select(symbol, phenotype_type = measure, n_hits_beta = n_hits, n_parcels,
                     medabs_beta_hit, medabs_beta_other, wilcox_p, wilcox_p_fdr,
                     distdecay_rho, distdecay_p_spin, distdecay_p_spin_fdr),
    by = c("symbol", "phenotype_type")
  )

n_unmatched <- sum(is.na(areal_merged$wilcox_p_fdr))
if (n_unmatched > 0) {
  cat(sprintf("\nNOTE: %d / %d arealisation_filtered.csv rows had no matching beta-validation result (not merged).\n",
              n_unmatched, nrow(areal_merged)))
}

write_csv(areal_merged, filtered_path)
cat(sprintf("\narealisation_filtered.csv updated in place with within/outside-parcel beta (medabs_beta_hit vs medabs_beta_other, wilcox_p_fdr) and distance-decay (distdecay_rho, distdecay_p_spin_fdr) columns.\n"))
cat(sprintf("Backup of the pre-merge file written to %s\n",
            str_replace(filtered_path, "\\.csv$", "_pre_hitparcelmerge_backup.csv")))
