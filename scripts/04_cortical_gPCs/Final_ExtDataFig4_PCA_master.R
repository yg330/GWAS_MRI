# ===========================================================================
# Final_ExtDataFig4_PCA_master.R
# Manuscript item: Figure 4 and Supplementary Tables 8-10
#
# What this does:
#   Main gPC script. Does PCA on 1 - rg within each cortical measure to get
#   gPC1-3, writes loadings and variance explained, draws the brain maps,
#   runs spin tests against reference maps (10,000 rotations), and tests the
#   arealisation-gene lead-SNP beta profiles against gPC1-3.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Master script -- rebuild of Extended Data Fig 4 (within-measure
# genetic PC1/PC2/PC3), plus:
#   - cross-measure similarity of PC1/PC2/PC3 loadings
#   - PC1/PC2/PC3 vs restricted Sydnor reference maps + AHBA C1-C3
#     transcriptomic gradients (spin test)
#   - arealization gene lead-SNP beta profiles vs PC1/PC2/PC3 and
#     the same reference maps (spin test)
#
# IMPORTANT -- "PC1/PC2/PC3" here means the pheno_pca_pc1()-style
# construction (prcomp(1 - rg, center=TRUE) on the within-measure,
# diagonal-zeroed genetic correlation submatrix), i.e. the SAME
# construction that generates the existing Extended Data Fig 4 brain
# maps. This is NOT the same object as "PC1_GWAS" (which, per this
# project's other diagnostics, was built via Furtjes'-style eigen()
# decomposition of the raw correlation matrix and is close to the
# general/shared genetic factor for most measures). Nothing in this
# script requires a new GWAS -- everything here operates on loadings
# (this script) and already-existing lead-SNP betas (arealization
# genes), never on new individual-level scores or summary statistics.
#
# "RESTRICTED" parcel set -- AHBA gene-expression coverage (Hansen/
# Vos et al. 2024 DME components, sheet "1a" of the supplementary
# workbook) only spans a subset of the 180 Glasser left-hemisphere
# parcels (donor-tissue-sampling coverage limits). To keep PC1/2/3,
# Sydnor, and AHBA comparisons on a like-for-like footing, every
# comparison in this script (PCs vs maps, arealization betas vs PCs/
# maps) is restricted to the intersection of parcels covered by ALL
# THREE families (genetic PCs x Sydnor x AHBA). That intersection is
# written out (restricted_parcel_list.csv) so you can sanity-check
# it -- if this isn't what you meant by "restricted Sydnor maps",
# the fix is a one-line change to `restricted_parcels` in STEP 5.
#
# ORDERING -- all fast/deterministic steps (redo PCs, lollipop,
# cross-measure similarity, brain maps, reference-map + arealization-
# beta prep) run FIRST. The slow spherical spin-permutation tests
# (STEP 7-9) are all at the END, so you can debug/rerun everything
# before them without waiting on permutations.
#
# Inputs:
#   ~/Downloads/rstudio-export (3)/rg_matrix.rds                    -- same as Final_Figure3_script.R
#   Data/Sydnor_maps_scripts/{glasser_regions.csv, brainmaps_glasser.csv,
#     Sensorimotor_Association_Axis_AverageRanks.csv,
#     glasser_spherical_coordinates.txt, rotate.parcellation.R, perm.sphere.p.R}
#   Data/AHBA_DME_Hansen2024/41593_2024_1624_MOESM3_ESM.xlsx (sheet "1a") --
#     Hansen et al. 2024 Nat Neurosci differential module expression components
#   Outputs/Arealisation_leadSNP/index_snps_per_gene_measure.csv
#   Outputs/Arealisation_leadSNP/Arealisation_beta_outputs/{measure}.csv
#
# Outputs (all under out_dir, see below):
#   pc123_loadings.csv, pc123_variance_explained.csv
#   panelA_lollipop.pdf
#   panelC_PC1/2/3_similarity.pdf, cross_measure_similarity.rds
#   mainFig_PC1_brains.pdf, suppFig_PC2_brains.pdf, suppFig_PC3_brains.pdf
#   reference_maps_restricted.csv, restricted_parcel_list.csv
#   arealisation_beta_profiles.csv
#   PC123_vs_reference_maps_spintest.csv
#   arealisation_beta_vs_PC123_and_maps_spintest.csv

source("config.R")

if (!requireNamespace("clue", quietly = TRUE)) install.packages("clue")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggsegGlasser)
  library(cowplot)
  library(sf)
})

# ---------------- USER INPUTS ----------------
rg_path              <- file.path(proj_dir, "Data", "rg_matrix.rds")
sydnor_dir           <- "Data/Sydnor_maps_scripts"
ahba_path            <- "Data/AHBA_DME_Hansen2024/41593_2024_1624_MOESM3_ESM.xlsx"
areal_index_snp_path <- "Outputs/Arealisation_leadSNP/index_snps_per_gene_measure.csv"
areal_beta_dir       <- "Outputs/Arealisation_leadSNP/Arealisation_beta_outputs"
out_dir              <- "Outputs/ExtDataFig4_PCA_master"

n_perm <- 10000
seed   <- 42
all_phenos <- c("CT","FA","FI","GC","GMV","IC","ICVF","ISOVF","MC","MD","OD","SA")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

glasser_regions     <- read_csv(file.path(sydnor_dir, "glasser_regions.csv"), show_col_types = FALSE)
glasser_regions_180 <- glasser_regions$region
stopifnot(length(glasser_regions_180) == 180, !anyDuplicated(glasser_regions_180))

# ============================================================
# STEP 1 -- Redo PC1/PC2/PC3 from rg_matrix, all 12 measures
# (identical construction to pheno_pca_pc1() in Final_Figure3_script.R,
# extended to keep components 1-3 instead of just component 1)
# ============================================================
cat("==== STEP 1: computing PC1-3 per measure ====\n")
cat("Loading rg_matrix...\n")
rg_matrix <- readRDS(rg_path)
clean_names <- function(x) gsub("_as_pheno.*", "", gsub("GWAS_", "", x))
rownames(rg_matrix) <- clean_names(rownames(rg_matrix))
colnames(rg_matrix) <- clean_names(colnames(rg_matrix))

trait_names <- colnames(rg_matrix)
get_pheno_traits <- function(pheno) trait_names[grepl(paste0("^regional_", pheno, "_"), trait_names)]
extract_parcel   <- function(trait) paste(strsplit(trait, "_")[[1]][-(1:2)], collapse = "_")

pheno_pca_pc123 <- function(pheno) {
  traits <- get_pheno_traits(pheno)
  if (length(traits) < 4) return(NULL)

  parcels_used <- sapply(traits, extract_parcel)
  bad <- setdiff(parcels_used, glasser_regions_180)
  if (length(bad) > 0) {
    stop(sprintf("pheno_pca_pc123('%s'): %d extracted parcel name(s) not in the 180 Glasser\n",
                  pheno, length(bad)),
         "cortical parcels: ", paste(bad, collapse = ", "))
  }

  sub_mat <- rg_matrix[traits, traits]
  diag(sub_mat) <- 0
  pca_res <- prcomp(1 - sub_mat, center = TRUE, scale. = FALSE)

  n_pc <- min(3, ncol(pca_res$x))
  var_exp <- summary(pca_res)$importance["Proportion of Variance", 1:n_pc]
  scores  <- pca_res$x[, 1:n_pc, drop = FALSE]
  colnames(scores) <- paste0("PC", 1:n_pc)
  for (j in 1:n_pc) if (mean(scores[, j]) < 0) scores[, j] <- -scores[, j]
  rownames(scores) <- sapply(traits, extract_parcel)

  list(scores = scores, var_exp = var_exp)
}

pc123_res <- setNames(lapply(all_phenos, pheno_pca_pc123), all_phenos)
for (ph in names(pc123_res)) {
  if (is.null(pc123_res[[ph]])) {
    cat(sprintf("  %s: skipped (fewer than 4 traits)\n", ph))
  } else {
    ve <- pc123_res[[ph]]$var_exp * 100
    cat(sprintf("  %s: PC1=%.1f%%  PC2=%.1f%%  PC3=%.1f%%\n", ph, ve[1], ve[2], ve[3]))
  }
}
pc123_res <- Filter(Negate(is.null), pc123_res)

pc123_long <- map_dfr(names(pc123_res), function(ph) {
  sc <- pc123_res[[ph]]$scores
  tibble(parcel = rownames(sc), measure = ph,
         PC1 = sc[, "PC1"], PC2 = sc[, "PC2"], PC3 = sc[, "PC3"])
})
stopifnot(all(pc123_long$parcel %in% glasser_regions_180))
write_csv(pc123_long, file.path(out_dir, "pc123_loadings.csv"))

var_exp_df <- map_dfr(names(pc123_res), function(ph) {
  ve <- pc123_res[[ph]]$var_exp
  tibble(measure = ph, PC = paste0("PC", 1:3), pct_variance = as.numeric(ve[1:3]) * 100)
})
write_csv(var_exp_df, file.path(out_dir, "pc123_variance_explained.csv"))
cat("Wrote pc123_loadings.csv, pc123_variance_explained.csv\n\n")

# ============================================================
# STEP 2 -- Lollipop plot: PC1/PC2/PC3 variance explained per measure
# (Panel A style: one row per measure, PC1=circle, PC2=triangle, PC3=square)
# ============================================================
cat("==== STEP 2: lollipop plot ====\n")
lollipop_wide <- var_exp_df %>%
  pivot_wider(names_from = PC, values_from = pct_variance) %>%
  mutate(measure = fct_reorder(measure, PC1))

p_lollipop <- ggplot(lollipop_wide) +
  geom_segment(aes(x = PC3, xend = PC1, y = measure, yend = measure), colour = "grey75") +
  geom_point(aes(x = PC1, y = measure, colour = measure), shape = 16, size = 3) +
  geom_point(aes(x = PC2, y = measure, colour = measure), shape = 17, size = 3) +
  geom_point(aes(x = PC3, y = measure, colour = measure), shape = 15, size = 3) +
  guides(colour = "none") +
  labs(x = "Variance explained (%)", y = NULL,
       title = "A  Within-measure PC1/PC2/PC3 variance explained",
       subtitle = "circle = PC1, triangle = PC2, square = PC3") +
  theme_classic(base_size = 11)

ggsave(file.path(out_dir, "panelA_lollipop.pdf"), p_lollipop, width = 6.5, height = 6)
cat("Wrote panelA_lollipop.pdf\n\n")

# ============================================================
# STEP 3 -- Cross-measure similarity of PC1/PC2/PC3 loadings
# (Spearman correlation of loadings between every pair of the 12
# measures, done separately for PC1, PC2, PC3; sign is arbitrary per
# measure/component, so magnitude/|r| is what's interpretable here)
# ============================================================
cat("==== STEP 3: cross-measure PC1/PC2/PC3 similarity ====\n")
cross_measure_similarity <- function(pc_name) {
  wide <- pc123_long %>%
    dplyr::select(parcel, measure, value = all_of(pc_name)) %>%
    pivot_wider(names_from = measure, values_from = value)
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$parcel
  cor(mat, method = "spearman", use = "pairwise.complete.obs")
}

sim_list <- setNames(lapply(c("PC1", "PC2", "PC3"), cross_measure_similarity), c("PC1", "PC2", "PC3"))
saveRDS(sim_list, file.path(out_dir, "cross_measure_similarity.rds"))

plot_sim_heatmap <- function(mat, title_str) {
  df <- as.data.frame(mat) %>%
    rownames_to_column("m1") %>%
    pivot_longer(-m1, names_to = "m2", values_to = "r")
  ggplot(df, aes(m1, m2, fill = r)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.6) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0, limits = c(-1, 1)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL, title = title_str, fill = "Spearman r")
}

for (pc in c("PC1", "PC2", "PC3")) {
  ggsave(file.path(out_dir, sprintf("panelC_%s_similarity.pdf", pc)),
         plot_sim_heatmap(sim_list[[pc]], sprintf("%s loading similarity across measures", pc)),
         width = 7, height = 6)
}
cat("Wrote panelC_PC1/2/3_similarity.pdf, cross_measure_similarity.rds\n\n")

# ============================================================
# STEP 4 -- Brain plots: PC1 as main figure, PC2/PC3 as separate
# supplementary figures (same rendering logic as Final_Figure3_script.R)
# ============================================================
cat("==== STEP 4: brain-map grids ====\n")
atlas_obj <- ggsegGlasser::glasser()
lh_sf_all <- atlas_obj$data$sf %>% dplyr::filter(grepl("^lh_L_", label))

make_pc_brain <- function(score_vec, view_name, title_str, lim = NULL) {
  sf_view <- lh_sf_all %>% dplyr::filter(view == view_name)
  sf_view$parcel <- gsub("^lh_L_(.*)_ROI$", "\\1", sf_view$label)
  sf_view$val <- score_vec[sf_view$parcel]
  if (is.null(lim)) lim <- max(abs(score_vec), na.rm = TRUE)
  bb    <- sf::st_bbox(sf_view$geometry)
  pad_x <- (bb["xmax"] - bb["xmin"]) * 0.02
  pad_y <- (bb["ymax"] - bb["ymin"]) * 0.02

  ggplot(sf_view) +
    geom_sf(aes(fill = val, geometry = geometry), colour = "#DDDDDD", linewidth = 0.05) +
    scale_fill_gradient2(low = "#3A6EA5", mid = "white", high = "#8B0000",
                          midpoint = 0, limits = c(-lim, lim), na.value = "#EEEEEE", name = NULL) +
    coord_sf(xlim = c(bb["xmin"] - pad_x, bb["xmax"] + pad_x),
             ylim = c(bb["ymin"] - pad_y, bb["ymax"] + pad_y), expand = FALSE) +
    theme_void(base_size = 7) +
    theme(legend.position = "none",
          plot.title = element_text(size = 7, hjust = 0.5),
          plot.margin = margin(1, 1, 1, 1, "pt")) +
    labs(title = title_str)
}

build_pc_grid <- function(pc_name) {
  panels <- lapply(names(pc123_res), function(ph) {
    sc <- pc123_res[[ph]]$scores[, pc_name]
    ve <- pc123_res[[ph]]$var_exp[pc_name]
    lim <- max(abs(sc), na.rm = TRUE)
    lab <- sprintf("%s  %s  (%.1f%%)", ph, pc_name, ve * 100)
    p_lat <- make_pc_brain(sc, "lateral", lab, lim = lim)
    p_med <- make_pc_brain(sc, "medial", "", lim = lim)
    plot_grid(p_lat, p_med, nrow = 1)
  })
  plot_grid(plotlist = panels, ncol = 3)
}

cat("  Building PC1 grid (main figure)...\n")
grid_PC1 <- build_pc_grid("PC1")
ggsave(file.path(out_dir, "mainFig_PC1_brains.pdf"), grid_PC1, width = 12, height = 14)

cat("  Building PC2/PC3 grids (supplementary)...\n")
grid_PC2 <- build_pc_grid("PC2")
grid_PC3 <- build_pc_grid("PC3")
ggsave(file.path(out_dir, "suppFig_PC2_brains.pdf"), grid_PC2, width = 12, height = 14)
ggsave(file.path(out_dir, "suppFig_PC3_brains.pdf"), grid_PC3, width = 12, height = 14)
cat("Wrote mainFig_PC1_brains.pdf, suppFig_PC2_brains.pdf, suppFig_PC3_brains.pdf\n\n")

# ============================================================
# STEP 5 -- Reference maps: restricted Sydnor maps + AHBA C1-C3
# ============================================================
cat("==== STEP 5: reference maps (Sydnor + AHBA) ====\n")
brainmaps <- read_csv(file.path(sydnor_dir, "brainmaps_glasser.csv"), show_col_types = FALSE)
sa_axis   <- read_csv(file.path(sydnor_dir, "Sensorimotor_Association_Axis_AverageRanks.csv"),
                       show_col_types = FALSE)
stopifnot(nrow(glasser_regions) == nrow(brainmaps))

sydnor_long <- brainmaps %>%
  mutate(parcel = glasser_regions$region) %>%
  pivot_longer(-parcel, names_to = "map_name", values_to = "value") %>%
  bind_rows(sa_axis %>% transmute(parcel = region, map_name = "SA_axis", value = final.rank))

cat("Reading AHBA C1-C3 (Hansen et al. 2024 Nat Neurosci, DME components, sheet '1a')...\n")
ahba_raw <- read_excel(ahba_path, sheet = "1a") %>% filter(!is.na(id))
stopifnot(all(c("label", "C1", "C2", "C3") %in% names(ahba_raw)))
ahba_long <- ahba_raw %>%
  transmute(parcel = label, C1, C2, C3) %>%
  pivot_longer(-parcel, names_to = "map_name", values_to = "value") %>%
  mutate(map_name = paste0("AHBA_", map_name))

n_ahba_parcels <- n_distinct(ahba_long$parcel)
cat(sprintf("AHBA C1-C3 covers %d / 180 Glasser parcels -- this is what defines the 'restricted'\n", n_ahba_parcels))
cat("parcel set used for every comparison below (see header note).\n")

restricted_parcels <- Reduce(intersect, list(
  glasser_regions_180,
  unique(pc123_long$parcel),
  unique(sydnor_long$parcel),
  unique(ahba_long$parcel)
))
cat(sprintf("Restricted parcel set (genetic PCs x Sydnor x AHBA): %d parcels\n", length(restricted_parcels)))
write_csv(tibble(parcel = restricted_parcels), file.path(out_dir, "restricted_parcel_list.csv"))

reference_maps <- bind_rows(sydnor_long, ahba_long) %>% filter(parcel %in% restricted_parcels)
write_csv(reference_maps, file.path(out_dir, "reference_maps_restricted.csv"))
cat(sprintf("Reference maps: %s\n\n", paste(unique(reference_maps$map_name), collapse = ", ")))

# ============================================================
# STEP 6 -- Arealization gene lead-SNP beta profiles
# (identical loading logic to Arealisation_beta_hitparcel_validation.R
# STEP 2: primary-locus lead SNP only, per gene x measure x parcel)
# ============================================================
cat("==== STEP 6: arealization beta profiles ====\n")
index_snps <- read_csv(areal_index_snp_path, show_col_types = FALSE)
primary_loci <- index_snps %>%
  filter(status == "ok", is_primary_locus == TRUE) %>%
  distinct(symbol, phenotype_type, locus_rank)

beta_cache <- list()
get_beta_profile <- function(pheno) {
  if (!is.null(beta_cache[[pheno]])) return(beta_cache[[pheno]])
  path <- file.path(areal_beta_dir, paste0(pheno, ".csv"))
  if (!file.exists(path)) return(NULL)
  df <- read_csv(path, show_col_types = FALSE)
  beta_cache[[pheno]] <<- df
  df
}

areal_beta_profiles <- primary_loci %>%
  pmap_dfr(function(symbol, phenotype_type, locus_rank) {
    df <- get_beta_profile(phenotype_type)
    if (is.null(df)) return(NULL)
    sub <- df %>% filter(symbol == !!symbol, locus_rank == !!locus_rank, status == "ok")
    if (nrow(sub) == 0) return(NULL)
    sub %>% transmute(symbol, phenotype_type, parcel, beta)
  })
write_csv(areal_beta_profiles, file.path(out_dir, "arealisation_beta_profiles.csv"))
cat(sprintf("Loaded beta profiles for %d gene x measure pairs\n\n",
            n_distinct(areal_beta_profiles[c("symbol", "phenotype_type")])))

cat("==== Fast steps done. Everything below this point (STEP 7-9) is the slow\n")
cat("spherical spin-permutation block -- rerun just this part once STEPS 1-6\n")
cat("look right. ====\n\n")

# ============================================================
# STEP 7 -- Spherical spin-test infrastructure (SLOW -- run last)
# ============================================================
cat("==== STEP 7: spin-test setup ====\n")
rotate_parcellation_path <- file.path(sydnor_dir, "rotate.parcellation.R")
perm_sphere_p_path       <- file.path(sydnor_dir, "perm.sphere.p.R")
sphere_path              <- file.path(sydnor_dir, "glasser_spherical_coordinates.txt")
stopifnot(file.exists(rotate_parcellation_path, perm_sphere_p_path, sphere_path))
source(rotate_parcellation_path)  # provides rotate.parcellation() -- requires library(clue)
source(perm_sphere_p_path)

sphere_full <- as.matrix(read.table(sphere_path, header = FALSE))
stopifnot(nrow(sphere_full) >= 180)  # sphere_HCP.txt: 360 rows, left hemi (1:180) then right (181:360)
sphere_lh <- sphere_full[1:180, , drop = FALSE]
rownames(sphere_lh) <- glasser_regions_180
sphere_lh_ordered <- sphere_lh[restricted_parcels, , drop = FALSE]

cat(sprintf("Generating %d spherical rotations for %d restricted parcels (this can take a\n",
            n_perm, length(restricted_parcels)))
cat("few minutes)...\n")
set.seed(seed)
# rotate.parcellation() requires both coord.l and coord.r; passing the same
# left-hemisphere coordinates as both and keeping only the first block gives
# a valid single-hemisphere spin (see PC1_spatial_correlation.R for the same
# approach / rationale).
perm_full <- rotate.parcellation(coord.l = sphere_lh_ordered, coord.r = sphere_lh_ordered,
                                  nrot = n_perm, method = "vasa")
spins <- perm_full[1:length(restricted_parcels), , drop = FALSE]
saveRDS(spins, file.path(out_dir, "spins_restricted.rds"))

spin_test <- function(x, y, spins, method = "spearman") {
  obs <- suppressWarnings(cor(x, y, method = method, use = "complete.obs"))
  null <- apply(spins, 2, function(idx) suppressWarnings(cor(x, y[idx], method = method, use = "complete.obs")))
  p_spin <- mean(abs(null) >= abs(obs), na.rm = TRUE)
  list(rho = obs, p_spin = p_spin)
}
cat("Spin-test infrastructure ready.\n\n")

# ============================================================
# STEP 8 -- PC1/PC2/PC3 vs restricted Sydnor + AHBA maps
# ============================================================
cat("==== STEP 8: PC1/2/3 vs reference maps (spin test) ====\n")
map_wide <- reference_maps %>% pivot_wider(names_from = map_name, values_from = value)
map_wide <- map_wide[match(restricted_parcels, map_wide$parcel), ]
stopifnot(identical(map_wide$parcel, restricted_parcels))
map_names <- setdiff(names(map_wide), "parcel")

pc_vs_maps_results <- expand_grid(measure = names(pc123_res), PC = c("PC1", "PC2", "PC3"),
                                   map_name = map_names) %>%
  pmap_dfr(function(measure, PC, map_name) {
    pc_sub  <- pc123_long %>% filter(measure == !!measure)
    pc_vec  <- pc_sub[[PC]][match(restricted_parcels, pc_sub$parcel)]
    map_vec <- map_wide[[map_name]]
    ok <- !is.na(pc_vec) & !is.na(map_vec)
    if (sum(ok) < 20) {
      return(tibble(measure, PC, map_name, rho = NA_real_, p_spin = NA_real_, n = sum(ok)))
    }
    r <- spin_test(pc_vec, map_vec, spins)
    tibble(measure, PC, map_name, rho = r$rho, p_spin = r$p_spin, n = sum(ok))
  }) %>%
  group_by(map_name, PC) %>%
  mutate(p_spin_fdr = p.adjust(p_spin, method = "BH")) %>%
  ungroup() %>%
  arrange(map_name, PC, desc(abs(rho)))

write_csv(pc_vs_maps_results, file.path(out_dir, "PC123_vs_reference_maps_spintest.csv"))
cat("Wrote PC123_vs_reference_maps_spintest.csv\n\n")

# ============================================================
# STEP 9 -- Arealisation beta profiles vs PC1/2/3 and reference maps
# (|beta| used throughout, matching Arealisation_beta_hitparcel_
# validation.R's Test A/B convention -- this tests whether effect
# MAGNITUDE aligns spatially with each map, not effect direction,
# since PC/Sydnor/AHBA sign is not a fixed, comparable reference)
# ============================================================
cat("==== STEP 9: arealization betas vs PC1/2/3 + reference maps (spin test) ====\n")
gene_measure_pairs <- areal_beta_profiles %>% distinct(symbol, phenotype_type)
cat(sprintf("Running spin tests for %d gene x measure pairs x (3 PCs + %d reference maps)...\n",
            nrow(gene_measure_pairs), length(map_names)))

areal_vs_pc_results <- gene_measure_pairs %>%
  pmap_dfr(function(symbol, phenotype_type) {
    beta_sub <- areal_beta_profiles %>% filter(symbol == !!symbol, phenotype_type == !!phenotype_type)
    beta_vec <- beta_sub$beta[match(restricted_parcels, beta_sub$parcel)]
    if (sum(!is.na(beta_vec)) < 20) return(NULL)
    abs_beta <- abs(beta_vec)

    pc_sub <- pc123_long %>% filter(measure == phenotype_type)
    out_pc <- map_dfr(c("PC1", "PC2", "PC3"), function(PC) {
      pc_vec <- pc_sub[[PC]][match(restricted_parcels, pc_sub$parcel)]
      ok <- !is.na(pc_vec) & !is.na(abs_beta)
      if (sum(ok) < 20) return(tibble(ref = "own_measure_PC", ref_name = PC, rho = NA_real_, p_spin = NA_real_, n = sum(ok)))
      r <- spin_test(abs_beta, pc_vec, spins)
      tibble(ref = "own_measure_PC", ref_name = PC, rho = r$rho, p_spin = r$p_spin, n = sum(ok))
    })

    out_maps <- map_dfr(map_names, function(mn) {
      map_vec <- map_wide[[mn]]
      ok <- !is.na(map_vec) & !is.na(abs_beta)
      if (sum(ok) < 20) return(tibble(ref = "reference_map", ref_name = mn, rho = NA_real_, p_spin = NA_real_, n = sum(ok)))
      r <- spin_test(abs_beta, map_vec, spins)
      tibble(ref = "reference_map", ref_name = mn, rho = r$rho, p_spin = r$p_spin, n = sum(ok))
    })

    bind_rows(out_pc, out_maps) %>% mutate(symbol = symbol, phenotype_type = phenotype_type, .before = 1)
  }) %>%
  group_by(ref_name) %>%
  mutate(p_spin_fdr = p.adjust(p_spin, method = "BH")) %>%
  ungroup() %>%
  arrange(ref_name, desc(abs(rho)))

write_csv(areal_vs_pc_results, file.path(out_dir, "arealisation_beta_vs_PC123_and_maps_spintest.csv"))
cat("Wrote arealisation_beta_vs_PC123_and_maps_spintest.csv\n\n")

cat(sprintf("All done. Outputs in %s/\n", out_dir))
