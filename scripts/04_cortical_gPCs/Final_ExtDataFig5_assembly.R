# ===========================================================================
# Final_ExtDataFig5_assembly.R
# Manuscript item: Figure 4 and Supplementary Figure 8
#
# What this does:
#   Assembles the gPC figure: variance explained, gPC1 maps by measure
#   (signs aligned for display only), within- vs between-family |Pearson r|
#   violins of the loadings (panel C), and the gPC2/gPC3 supplementary brain
#   maps and violins.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Assembles the final Extended Data Figure 5 (genetic topographical
# organisation of the cortex) and its associated supplementary
# figures.
#
#   Fig 5A -- PC1-3 variance explained (lollipop, with legend)
#   Fig 5B -- PC1 spatial projection on the Glasser parcellation, by
#             measure, arranged 4-across x 3-down (was 3x4 - flattened
#             to cut vertical whitespace, see STEP 4/7). Each measure's
#             PC1 sign is realigned (see align_signs() below) so that,
#             as far as a single global sign flip per measure can
#             achieve it, PC1 maps are POSITIVELY correlated with each
#             other across measures -- purely a visualisation choice;
#             sign is otherwise arbitrary and doesn't affect any
#             reported statistic. A single shared colour key (red =
#             positive loading, blue = negative, per-measure scaled to
#             its own range since limits differ by measure) sits below
#             the grid.
#   Fig 5C -- violin plot of |Pearson r| between measures' PC1-3
#             loadings, split into within-family vs between-family
#             pairs (family = the paper's 6 canonical families: brain
#             size, cortical thickness, curvature, microstructural
#             coherence, microstructural diffusivity, orientation
#             dispersion -- see family_map below). Compressed
#             (narrower, smaller legend/fonts) so it sits comfortably
#             in the top row next to A and D -- see STEP 5.
#   Fig 5D -- PC1-3 vs a restricted set of 7 reference maps (S-A axis,
#             allometric scaling, evolutionary expansion, T1w/T2w
#             myelin, AHBA C1-3 [Dear et al.]), FDR-corrected within
#             that restricted set (see ASSUMPTION note below). Only 3
#             of those 7 (S-A axis, allometric scaling, AHBA C1) are
#             actually plotted, as |rho| fill (now RED scale, was
#             blue -- see STEP 6) with a bare "*" for FDR<0.05 -- no
#             numbers shown, per instruction. Moved into the top row
#             alongside A and C (previously its own row) -- see STEP 7.
#   Supp Fig 6 -- PC2 and PC3 spatial projections, by measure (same
#             per-PC sign realignment and 4x3 layout as Fig 5B).
#   Supp Fig 7 -- same within/between-family |r| violin as Fig 5C,
#             for PC2 and PC3.
#
# ASSUMPTION (flagged in chat): FDR correction for Fig 5D is done
# within each gPC, across all 7 restricted reference maps x 12
# measures together (not per-map, and not restricted to just the 3
# AHBA axes). Change the `group_by()` call in STEP 6 if you meant
# something narrower.
#
# CHANGELOG (2026-09-01, per feedback on the first assembled draft):
#   - Top row is now A + C + D together (was A+C, then B, then D as
#     three separate stacked rows) -- see STEP 7.
#   - Panel C narrowed and its legend/fonts shrunk so it fits
#     comfortably in that row -- see STEP 5.
#   - Panel B's brain grid changed from 3-across x 4-down to
#     4-across x 3-down, and inter-row spacing tightened, to cut the
#     large vertical whitespace it had at full page height -- STEP 4.
#   - Panel B gained a shared colour-direction legend (it had none
#     before -- STEP 4).
#   - Panel D's fill scale changed from white->blue to white->red --
#     STEP 6.
#
# This script does NOT rerun the spherical spin-permutation test --
# it reads the already-computed PC123_vs_reference_maps_spintest.csv
# from Final_ExtDataFig4_PCA_master.R (10,000 rotations, already run)
# and just refilters/re-corrects on the restricted map set. Everything
# else (PC1-3 loadings, lollipop, family violins, brain maps) is
# fast/deterministic and computed fresh here.
#
# Inputs:
#   ~/Downloads/rstudio-export (3)/rg_matrix.rds
#   Data/Sydnor_maps_scripts/glasser_regions.csv
#   Outputs/ExtDataFig4_PCA_master/PC123_vs_reference_maps_spintest.csv
#
# Outputs (all under out_dir):
#   ExtDataFig5_main.pdf          -- Panels A, C, D (top row) + B (below)
#   SuppFig6_PC2_PC3_brains.pdf   -- PC2 (top) and PC3 (bottom) brain grids
#   SuppFig7_PC2_PC3_family_violin.pdf -- within/between-family |r| violins, PC2 & PC3

source("config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsegGlasser)
  library(cowplot)
  library(sf)
})

# ---------------- USER INPUTS ----------------
rg_path       <- file.path(proj_dir, "Data", "rg_matrix.rds")
sydnor_dir    <- "Data/Sydnor_maps_scripts"
spintest_path <- "Outputs/ExtDataFig4_PCA_master/PC123_vs_reference_maps_spintest.csv"
out_dir       <- "Outputs/ExtDataFig5_assembly"
all_phenos    <- c("CT","FA","FI","GC","GMV","IC","ICVF","ISOVF","MC","MD","OD","SA")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(spintest_path)) {
  stop("Missing ", spintest_path, " -- run Final_ExtDataFig4_PCA_master.R (through STEP 8) first.")
}

glasser_regions     <- read_csv(file.path(sydnor_dir, "glasser_regions.csv"), show_col_types = FALSE)
glasser_regions_180 <- glasser_regions$region

# The paper's 6 canonical families (abstract): brain size, cortical
# thickness, curvature, microstructural coherence, microstructural
# diffusivity, orientation dispersion. FI/IC assigned to curvature
# (not brain size) per the already-established cross-measure
# similarity result (|r| up to 0.90 with GC, 0.82 with IC).
family_map <- c(SA = "BrainSize", GMV = "BrainSize",
                 GC = "Curvature", MC = "Curvature", FI = "Curvature", IC = "Curvature",
                 CT = "Thickness",
                 FA = "Coherence", ICVF = "Coherence",
                 MD = "Diffusivity", ISOVF = "Diffusivity",
                 OD = "OrientationDispersion")

# ============================================================
# STEP 1 -- PC1-3 loadings (fast; identical construction to
# Final_ExtDataFig4_PCA_master.R STEP 1)
# ============================================================
cat("==== STEP 1: PC1-3 loadings ====\n")
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

pc123_res <- Filter(Negate(is.null), setNames(lapply(all_phenos, pheno_pca_pc123), all_phenos))

pc123_long <- map_dfr(names(pc123_res), function(ph) {
  sc <- pc123_res[[ph]]$scores
  tibble(parcel = rownames(sc), measure = ph, PC1 = sc[, "PC1"], PC2 = sc[, "PC2"], PC3 = sc[, "PC3"])
})

var_exp_df <- map_dfr(names(pc123_res), function(ph) {
  ve <- pc123_res[[ph]]$var_exp
  tibble(measure = ph, PC = paste0("PC", 1:3), pct_variance = as.numeric(ve[1:3]) * 100)
})

pc1_order <- var_exp_df %>% filter(PC == "PC1") %>% arrange(pct_variance) %>% pull(measure)
var_exp_df <- var_exp_df %>% mutate(measure = factor(measure, levels = pc1_order))
cat("Done.\n\n")

# ============================================================
# STEP 2 -- Panel A: lollipop, with a proper shape+colour legend
# ============================================================
cat("==== STEP 2: Panel A ====\n")
seg_df <- var_exp_df %>%
  group_by(measure) %>%
  summarise(xmin = min(pct_variance), xmax = max(pct_variance), .groups = "drop")

pc_colours <- c(PC1 = "#0C447C", PC2 = "#993C1D", PC3 = "#085041")
pc_shapes  <- c(PC1 = 16, PC2 = 17, PC3 = 15)

panel_A <- ggplot() +
  geom_segment(data = seg_df, aes(x = xmin, xend = xmax, y = measure, yend = measure), colour = "grey75") +
  geom_point(data = var_exp_df, aes(x = pct_variance, y = measure, shape = PC, colour = PC), size = 2.6) +
  scale_shape_manual(values = pc_shapes, name = "Component") +
  scale_colour_manual(values = pc_colours, name = "Component") +
  labs(x = "Variance explained (%)", y = NULL, title = "A") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7.5))

# ============================================================
# STEP 3 -- Sign alignment for brain-map visualisation only.
# For each PC, greedily flip each measure's sign to maximise total
# pairwise agreement (sum of signed correlations) with the other
# measures' loadings -- a purely cosmetic step so that, e.g., "red"
# means the same direction across the panel-B grid as far as a single
# global sign per measure allows. Does not touch pc123_long/var_exp_df
# (used for Panels A, C, D), which keep the original mean<0-flip
# convention already reported in the text.
# ============================================================
cat("==== STEP 3: sign alignment for brain maps ====\n")
align_signs <- function(pc_name) {
  wide <- pc123_long %>%
    dplyr::select(parcel, measure, value = all_of(pc_name)) %>%
    pivot_wider(names_from = measure, values_from = value)
  mat <- as.matrix(wide[, -1]); rownames(mat) <- wide$parcel
  measures <- colnames(mat); n <- length(measures)
  cm <- matrix(NA_real_, n, n, dimnames = list(measures, measures))
  for (i in 1:n) for (j in 1:n) {
    ok <- stats::complete.cases(mat[, c(i, j)])
    cm[i, j] <- if (i == j) 1 else suppressWarnings(cor(mat[ok, i], mat[ok, j]))
  }
  cm[is.na(cm)] <- 0
  signs <- setNames(rep(1, n), measures)
  for (iter in 1:50) {
    improved <- FALSE
    for (i in 1:n) {
      cur <- sum(signs[i] * signs * cm[i, ])
      alt <- signs; alt[i] <- -alt[i]
      new <- sum(alt[i] * alt * cm[i, ])
      if (new > cur + 1e-9) { signs[i] <- -signs[i]; improved <- TRUE }
    }
    if (!improved) break
  }
  signs
}

signs_PC1 <- align_signs("PC1")
signs_PC2 <- align_signs("PC2")
signs_PC3 <- align_signs("PC3")
cat("Sign flips applied (PC1):", paste(names(signs_PC1)[signs_PC1 < 0], collapse = ", "), "\n")
cat("Sign flips applied (PC2):", paste(names(signs_PC2)[signs_PC2 < 0], collapse = ", "), "\n")
cat("Sign flips applied (PC3):", paste(names(signs_PC3)[signs_PC3 < 0], collapse = ", "), "\n\n")

# ============================================================
# STEP 4 -- Panel B (PC1) and Supp Fig 6 (PC2, PC3) brain grids,
# using the aligned signs from STEP 3.
# Grid is now 4-across x 3-down (was 3-across x 4-down) to cut the
# vertical whitespace the taller layout had, and a single shared
# colour-direction legend is appended below the grid: since each
# measure's fill scale is independently limited to its own max |PC|
# (magnitudes aren't comparable across measures), the legend shows
# direction only (red = positive loading, blue = negative), not a
# shared numeric scale.
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
  # Pin the PANEL (not just the data) to the brain's true height/width
  # ratio via aspect.ratio. Without this, coord_sf's forced true-aspect
  # rendering still lets the surrounding panel area be any shape, so
  # cowplot's equal-height row packing leaves the actual brain sitting
  # inside a taller-than-needed box -- visually reads as "space between
  # rows" even though it's really unused vertical space WITHIN each
  # panel. Pinning aspect.ratio here (combined with the tighter overall
  # ggsave height below) removes that slack.
  asp <- as.numeric((bb["ymax"] - bb["ymin"]) / (bb["xmax"] - bb["xmin"]))
  ggplot(sf_view) +
    geom_sf(aes(fill = val, geometry = geometry), colour = "#DDDDDD", linewidth = 0.05) +
    scale_fill_gradient2(low = "#3A6EA5", mid = "white", high = "#8B0000",
                          midpoint = 0, limits = c(-lim, lim), na.value = "#EEEEEE", name = NULL) +
    coord_sf(xlim = c(bb["xmin"] - pad_x, bb["xmax"] + pad_x),
             ylim = c(bb["ymin"] - pad_y, bb["ymax"] + pad_y), expand = FALSE) +
    theme_void(base_size = 7) +
    theme(legend.position = "none",
          aspect.ratio = asp,
          plot.title = element_text(size = 7, hjust = 0.5, margin = margin(b = 1)),
          plot.margin = margin(0, 0.5, 0, 0.5, "pt")) +
    labs(title = title_str)
}

# Small standalone legend: a continuous red/white/blue bar labelled
# "+" / "0" / "-" only (no numbers, since scale differs per measure).
brain_direction_legend <- function() {
  df <- tibble(x = seq(-1, 1, length.out = 100), y = 1)
  p <- ggplot(df, aes(x, y, fill = x)) +
    geom_tile(height = 1) +
    scale_fill_gradient2(low = "#3A6EA5", mid = "white", high = "#8B0000", midpoint = 0) +
    scale_x_continuous(breaks = c(-1, 0, 1), labels = c("-", "0", "+")) +
    theme_void(base_size = 8) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8),
          axis.ticks.x = element_line(linewidth = 0.2),
          axis.line.x = element_line(linewidth = 0.2)) +
    labs(x = NULL, y = NULL)
  cowplot::plot_grid(
    ggdraw() + draw_label("PC loading (sign-aligned, own-measure scale):", size = 7.5, hjust = 0, x = 0.02),
    p, ncol = 1, rel_heights = c(0.4, 0.6)
  )
}

build_pc_grid <- function(pc_name, signs, ncol_grid = 4, with_legend = FALSE) {
  panels <- lapply(names(pc123_res), function(ph) {
    sc <- pc123_res[[ph]]$scores[, pc_name] * signs[[ph]]
    ve <- pc123_res[[ph]]$var_exp[pc_name]
    lim <- max(abs(sc), na.rm = TRUE)
    lab <- sprintf("%s  %s  (%.1f%%)", ph, pc_name, ve * 100)
    p_lat <- make_pc_brain(sc, "lateral", lab, lim = lim)
    p_med <- make_pc_brain(sc, "medial", "", lim = lim)
    plot_grid(p_lat, p_med, nrow = 1)
  })
  grid <- plot_grid(plotlist = panels, ncol = ncol_grid)
  if (!with_legend) return(grid)
  plot_grid(grid, brain_direction_legend(), ncol = 1, rel_heights = c(1, 0.06))
}

cat("  PC1 grid (Panel B, 4x3, with legend)...\n")
panel_B <- build_pc_grid("PC1", signs_PC1, ncol_grid = 4, with_legend = TRUE)

cat("  PC2 / PC3 grids (Supp Fig 6, 4x3)...\n")
grid_PC2 <- build_pc_grid("PC2", signs_PC2, ncol_grid = 4)
grid_PC3 <- build_pc_grid("PC3", signs_PC3, ncol_grid = 4)
supp_fig6 <- plot_grid(
  ggdraw() + draw_label("PC2", fontface = "bold", size = 12, hjust = 0),
  grid_PC2,
  ggdraw() + draw_label("PC3", fontface = "bold", size = 12, hjust = 0),
  grid_PC3,
  ncol = 1, rel_heights = c(0.04, 1, 0.04, 1)
)
ggsave(file.path(out_dir, "SuppFig6_PC2_PC3_brains.pdf"), supp_fig6, width = 14, height = 20, limitsize = FALSE)
cat("Done.\n\n")

# ============================================================
# STEP 5 -- Panel C: within/between-family |r| violin (PC1-3), and
# Supp Fig 7: the same for PC2/PC3 only.
# Panel C is now built COMPACT (narrow, small fonts/legend) so it fits
# in the top row alongside A and D -- see STEP 7. Supp Fig 7 keeps the
# original, larger sizing since it's a standalone figure.
# ============================================================
cat("==== STEP 5: family |r| violin plots ====\n")
build_pairs <- function(pc_name) {
  wide <- pc123_long %>%
    dplyr::select(parcel, measure, value = all_of(pc_name)) %>%
    pivot_wider(names_from = measure, values_from = value)
  mat <- as.matrix(wide[, -1]); rownames(mat) <- wide$parcel
  cm <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  measures <- colnames(cm)
  combn(measures, 2) %>% t() %>% as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("m1", "m2")) %>%
    as_tibble() %>%
    mutate(abs_r = abs(cm[cbind(m1, m2)]),
           relation = ifelse(family_map[m1] == family_map[m2], "Within family", "Between family"),
           PC = pc_name)
}
pairs_all <- bind_rows(lapply(c("PC1", "PC2", "PC3"), build_pairs)) %>%
  mutate(PC = factor(PC, levels = c("PC1", "PC2", "PC3")))

family_violin <- function(df, title_str, compact = FALSE) {
  p <- ggplot(df, aes(x = PC, y = abs_r, fill = relation)) +
    geom_violin(position = position_dodge(0.8), width = 0.7, colour = NA, alpha = 0.6) +
    geom_boxplot(position = position_dodge(0.8), width = 0.12, outlier.size = 0.5) +
    scale_fill_manual(values = c("Within family" = "#0F6E56", "Between family" = "#993C1D"), name = NULL) +
    labs(x = NULL, y = "|Pearson r|", title = title_str) +
    theme_classic(base_size = if (compact) 7.5 else 10)
  if (compact) {
    p + theme(legend.position = "bottom", legend.direction = "horizontal",
              legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 6.5),
              plot.title = element_text(face = "bold", size = 11),
              axis.text = element_text(size = 6.5), axis.title = element_text(size = 7.5),
              plot.margin = margin(2, 2, 2, 2, "pt"))
  } else {
    p + theme(legend.position = "top")
  }
}

panel_C <- family_violin(pairs_all, "C", compact = TRUE)
supp_fig7 <- family_violin(pairs_all %>% filter(PC %in% c("PC2", "PC3")) %>% mutate(PC = droplevels(PC)), NULL)
ggsave(file.path(out_dir, "SuppFig7_PC2_PC3_family_violin.pdf"), supp_fig7, width = 6, height = 5)
cat("Done.\n\n")

# ============================================================
# STEP 6 -- Panel D: PC1-3 vs restricted reference maps.
# FDR-corrected across the 7-map restricted set within each PC (see
# ASSUMPTION note in the header); only 3 of those 7 maps are plotted.
# Reuses the already-computed spin-test p-values -- no permutation
# rerun. Fill scale changed white->RED (was white->blue) per feedback;
# compact sizing (smaller fonts/tiles) so it fits in the top row.
# ============================================================
cat("==== STEP 6: Panel D ====\n")
spin_raw <- read_csv(spintest_path, show_col_types = FALSE)

fdr_maps     <- c("SA_axis", "AllometricScaling.PNC20mm", "Evolution.Expansion", "T1T2ratio",
                   "AHBA_C1", "AHBA_C2", "AHBA_C3")
display_maps <- c("SA_axis", "AllometricScaling.PNC20mm", "AHBA_C1")
display_labels <- c(SA_axis = "S-A axis", AllometricScaling.PNC20mm = "Allometric scaling",
                     AHBA_C1 = "AHBA C1 (Dear et al.)")

spin_restricted <- spin_raw %>%
  filter(map_name %in% fdr_maps) %>%
  group_by(PC) %>%                       # ASSUMPTION: FDR within each PC across all 7 maps together
  mutate(p_spin_fdr = p.adjust(p_spin, method = "BH")) %>%
  ungroup()

spin_display <- spin_restricted %>%
  filter(map_name %in% display_maps) %>%
  mutate(abs_rho = abs(rho),
         sig = !is.na(p_spin_fdr) & p_spin_fdr < 0.05,
         map_name = factor(map_name, levels = display_maps, labels = display_labels[display_maps]),
         measure  = factor(measure, levels = as.character(pc1_order)),
         PC       = factor(PC, levels = c("PC1", "PC2", "PC3")))

cat(sprintf("Panel D (display set): %d cells, %d FDR-significant\n",
            nrow(spin_display), sum(spin_display$sig, na.rm = TRUE)))

panel_D <- ggplot(spin_display, aes(x = map_name, y = measure, fill = abs_rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(sig, "*", "")), size = 3, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, 1), name = "|rho|") +
  facet_wrap(~PC, ncol = 3) +
  theme_minimal(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5),
        panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, title = "D",
       subtitle = "* FDR<0.05 (restricted 7-map set, spin-test permutation)")
cat("Done.\n\n")

# ============================================================
# STEP 7 -- Assemble Fig 5: top row = A + C + D together, Panel B
# (brain grid) spans the full width below. (Previously: A+C on one
# row, B on its own row, D on its own row -- three stacked rows before
# reaching B; now compressed to two.)
# ============================================================
cat("==== STEP 7: assembling Fig 5 ====\n")
# A narrowed further (0.65, was 1); D widened slightly to compensate
# since it carries 3 facets. Panel B given less relative height (1.6,
# was 2.3) and the whole figure's saved height cut (11, was 15) --
# together with the per-panel aspect.ratio fix in STEP 4, this is what
# actually compresses the brain grid rather than just reallocating
# proportions within the same oversized canvas.
top_row <- plot_grid(panel_A, panel_C, panel_D, ncol = 3, rel_widths = c(0.65, 0.8, 1.45))
fig5_main <- plot_grid(top_row, panel_B, ncol = 1, rel_heights = c(1, 1.6))

ggsave(file.path(out_dir, "ExtDataFig5_main.pdf"), fig5_main, width = 14, height = 11, limitsize = FALSE)

cat(sprintf("\nDone. Outputs in %s/:\n", out_dir))
cat("  ExtDataFig5_main.pdf                -- Panels A, C, D (top row) + B (below)\n")
cat("  SuppFig6_PC2_PC3_brains.pdf         -- PC2 and PC3 brain maps (sign-aligned, 4x3)\n")
cat("  SuppFig7_PC2_PC3_family_violin.pdf  -- within/between-family |r| violins, PC2 & PC3\n")
