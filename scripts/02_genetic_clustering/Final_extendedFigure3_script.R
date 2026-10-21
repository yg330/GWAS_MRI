# ===========================================================================
# Final_extendedFigure3_script.R
# Manuscript item: Supplementary Figure 4A and Extended Data Figure 4B-C
#
# What this does:
#   Fine-scale structure of the genetic correlation matrix: (A) within-
#   measure rg densities for the 12 cortical measures, (B) rg heatmap of
#   subcortical volumes / corpus callosum / ventricles, (C) how isolated
#   each tract is from the other tracts, per microstructural measure.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Extended Data Figure 3: Fine-scale genetic architecture
#
# Panel A: Faceted kernel density plots of within-measure
#          genetic correlations for all 12 cortical measures
# Panel B: Genetic correlation heatmap of subcortical volumes,
#          corpus callosum, and ventricles
# Panel C: Genetic isolation plot - mean rg of each tract to
#          all other tracts, shown separately for all 5
#          microstructural measures (Cleveland dot plot)
#
# Input:  rg_matrix.rds - full pairwise LDSC genetic correlation matrix
# Output: ExtDataFig3_A.pdf/.png, ExtDataFig3_B.pdf/.png,
#         ExtDataFig3_C.pdf/.png, ExtDataFig3_combined.pdf/.png

source("config.R")

library(tidyverse)
library(igraph)
library(ggplot2)
library(cowplot)
library(ggdendro)
library(reshape2)
library(scales)

# ---------------------
# Global typography (shared with Figure 2)
# ---------------------
BASE_SIZE   <- 10
TITLE_SIZE  <- 10
LABEL_SIZE  <- 9
TILE_SIZE   <- 2.8
BASE_FAMILY <- ""

theme_fig <- function() {
  theme_classic(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    theme(
      plot.title      = element_text(size = TITLE_SIZE, face = "bold",
                                     hjust = 0, margin = margin(b = 4)),
      axis.text       = element_text(size = LABEL_SIZE),
      axis.title      = element_text(size = LABEL_SIZE),
      legend.text     = element_text(size = LABEL_SIZE),
      legend.title    = element_text(size = LABEL_SIZE),
      legend.key.size = unit(0.35, "cm"),
      strip.text      = element_text(size = LABEL_SIZE)
    )
}

# ---------------------
# Colour palettes (shared with Figure 2)
# ---------------------
colors_measure <- c(
  CT             = "#0072B2",
  FA             = "#C0006E",
  FI             = "#7DDCB8",
  GC             = "#6B4800",
  GMV            = "#003D28",
  IC             = "#2DC48E",
  ICVF           = "#E8A8CB",
  ISOVF          = "#A8D8F0",
  MC             = "#F5C842",
  MD             = "#0D3F5F",
  OD             = "#D55E00",
  SA             = "#007A52",
  CorpusCallosum = "#8B0000",
  Ventricle      = "#4169E1",
  SubCorVol      = "#9370DB"
)

colors_micro <- c(
  FA    = "#C0006E",
  ICVF  = "#E8A8CB",
  ISOVF = "#A8D8F0",
  MD    = "#0D3F5F",
  OD    = "#D55E00"
)

colors_subcor_group <- c(
  "Ventricles"             = "#4169E1",
  "Corpus callosum"        = "#8B0000",
  "Striatal/limbic"        = "#6A3D9A",
  "Hindbrain/diencephalic" = "#9370DB"
)

# ---------------------
# Load and clean data
# ---------------------
rg_matrix_raw <- readRDS(file.path(proj_dir, "Data", "rg_matrix.rds"))

clean_names <- function(x) gsub("_as_pheno.*", "", gsub("GWAS_", "", x))
rownames(rg_matrix_raw) <- clean_names(rownames(rg_matrix_raw))
colnames(rg_matrix_raw) <- clean_names(colnames(rg_matrix_raw))

# Exclude non-brain subcortical phenotypes (shared with Figure 2)
exclude_patterns <- c(
  "SubCorVol_choroid.plexus", "SubCorVol_CSF", "SubCorVol_Optic.Chiasm",
  "SubCorVol_SubCortGrayVol", "SubCorVol_vessel", "SubCorVol_WM.hypointensities"
)
keep_idx    <- !rownames(rg_matrix_raw) %in% exclude_patterns
rg_matrix   <- rg_matrix_raw[keep_idx, keep_idx]
trait_names <- colnames(rg_matrix)

# ---------------------
# Trait classification helper
# ---------------------
get_measure_type <- function(name) {
  if (grepl("^SubCorVol_CC_", name))                 return("CorpusCallosum")
  if (grepl("^SubCorVol_.*(Ventricle|Vent)$", name)) return("Ventricle")
  parts <- strsplit(name, "_")[[1]]
  first <- parts[1]
  if (first %in% c("global", "regional", "tract") && length(parts) >= 2) {
    if (parts[2] == "GMV" || first == "SubCorVol") return("GMV") else return(parts[2])
  } else if (first == "SubCorVol") return("SubCorVol")
  return(first)
}

measure_types <- sapply(trait_names, get_measure_type)

# ============================================================
# PANEL A: Faceted density plots - within-measure rg
#          for all 12 cortical measures
# ============================================================

cortical_measures <- c("CT", "FA", "FI", "GC", "IC", "ICVF",
                       "ISOVF", "MC", "MD", "OD", "SA", "GMV")

# Order facets by phenotype family
measure_display_order <- c("SA", "GMV", "FI", "IC",
                           "CT",
                           "GC", "MC",
                           "FA", "ICVF",
                           "MD", "ISOVF",
                           "OD")

density_rows <- list()

for (m in cortical_measures) {
  if (m == "GMV") {
    idx <- which(grepl("^regional_GMV_", trait_names))
  } else {
    idx <- which(grepl(paste0("^regional_", m, "_"), trait_names))
  }
  if (length(idx) < 2) next
  
  sub_mat <- rg_matrix[idx, idx]
  vals    <- sub_mat[upper.tri(sub_mat)]
  vals    <- vals[!is.na(vals)]
  
  density_rows[[m]] <- data.frame(
    rg      = vals,
    measure = m
  )
}

density_df <- bind_rows(density_rows) %>%
  mutate(measure = factor(measure, levels = measure_display_order))

mean_rg_df <- density_df %>%
  group_by(measure) %>%
  summarise(mean_rg = mean(rg, na.rm = TRUE), .groups = "drop")

panel_A <- ggplot(density_df, aes(x = rg, fill = measure, colour = measure)) +
  geom_density(alpha = 0.6, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_vline(data = mean_rg_df,
             aes(xintercept = mean_rg, colour = measure),
             linetype = "solid", linewidth = 0.6, show.legend = FALSE) +
  scale_fill_manual(values   = colors_measure, guide = "none") +
  scale_colour_manual(values = colors_measure, guide = "none") +
  scale_x_continuous(limits = c(-0.6, 1.05), breaks = c(-0.5, 0, 0.5, 1)) +
  facet_wrap(~ measure, ncol = 3, strip.position = "top",
             scales = "free_y") +
  theme_fig() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = LABEL_SIZE, face = "bold"),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.spacing    = unit(0.2, "cm")
  ) +
  labs(
    title = "A",
    x     = "Genetic correlation (rg)",
    y     = "Density"
  )

# ============================================================
# PANEL B: Subcortical heatmap
#          Ventricles / CC / Striatal-limbic / Hindbrain-diencephalic
# ============================================================

get_subcor_group <- function(name) {
  if (!grepl("^SubCorVol_", name))   return(NA_character_)
  if (grepl("^SubCorVol_CC_", name)) return("Corpus callosum")
  if (grepl("^SubCorVol_.*(Ventricle|Vent)$", name, ignore.case = TRUE))
    return("Ventricles")
  if (grepl("Caudate|Accumbens|Putamen|Amygdala|Hippocampus",
            name, ignore.case = TRUE))
    return("Striatal/limbic")
  if (grepl("Thalamus|Pallidum|VentralDC|Brain.Stem|Cerebellum",
            name, ignore.case = TRUE))
    return("Hindbrain/diencephalic")
  return(NA_character_)
}

shorten_subcor <- function(x) {
  x <- gsub("^SubCorVol_CC_", "CC ", x)
  x <- gsub("^SubCorVol_", "", x)
  x <- gsub("\\.", " ", x)
  x <- gsub("Cerebellum White Matter", "Cereb. WM", x)
  x <- gsub("Cerebellum Cortex", "Cereb. Ctx", x)
  x <- gsub("Brain Stem", "Brainstem", x)
  x <- gsub("VentralDC", "Ventral DC", x)
  x <- gsub("Thalamus Proper", "Thalamus", x)
  x <- gsub("Accumbens area", "Accumbens", x)
  x <- gsub("Lateral Ventricle", "Lat. Ventricle", x)
  x <- gsub("Inf Lat Vent", "Inf. Lat. Vent.", x)
  x <- gsub("X3rd Ventricle", "3rd Ventricle", x)
  x <- gsub("X4th Ventricle", "4th Ventricle", x)
  x
}

subcor_idx    <- grepl("^SubCorVol_", trait_names)
subcor_names  <- trait_names[subcor_idx]
subcor_groups <- sapply(subcor_names, get_subcor_group)

valid         <- !is.na(subcor_groups)
subcor_names  <- subcor_names[valid]
subcor_groups <- subcor_groups[valid]

group_order_B <- c("Ventricles", "Corpus callosum",
                   "Striatal/limbic", "Hindbrain/diencephalic")

subcor_order <- unlist(lapply(group_order_B, function(g) {
  nms <- subcor_names[subcor_groups == g]
  if (length(nms) < 2) return(nms)
  sub <- rg_matrix[nms, nms]
  hc  <- hclust(as.dist(1 - sub), method = "average")
  nms[hc$order]
}))

rg_subcor   <- rg_matrix[subcor_order, subcor_order]
short_names <- shorten_subcor(subcor_order)
rownames(rg_subcor) <- short_names
colnames(rg_subcor) <- short_names

melt_subcor <- melt(rg_subcor, varnames = c("Var1", "Var2"),
                    value.name = "rg") %>%
  mutate(
    Var1 = factor(Var1, levels = rev(short_names)),
    Var2 = factor(Var2, levels = short_names)
  )

melt_subcor <- melt_subcor %>%
  mutate(label = ifelse(!is.na(rg), sprintf("%.2f", rg), ""))

panel_B <- ggplot(melt_subcor, aes(x = Var2, y = Var1, fill = rg)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = TILE_SIZE - 0.8,
            colour = "black", family = BASE_FAMILY) +
  scale_fill_gradient2(
    low      = "#4575B4",
    mid      = "white",
    high     = "#D73027",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "rg"
  ) +
  theme_fig() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = LABEL_SIZE - 1),
    axis.text.y  = element_text(size = LABEL_SIZE - 1),
    axis.title   = element_blank(),
    panel.border = element_blank(),
    panel.grid   = element_blank()
  ) +
  labs(title = "B")

# ============================================================
# PANEL C: Genetic isolation plot
#          Mean rg of each tract to all other tracts,
#          shown separately for all 5 microstructural measures
#          Tracts ordered by median mean_rg (lowest = most isolated)
# ============================================================

micro_measures <- c("FA", "ICVF", "ISOVF", "MD", "OD")

isolation_rows <- list()

for (m in micro_measures) {
  tract_idx <- which(grepl(paste0("^tract_", m, "_"), trait_names) |
                       grepl(paste0("^tr_", m, "_"),    trait_names))
  if (length(tract_idx) < 2) next
  
  tract_nms <- trait_names[tract_idx]
  sub_mat   <- rg_matrix[tract_nms, tract_nms]
  diag(sub_mat) <- NA
  mean_to_others <- rowMeans(sub_mat, na.rm = TRUE)
  
  short_nms <- gsub(paste0("^tract_", m, "_|^tr_", m, "_"),
                    "", tract_nms, ignore.case = TRUE)
  short_nms <- gsub("_", " ", short_nms)
  short_nms <- toupper(short_nms)
  
  isolation_rows[[m]] <- data.frame(
    tract      = short_nms,
    tract_full = tract_nms,
    mean_rg    = mean_to_others,
    measure    = m,
    stringsAsFactors = FALSE
  )
}

isolation_df <- bind_rows(isolation_rows)

# Flag known outlier tracts
outlier_patterns <- c("MCP", "ML L", "ML R", "CGH L", "CGH R")

isolation_df <- isolation_df %>%
  mutate(is_outlier = tract %in% outlier_patterns)

# Order tracts by median mean_rg across measures (lowest to highest)
tract_order <- isolation_df %>%
  group_by(tract) %>%
  summarise(median_rg = median(mean_rg, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_rg) %>%
  pull(tract)

isolation_df <- isolation_df %>%
  mutate(tract = factor(tract, levels = tract_order))

panel_C <- ggplot(isolation_df,
                  aes(x = tract, y = mean_rg, colour = measure)) +
  # Grey lines connecting dots across measures for each tract
  geom_line(aes(group = tract), colour = "grey80", linewidth = 0.4) +
  # Uniform dot size for all tracts
  geom_point(size = 2.0, shape = 16, alpha = 0.85) +
  # Reference line at overall median
  geom_hline(
    yintercept = median(isolation_df$mean_rg, na.rm = TRUE),
    linetype = "dashed", colour = "grey50", linewidth = 0.4
  ) +
  scale_colour_manual(values = colors_micro, name = "Measure") +
  scale_y_continuous(
    limits = c(NA, 1),
    breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)
  ) +
  theme_fig() +
  theme(
    axis.text.x        = element_text(size = LABEL_SIZE - 1, angle = 45,
                                      hjust = 1),
    panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3),
    legend.position    = "right"
  ) +
  labs(
    title = "C",
    x     = NULL,
    y     = "Mean rg to all other tracts"
  )

# ============================================================
# Save individual panels
# ============================================================
out_dir <- file.path(proj_dir, "Outputs/ExtDataFig3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "ExtDataFig3_A.pdf"),
       plot = panel_A, width = 6, height = 12, device = pdf)
ggsave(file.path(out_dir, "ExtDataFig3_A.png"),
       plot = panel_A, width = 6, height = 12, dpi = 300)

ggsave(file.path(out_dir, "ExtDataFig3_B.pdf"),
       plot = panel_B, width = 8, height = 8, device = pdf)
ggsave(file.path(out_dir, "ExtDataFig3_B.png"),
       plot = panel_B, width = 8, height = 8, dpi = 300)

ggsave(file.path(out_dir, "ExtDataFig3_C.pdf"),
       plot = panel_C, width = 12, height = 5, device = pdf)
ggsave(file.path(out_dir, "ExtDataFig3_C.png"),
       plot = panel_C, width = 12, height = 5, dpi = 300)

# ============================================================
# Combined figure - landscape layout
# Top row:    A (3x4 densities) | B (square heatmap)
# Bottom row: C (27 tracts horizontal, full width)
# ============================================================
top_row <- plot_grid(
  panel_A, panel_B,
  nrow = 1,
  rel_widths = c(1.1, 1),
  align = "h", axis = "tb"
)

combined <- plot_grid(
  top_row, panel_C,
  ncol = 1,
  rel_heights = c(1.2, 0.6)
)

ggsave(file.path(out_dir, "ExtDataFig3_combined.pdf"),
       plot = combined, width = 12, height = 8, device = pdf)
ggsave(file.path(out_dir, "ExtDataFig3_combined.png"),
       plot = combined, width = 12, height = 8, dpi = 300)

message("All panels saved to: ", out_dir)