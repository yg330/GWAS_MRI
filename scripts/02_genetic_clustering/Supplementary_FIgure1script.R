# ===========================================================================
# Supplementary_FIgure1script.R
# Manuscript item: Supplementary Figure 2 (clustering after pruning correlated IDPs)
#
# What this does:
#   Checks the clustering is not driven by near-duplicate IDPs. Greedily
#   prunes IDPs with |rg| above 0.9 and 0.7, re-clusters, and reports the
#   adjusted Rand index against the full Figure 2 Leiden solution.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Supplementary Figure 1: Clustering robustness after removing
# highly genetically correlated IDPs
#
# Demonstrates that the hierarchical clustering structure shown
# in Figure 2A is not driven by redundant, highly correlated
# IDPs. Pruning is performed separately at two thresholds:
#   Panel A: IDPs with pairwise |rg| > 0.9 pruned
#   Panel B: IDPs with pairwise |rg| > 0.7 pruned
#
# Pruning method: greedy redundancy removal. At each iteration,
# among all IDP pairs whose |rg| still exceeds the threshold,
# the IDP involved in the most exceedances is dropped; this
# repeats until no pair exceeds the threshold. This is the same
# logic as caret::findCorrelation() and keeps a maximal
# near-independent subset rather than removing one arbitrary
# member of every over-threshold pair.
#
# Agreement between each pruned solution's Leiden clusters and
# the original Figure 2 Leiden clusters (subset to the retained
# IDPs) is summarised with the Adjusted Rand Index (ARI) and
# printed in each panel subtitle.
#
# Input:  rg_matrix.rds (same file used in Figure 2)
# Output: FigureS1_clustering_robustness.pdf/.png

source("config.R")

library(tidyverse)
library(igraph)
library(ggdendro)
library(cowplot)
library(mclust)   # adjustedRandIndex

# ---------------------
# Global typography (matches Figure 2)
# ---------------------
BASE_SIZE   <- 10
TITLE_SIZE  <- 10
LABEL_SIZE  <- 9
BASE_FAMILY <- ""

theme_figS1 <- function() {
  theme_classic(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    theme(
      plot.title    = element_text(size = TITLE_SIZE, face = "bold",
                                   hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = LABEL_SIZE, hjust = 0,
                                   margin = margin(b = 4)),
      axis.text     = element_text(size = LABEL_SIZE),
      axis.title    = element_text(size = LABEL_SIZE)
    )
}

# ---------------------
# Load and clean data (identical to Figure 2)
# ---------------------
rg_matrix_raw <- readRDS(file.path(proj_dir, "Data", "rg_matrix.rds"))

clean_names <- function(x) gsub("_as_pheno.*", "", gsub("GWAS_", "", x))
rownames(rg_matrix_raw) <- clean_names(rownames(rg_matrix_raw))
colnames(rg_matrix_raw) <- clean_names(colnames(rg_matrix_raw))

exclude_patterns <- c(
  "SubCorVol_choroid.plexus", "SubCorVol_CSF", "SubCorVol_Optic.Chiasm",
  "SubCorVol_SubCortGrayVol", "SubCorVol_vessel", "SubCorVol_WM.hypointensities"
)
keep_idx  <- !rownames(rg_matrix_raw) %in% exclude_patterns
rg_matrix <- rg_matrix_raw[keep_idx, keep_idx]

# ---------------------
# Leiden clustering helper (identical to Figure 2)
# ---------------------
run_leiden <- function(mat, resolution = 0.56, seed = 1, min_size = 10) {
  diag(mat) <- 0
  adj <- (mat + 1) / 2
  adj[adj < 0.3] <- 0
  g <- graph_from_adjacency_matrix(adj, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  set.seed(seed)
  res          <- cluster_leiden(g, weights = E(g)$weight,
                                 resolution_parameter = resolution)
  clusters_raw <- membership(res)
  sizes        <- table(clusters_raw)
  clusters     <- as.character(clusters_raw)
  names(clusters) <- names(clusters_raw)
  for (cid in names(sizes))
    if (sizes[cid] < min_size)
      clusters[clusters_raw == as.numeric(cid)] <- "N/A"
  clusters[is.na(clusters)] <- "N/A"
  non_na <- sort(as.numeric(unique(clusters[clusters != "N/A"])))
  remap  <- setNames(as.character(seq_along(non_na)), as.character(non_na))
  ifelse(clusters == "N/A", "N/A", remap[clusters])
}

# ---------------------
# Original (unpruned) clustering - reference solution
# ---------------------
leiden_original <- run_leiden(rg_matrix)

cluster_palette <- c(
  "1" = "#0072B2", "2" = "#E69F00", "3" = "#009E73",
  "4" = "#CC79A7", "5" = "#56B4E9", "6" = "#D55E00",
  "7" = "#4B0082", "8" = "#666666"
)

# ---------------------
# Greedy redundancy pruning at a given |rg| threshold
# ---------------------
prune_correlated <- function(mat, threshold) {
  mat_prune <- mat
  repeat {
    m <- abs(mat_prune)
    diag(m) <- 0
    over <- which(m > threshold, arr.ind = TRUE)
    if (nrow(over) == 0) break
    tbl       <- table(c(over[, 1], over[, 2]))
    drop_idx  <- as.integer(names(tbl)[which.max(tbl)])
    mat_prune <- mat_prune[-drop_idx, -drop_idx]
  }
  mat_prune
}

# ---------------------
# Build one dendrogram panel for a given threshold
# ---------------------
build_panel <- function(threshold, panel_letter) {
  
  mat_p <- prune_correlated(rg_matrix, threshold)
  n_removed <- nrow(rg_matrix) - nrow(mat_p)
  
  leiden_p <- run_leiden(mat_p)
  
  # Agreement with original clustering, restricted to retained IDPs
  common   <- intersect(rownames(mat_p), names(leiden_original))
  ari_val  <- mclust::adjustedRandIndex(
    leiden_original[common], leiden_p[common]
  )
  
  # Dendrogram on pruned matrix
  dist_p  <- as.dist(1 - mat_p)
  tree_p  <- hclust(dist_p, method = "average")
  dendro  <- ggdendro::dendro_data(tree_p)
  dendro$segments$y    <- pmax(dendro$segments$y, 0)
  dendro$segments$yend <- pmax(dendro$segments$yend, 0)
  
  order_labels <- dendro$labels$label
  
  # Colour tips by ORIGINAL Leiden cluster membership, so panels are
  # visually comparable to Figure 2A and to each other
  orig_cluster_for_tip <- leiden_original[order_labels]
  orig_cluster_for_tip[is.na(orig_cluster_for_tip)] <- "N/A"
  bar_colors <- c(cluster_palette, "N/A" = "#BBBBBB")[orig_cluster_for_tip]
  bar_colors[is.na(bar_colors)] <- "#BBBBBB"
  
  bar_df <- data.frame(x = seq_along(bar_colors), y = 1, fill = bar_colors)
  
  p_dend <- ggplot(segment(dendro)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_figS1() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid   = element_blank(),
          plot.margin  = margin(5, 5, 0, 5, "pt")) +
    labs(
      title    = sprintf("%s  IDPs pruned at |rg| > %.1f", panel_letter, threshold),
      subtitle = sprintf("n = %d IDPs retained (%d removed); Adjusted Rand Index vs. Original cluster = %.2f",
                         nrow(mat_p), n_removed, ari_val),
      x = NULL, y = "1 \u2212 rg"
    ) +
    coord_cartesian(ylim = c(0, max(dendro$segments$y) * 1.02))
  
  p_bar <- ggplot(bar_df, aes(x = x, y = y, fill = fill)) +
    geom_tile() +
    scale_fill_identity() +
    theme_void(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    theme(plot.margin = margin(0, 5, 0, 5, "pt")) +
    labs(title = "Original Leiden cluster")
  
  plot_grid(p_dend, p_bar, ncol = 1, rel_heights = c(5, 0.4),
            align = "v", axis = "lr")
}

panel_A <- build_panel(0.9, "A")
panel_B <- build_panel(0.8, "B")
panel_C <- build_panel(0.7, "C")
panel_D <- build_panel(0.6, "C")

final_plot <- plot_grid(panel_A, panel_B, panel_C, panel_D, ncol = 2)

# ---------------------
# Save
# ---------------------
out_dir <- file.path(proj_dir, "Outputs/Figure_S1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "FigureS1_clustering_robustness.pdf"),
       plot = final_plot, width = 12, height = 10, device = pdf)
ggsave(file.path(out_dir, "FigureS1_clustering_robustness.png"),
       plot = final_plot, width = 12, height = 10, dpi = 300)

message("Figure saved to: ", out_dir)