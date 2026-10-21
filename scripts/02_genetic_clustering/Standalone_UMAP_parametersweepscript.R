# ===========================================================================
# Standalone_UMAP_parametersweepscript.R
# Manuscript item: Figure 2 (UMAP settings)
#
# What this does:
#   Sweeps UMAP n_neighbors and min_dist to see how sensitive the embedding
#   is; used to choose the settings in Figure 2.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# UMAP exploration script
# Standalone - loads and classifies data from scratch
# Sweeps n_neighbors and min_dist to find best separation

source("config.R")

library(tidyverse)
library(umap)
library(ggplot2)
library(cowplot)
library(igraph)

# ---------------------
# Load and clean data
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
cat(sprintf("Matrix: %d x %d\n", nrow(rg_matrix), ncol(rg_matrix)))

# ---------------------
# Classify traits
# ---------------------
get_measure_type <- function(name) {
  if (grepl("^SubCorVol_CC_", name))                 return("CorpusCallosum")
  if (grepl("^SubCorVol_.*(Ventricle|Vent)$", name)) return("Ventricle")
  parts <- strsplit(name, "_")[[1]]
  first <- parts[1]
  if (first %in% c("global", "regional", "tract") && length(parts) >= 2) {
    if (parts[2] == "GMV") return("GMV_SubCorVol") else return(parts[2])
  } else if (first == "SubCorVol") return("SubCorVol")
  return(first)
}

trait_names   <- colnames(rg_matrix)
measure_types <- sapply(trait_names, get_measure_type)

# ---------------------
# Leiden clustering (res=0.55, min_size=10)
# ---------------------
run_leiden <- function(mat, resolution = 0.55, seed = 1, min_size = 10) {
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

leiden_clusters <- run_leiden(rg_matrix)
cat("Leiden cluster sizes:\n")
print(sort(table(leiden_clusters), decreasing = TRUE))

# ---------------------
# Colour palettes
# ---------------------
colors_measure <- c(
  CT             = "#0072B2",
  FA             = "#C0006E",
  FI             = "#7DDCB8",
  GC             = "#6B4800",
  GMV_SubCorVol  = "#003D28",
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

cluster_palette <- c(
  "1" = "#0072B2", "2" = "#E69F00", "3" = "#009E73",
  "4" = "#CC79A7", "5" = "#56B4E9", "6" = "#D55E00",
  "7" = "#4B0082", "8" = "#666666", "N/A" = "#BBBBBB"
)

# ---------------------
# UMAP parameter sweep
# ---------------------
dist_mat <- 1 - rg_matrix

neighbors_vals <- c(5, 10, 15, 30)
min_dist_vals  <- c(0.01, 0.05, 0.1, 0.3)

plots_measure <- list()
plots_cluster <- list()

for (nn in neighbors_vals) {
  for (md in min_dist_vals) {
    
    cfg            <- umap.defaults
    cfg$n_neighbors <- nn
    cfg$min_dist    <- md
    cfg$random_state <- 42
    
    set.seed(42)
    layout <- umap(dist_mat, config = cfg)$layout
    df <- data.frame(
      UMAP1   = layout[, 1],
      UMAP2   = layout[, 2],
      Measure = measure_types[trait_names],
      Cluster = leiden_clusters[trait_names]
    )
    
    label <- sprintf("nn=%d  md=%.2f", nn, md)
    
    p_m <- ggplot(df, aes(x = UMAP1, y = UMAP2, colour = Measure)) +
      geom_point(size = 0.8, alpha = 0.8) +
      scale_colour_manual(values = colors_measure) +
      theme_classic(base_size = 9) +
      theme(legend.position = "none",
            plot.title = element_text(size = 8)) +
      ggtitle(label)
    
    p_c <- ggplot(df, aes(x = UMAP1, y = UMAP2, colour = Cluster)) +
      geom_point(size = 0.8, alpha = 0.8) +
      scale_colour_manual(values = cluster_palette) +
      theme_classic(base_size = 9) +
      theme(legend.position = "none",
            plot.title = element_text(size = 8)) +
      ggtitle(label)
    
    key <- paste0("nn", nn, "_md", md)
    plots_measure[[key]] <- p_m
    plots_cluster[[key]] <- p_c
    
    cat(sprintf("Done: %s\n", label))
  }
}

# ---------------------
# Assemble sweep grids
# ---------------------
grid_measure <- plot_grid(plotlist = plots_measure,
                          ncol = length(min_dist_vals),
                          nrow = length(neighbors_vals))

grid_cluster <- plot_grid(plotlist = plots_cluster,
                          ncol = length(min_dist_vals),
                          nrow = length(neighbors_vals))

# ---------------------
# Save
# ---------------------
out_dir <- file.path(proj_dir, "Outputs/Figure_2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "umap_sweep_measure.pdf"),
       plot = grid_measure, width = 16, height = 14, device = pdf)
ggsave(file.path(out_dir, "umap_sweep_measure.png"),
       plot = grid_measure, width = 16, height = 14, dpi = 200)

ggsave(file.path(out_dir, "umap_sweep_cluster.pdf"),
       plot = grid_cluster, width = 16, height = 14, device = pdf)
ggsave(file.path(out_dir, "umap_sweep_cluster.png"),
       plot = grid_cluster, width = 16, height = 14, dpi = 200)

cat("\nSweep grids saved — check umap_sweep_measure.png and umap_sweep_cluster.png\n")

# ---------------------
# Single best-guess plot with legend
# (nn=10, md=0.05 - good starting point)
# ---------------------
cfg_best             <- umap.defaults
cfg_best$n_neighbors <- 10
cfg_best$min_dist    <- 0.05
cfg_best$random_state <- 42

set.seed(42)
layout_best <- umap(dist_mat, config = cfg_best)$layout
df_best <- data.frame(
  UMAP1   = layout_best[, 1],
  UMAP2   = layout_best[, 2],
  Measure = measure_types[trait_names],
  Cluster = leiden_clusters[trait_names]
)

p_best_measure <- ggplot(df_best, aes(x = UMAP1, y = UMAP2, colour = Measure)) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_colour_manual(values = colors_measure) +
  theme_classic(base_size = 12) +
  labs(title = "UMAP by phenotype (nn=10, min_dist=0.05)")

p_best_cluster <- ggplot(df_best, aes(x = UMAP1, y = UMAP2, colour = Cluster)) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_colour_manual(values = cluster_palette) +
  theme_classic(base_size = 12) +
  labs(title = "UMAP by Leiden cluster (nn=10, min_dist=0.05)")

p_best <- plot_grid(p_best_measure, p_best_cluster, nrow = 1)

ggsave(file.path(out_dir, "umap_best.pdf"),
       plot = p_best, width = 16, height = 7, device = pdf)
ggsave(file.path(out_dir, "umap_best.png"),
       plot = p_best, width = 16, height = 7, dpi = 300)

cat("Best single plot saved — umap_best.png\n")
cat("All outputs in: ", out_dir, "\n")