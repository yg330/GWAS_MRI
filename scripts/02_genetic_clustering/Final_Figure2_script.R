# ===========================================================================
# Final_Figure2_script.R
# Manuscript item: Figure 2
#
# What this does:
#   Clusters the IDPs on genetic correlation distance (1 - rg): average-
#   linkage dendrogram, UMAP embeddings coloured by region / measure /
#   Leiden community, and the mean-rg heatmap between phenotype groups
#   (panel E).
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Figure 2: Genetic clustering of 2,332 structural MRI IDPs
#
# Generates a multi-panel figure comprising:
#   A. Average-linkage hierarchical clustering dendrogram with
#      annotation bars for region type, phenotype, and Leiden cluster
#   B-D. UMAP embeddings coloured by region type, phenotype measure,
#      and Leiden community membership (6 clusters, res=0.55)
#   E. Mean pairwise genetic correlation heatmap between phenotype
#      groups, with subcortical structures stratified into ventricles,
#      corpus callosum, striatal/limbic, and hindbrain/diencephalic
#
# Input:  rg_matrix.rds - full pairwise LDSC genetic correlation matrix
# Output: Figure2_with_heatmap.pdf/.png

source("config.R")

library(tidyverse)
library(RColorBrewer)
library(igraph)
library(umap)
library(ggplot2)
library(grid)
library(cowplot)
library(ggdendro)
library(gridExtra)
library(scales)
library(reshape2)

# ---------------------
# Global typography
# All text in the figure uses these values for consistency
# ---------------------
BASE_SIZE   <- 10     # base font size (pts) for all panels
TITLE_SIZE  <- 10     # panel title size
LABEL_SIZE  <- 9      # axis tick labels, legend text
TILE_SIZE   <- 2.8    # heatmap cell label size (ggplot mm units)
BASE_FAMILY <- ""     # set to e.g. "Helvetica" if needed

theme_fig2 <- function() {
  theme_classic(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    theme(
      plot.title       = element_text(size = TITLE_SIZE, face = "bold",
                                      hjust = 0, margin = margin(b = 4)),
      axis.text        = element_text(size = LABEL_SIZE),
      axis.title       = element_text(size = LABEL_SIZE),
      legend.text      = element_text(size = LABEL_SIZE),
      legend.title     = element_text(size = LABEL_SIZE),
      legend.key.size  = unit(0.35, "cm"),
      strip.text       = element_text(size = LABEL_SIZE)
    )
}

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

# ---------------------
# Trait classification
# ---------------------
get_main_category <- function(name) {
  if (grepl("^SubCorVol_CC_", name))                 return("CorpusCallosum")
  if (grepl("^SubCorVol_.*(Ventricle|Vent)$", name)) return("Ventricle")
  strsplit(name, "_")[[1]][1]
}

get_measure_type <- function(name) {
  if (grepl("^SubCorVol_CC_", name))                 return("CorpusCallosum")
  if (grepl("^SubCorVol_.*(Ventricle|Vent)$", name)) return("Ventricle")
  parts <- strsplit(name, "_")[[1]]
  first <- parts[1]
  if (first %in% c("global", "regional", "tract") && length(parts) >= 2) {
    if (parts[2] == "GMV" || first == "SubCorVol") return("GMV_SubCorVol") else return(parts[2])
  } else if (first == "SubCorVol") return("SubCorVol")
  return(first)
}

# ---------------------
# Subcortical group assignment
# Renamed per paragraph:
#   Striatal/Limbic: caudate, accumbens, putamen, amygdala, hippocampus
#   Hindbrain/Diencephalic: thalamus, pallidum, ventral diencephalon, brainstem, cerebellum
# ---------------------
get_subcor_group <- function(name) {
  if (!grepl("^SubCorVol_", name)) return(NA)
  if (grepl("^SubCorVol_CC_", name))                 return("Corpus callosum")
  if (grepl("^SubCorVol_.*(Ventricle|Vent)$", name)) return("Ventricles")
  # Striatal/limbic
  if (grepl("Caudate|Accumbens|Putamen|Amygdala|Hippocampus",
            name, ignore.case = TRUE))                return("Striatal/limbic")
  # Hindbrain/diencephalic
  if (grepl("Thalamus|Pallidum|VentralDC|Brain.Stem|Cerebellum",
            name, ignore.case = TRUE))                return("Hindbrain/\ndiencephalic")
  return("Other subcortical")
}

trait_names     <- colnames(rg_matrix)
main_categories <- sapply(trait_names, get_main_category)
measure_types   <- sapply(trait_names, get_measure_type)

# ---------------------
# Between-groups assignment for heatmap
# ---------------------
get_phenotype_group <- function(name, measure) {
  # Subcortical subgroups
  if (grepl("^SubCorVol_CC_", name))                 return("Corpus\ncallosum")
  if (grepl("^SubCorVol_.*(Ventricle|Vent)$", name)) return("Ventricles")
  if (grepl("^SubCorVol_", name)) {
    if (grepl("Caudate|Accumbens|Putamen|Amygdala|Hippocampus",
              name, ignore.case = TRUE)) return("Striatal/\nlimbic")
    if (grepl("Thalamus|Pallidum|VentralDC|Brain.Stem|Cerebellum",
              name, ignore.case = TRUE)) return("Hindbrain/\ndiencephalic")
    return("Other\nsubcortical")
  }
  # Cortical/tract measures - use measure type
  measure
}

phenotype_groups <- mapply(get_phenotype_group, trait_names, measure_types)

# Desired order of groups in heatmap
group_order <- c(
  "CT", "FA", "ICVF", "Ventricles", "ISOVF", "MD",
  "Corpus\ncallosum", "IC", "FI", "GMV_SubCorVol", "SA",
  "Hindbrain/\ndiencephalic", "Striatal/\nlimbic",
  "OD", "GC", "MC"
)

# Display labels (cleaner names for axis)
group_labels <- c(
  "CT"                    = "CT",
  "FA"                    = "FA",
  "ICVF"                  = "ICVF",
  "Ventricles"            = "Ventricles",
  "ISOVF"                 = "ISOVF",
  "MD"                    = "MD",
  "Corpus\ncallosum"      = "Corpus\ncallosum",
  "IC"                    = "IC",
  "FI"                    = "FI",
  "GMV_SubCorVol"         = "GMV",
  "SA"                    = "SA",
  "Hindbrain/\ndiencephalic" = "Hindbrain/\ndiencephalic",
  "Striatal/\nlimbic"     = "Striatal/\nlimbic",
  "OD"                    = "OD",
  "GC"                    = "GC",
  "MC"                    = "MC"
)

# Compute mean rg between groups
groups_present <- unique(phenotype_groups)
groups_present <- groups_present[groups_present %in% group_order]
groups_ordered <- group_order[group_order %in% groups_present]

n_groups <- length(groups_ordered)
mean_rg_mat <- matrix(NA, nrow = n_groups, ncol = n_groups,
                      dimnames = list(groups_ordered, groups_ordered))

for (g1 in groups_ordered) {
  for (g2 in groups_ordered) {
    idx1 <- which(phenotype_groups == g1)
    idx2 <- which(phenotype_groups == g2)
    if (length(idx1) == 0 || length(idx2) == 0) next
    if (g1 == g2) {
      # Within-group: upper triangle only
      sub  <- rg_matrix[idx1, idx2]
      vals <- sub[upper.tri(sub)]
    } else {
      vals <- as.vector(rg_matrix[idx1, idx2])
    }
    mean_rg_mat[g1, g2] <- mean(vals, na.rm = TRUE)
  }
}

# Melt for ggplot
heatmap_df <- melt(mean_rg_mat, varnames = c("Group1", "Group2"),
                   value.name = "mean_rg") %>%
  mutate(
    Group1 = factor(Group1, levels = rev(groups_ordered)),
    Group2 = factor(Group2, levels = groups_ordered),
    label  = ifelse(!is.na(mean_rg), sprintf("%.2f", mean_rg), "")
  )

# Colour blocks for subcortical annotation
subcor_fill <- c(
  "CT" = NA, "FA" = NA, "ICVF" = NA, "ISOVF" = NA, "MD" = NA,
  "IC" = NA, "FI" = NA, "GMV_SubCorVol" = NA, "SA" = NA,
  "OD" = NA, "GC" = NA, "MC" = NA,
  "Ventricles"               = "#4169E1",
  "Corpus\ncallosum"         = "#8B0000",
  "Striatal/\nlimbic"        = "#6A3D9A",
  "Hindbrain/\ndiencephalic" = "#9370DB"
)

# Heatmap panel
p_heatmap <- ggplot(heatmap_df, aes(x = Group2, y = Group1, fill = mean_rg)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = TILE_SIZE, colour = "black",
            family = BASE_FAMILY) +
  scale_fill_gradient2(
    low      = "#4575B4",
    mid      = "white",
    high     = "#D73027",
    midpoint = 0,
    limits   = c(-0.5, 1),
    na.value = "grey90",
    name     = "Mean rg"
  ) +
  scale_x_discrete(labels = group_labels) +
  scale_y_discrete(labels = group_labels) +
  theme_fig2() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = LABEL_SIZE),
    axis.text.y      = element_text(size = LABEL_SIZE),
    axis.title       = element_blank(),
    panel.border     = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "right"
  ) +
  labs(title = "E  Mean pairwise genetic correlations between phenotype groups")

# ---------------------
# Leiden clustering
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

leiden_clusters <- run_leiden(rg_matrix)

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

all_main_cats <- unique(main_categories)
base_palette  <- c(
  "#E31A1C","#FDBF6F","#CAB2D6","#33A02C","#FF7F00","#6A3D9A","#B15928",
  "#A6CEE3","#FFFF99","#B2DF8A","#FB9A99","#8B0000","#4169E1"
)
colors_main <- setNames(base_palette[seq_along(all_main_cats)], all_main_cats)
colors_main["CorpusCallosum"] <- "#8B0000"
colors_main["Ventricle"]      <- "#4169E1"

cluster_palette <- c(
  "1" = "#0072B2", "2" = "#E69F00", "3" = "#009E73",
  "4" = "#CC79A7", "5" = "#56B4E9", "6" = "#D55E00",
  "7" = "#4B0082", "8" = "#666666"
)
unique_clusters <- c(as.character(sort(as.numeric(
  unique(leiden_clusters[leiden_clusters != "N/A"])))), "N/A")
clusters_no_na  <- unique_clusters[unique_clusters != "N/A"]
colors_cluster  <- c(cluster_palette[clusters_no_na], "N/A" = "#BBBBBB")

# ---------------------
# Dendrogram (Panel A)
# ---------------------
rg.distance <- as.dist(1 - rg_matrix)
rg.tree     <- hclust(rg.distance, method = "average")
dendro_data <- ggdendro::dendro_data(rg.tree)
order_labels <- dendro_data$labels$label

dendro_data$segments$y    <- pmax(dendro_data$segments$y, 0)
dendro_data$segments$yend <- pmax(dendro_data$segments$yend, 0)

ordered_main_cat     <- main_categories[order_labels]
ordered_measure_type <- measure_types[order_labels]
ordered_leiden       <- leiden_clusters[order_labels]
ordered_leiden[is.na(ordered_leiden)] <- "N/A"

bar_colors_main    <- colors_main[ordered_main_cat]
bar_colors_measure <- colors_measure[ordered_measure_type]
bar_colors_cluster <- colors_cluster[ordered_leiden]

bar_colors_main[is.na(bar_colors_main)]       <- "#CCCCCC"
bar_colors_measure[is.na(bar_colors_measure)] <- "#CCCCCC"
bar_colors_cluster[is.na(bar_colors_cluster)] <- "#333333"

p_dend <- ggplot(segment(dendro_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_fig2() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = LABEL_SIZE),
        panel.grid   = element_blank(),
        plot.margin  = margin(5, 5, 0, 5, "pt")) +
  labs(title = "A  Hierarchical clustering dendrogram", x = NULL, y = "1 \u2212 rg") +
  coord_cartesian(ylim = c(0, max(dendro_data$segments$y) * 1.02))

make_color_bar <- function(colors, label) {
  df <- data.frame(x = seq_along(colors), y = 1, fill = colors)
  ggplot(df, aes(x = x, y = y, fill = fill)) +
    geom_tile() +
    scale_fill_identity() +
    theme_void(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    theme(legend.position = "none",
          plot.title  = element_text(hjust = 0, size = LABEL_SIZE,
                                     family = BASE_FAMILY),
          plot.margin = margin(0, 5, 0, 5, "pt")) +
    labs(title = label)
}

p_bar_main    <- make_color_bar(bar_colors_main,    "Region type")
p_bar_measure <- make_color_bar(bar_colors_measure, "Phenotype")
p_bar_cluster <- make_color_bar(bar_colors_cluster, "Leiden")

panel_A <- plot_grid(
  p_dend, p_bar_main, p_bar_measure, p_bar_cluster,
  ncol = 1, rel_heights = c(5, 0.4, 0.4, 0.4),
  align = "v", axis = "lr")

# ---------------------
# UMAP - increased transparency (alpha 0.55)
# ---------------------
cfg              <- umap.defaults
cfg$n_neighbors  <- 10
cfg$min_dist     <- 0.05
cfg$random_state <- 42
set.seed(42)
umap_result <- umap(1 - rg_matrix, config = cfg)

umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df <- cbind(
  umap_df,
  Trait        = trait_names,
  Cluster      = factor(leiden_clusters),
  MainCategory = factor(main_categories),
  MeasureType  = factor(measure_types)
)

plot_umap <- function(data, color_by, palette, label) {
  ggplot(data, aes(x = UMAP1, y = UMAP2, color = .data[[color_by]])) +
    geom_point(size = 1.2, alpha = 0.75) +
    scale_color_manual(values = palette, name = NULL) +
    theme_fig2() +
    theme(legend.position = "right") +
    labs(title = label, x = "UMAP 1", y = "UMAP 2")
}

panel_B <- plot_umap(umap_df, "MainCategory", colors_main,
                     "B  UMAP by region type")
panel_C <- plot_umap(umap_df, "MeasureType",  colors_measure,
                     "C  UMAP by phenotype")
panel_D <- plot_umap(umap_df, "Cluster",      colors_cluster,
                     "D  UMAP by Leiden cluster")

# ---------------------
# Assemble figure - two-column layout
#
# Left column:  A (dendrogram) stacked above E (heatmap)
# Right column: B, C, D (UMAPs) stacked vertically
# ---------------------

# Left column: dendrogram on top, heatmap below
left_col <- plot_grid(
  panel_A, p_heatmap,
  ncol = 1,
  rel_heights = c(1, 1.4),
  align = "v", axis = "lr"
)

# Right column: three UMAPs stacked
right_col <- plot_grid(
  panel_B, panel_C, panel_D,
  ncol = 1,
  rel_heights = c(1, 1, 1),
  align = "v", axis = "lr"
)

# Combine columns - left slightly wider to accommodate heatmap labels
final_plot <- plot_grid(
  left_col, right_col,
  nrow = 1,
  rel_widths = c(1.3, 1)
)

# ---------------------
# Save
# ---------------------
out_dir <- file.path(proj_dir, "Outputs/Figure_2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "Figure2_with_heatmap.pdf"),
       plot = final_plot, width = 14, height = 10, device = pdf)
ggsave(file.path(out_dir, "Figure2_with_heatmap.png"),
       plot = final_plot, width = 14, height = 10, dpi = 300)

message("Figure saved to: ", out_dir)