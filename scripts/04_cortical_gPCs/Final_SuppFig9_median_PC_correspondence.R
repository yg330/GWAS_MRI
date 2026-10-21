# ===========================================================================
# Final_SuppFig9_median_PC_correspondence.R
# Manuscript item: Supplementary Figure 9
#
# What this does:
#   Clustered similarity of the gPC1-3 loadings across measures.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Supplementary Figure 9 -- cross-measure correspondence within
# each of gPC1, gPC2, gPC3, clustered
#
# Faceted by component, not by measure: THREE panels (gPC1, gPC2,
# gPC3), each a 12x12 measure x measure heatmap of |Spearman rho|
# between that SAME component's parcel-wise loadings across every
# pair of the 12 measures (sign is arbitrary per measure/component --
# see Final_ExtDataFig4_PCA_master.R STEP 3 -- so |rho| is what's
# interpretable). Within each panel, the 12 measures are
# independently reordered via hierarchical clustering (hclust,
# average linkage, distance = 1 - |rho|) on THAT component's own
# similarity matrix -- so gPC1's measure ordering can legitimately
# differ from gPC2's or gPC3's, rather than forcing one fixed
# measure order across all three panels.
#
# Inputs:
#   ~/Downloads/rstudio-export (3)/rg_matrix.rds   (PC1-3 recomputed
#   fresh here, identical construction to Final_ExtDataFig4_PCA_
#   master.R STEP 1, so this script is fully standalone)
#
# Outputs (under out_dir):
#   PC1_PC2_PC3_measure_similarity_long.csv  -- long-format table,
#     one row per gPC x measure_1 x measure_2, with |rho|
#   SuppFig9_gPC_clustered_similarity.pdf / .png

source("config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
})

ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# ---------------- USER INPUTS ----------------
rg_path    <- file.path(proj_dir, "Data", "rg_matrix.rds")
base_dir   <- proj_dir
out_dir    <- file.path(base_dir, "Outputs/SuppFig9_median_PC_correspondence")
all_phenos <- c("CT", "FA", "FI", "GC", "GMV", "IC", "ICVF", "ISOVF", "MC", "MD", "OD", "SA")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ================================================================
# STEP 1 -- PC1-3 loadings (identical construction to
# Final_ExtDataFig4_PCA_master.R STEP 1 / Final_Figure3_script.R)
# ================================================================
ts("==== STEP 1: PC1-3 loadings ====")
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
  scores <- pca_res$x[, 1:n_pc, drop = FALSE]
  colnames(scores) <- paste0("PC", 1:n_pc)
  for (j in 1:n_pc) if (mean(scores[, j]) < 0) scores[, j] <- -scores[, j]
  rownames(scores) <- sapply(traits, extract_parcel)
  scores
}

pc123_res <- Filter(Negate(is.null), setNames(lapply(all_phenos, pheno_pca_pc123), all_phenos))

pc123_long <- map_dfr(names(pc123_res), function(ph) {
  sc <- pc123_res[[ph]]
  tibble(parcel = rownames(sc), measure = ph, PC1 = sc[, "PC1"], PC2 = sc[, "PC2"], PC3 = sc[, "PC3"])
})

measures_present <- sort(unique(pc123_long$measure))
cat(sprintf("  %d measures with PC1-3: %s\n", length(measures_present), paste(measures_present, collapse = ", ")))

# ================================================================
# STEP 2 -- for each gPC, a 12x12 measure x measure |rho| matrix
# ================================================================
ts("==== STEP 2: per-gPC cross-measure |rho| ====")

pc_names <- c("PC1", "PC2", "PC3")

cross_measure_similarity <- function(pc_name) {
  wide <- pc123_long %>%
    dplyr::select(parcel, measure, value = all_of(pc_name)) %>%
    pivot_wider(names_from = measure, values_from = value)
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$parcel
  abs(cor(mat, method = "spearman", use = "pairwise.complete.obs"))
}

sim_list <- setNames(lapply(pc_names, cross_measure_similarity), pc_names)

sim_long <- imap_dfr(sim_list, function(mat, pc) {
  as_tibble(mat, rownames = "measure_1") %>%
    pivot_longer(-measure_1, names_to = "measure_2", values_to = "abs_rho") %>%
    mutate(gPC = pc)
})
write.csv(sim_long, file.path(out_dir, "PC1_PC2_PC3_measure_similarity_long.csv"), row.names = FALSE)
cat("  Wrote PC1_PC2_PC3_measure_similarity_long.csv\n")

# ================================================================
# STEP 3 -- within EACH gPC, cluster the 12 measures independently
# (hclust, average linkage, distance = 1 - |rho|) -- gPC1's order
# need not match gPC2's or gPC3's
# ================================================================
ts("==== STEP 3: clustering measures within each gPC ====")

measure_order_by_pc <- setNames(lapply(pc_names, function(pc) {
  mat <- sim_list[[pc]]
  hc <- hclust(as.dist(1 - mat), method = "average")
  ord <- hc$labels[hc$order]
  cat(sprintf("  %s clustered order: %s\n", pc, paste(ord, collapse = " -> ")))
  ord
}), pc_names)

# ================================================================
# STEP 4 -- one 12x12 heatmap per gPC, rows/cols in that gPC's own
# clustered order
# ================================================================
ts("==== STEP 4: building panels ====")

plot_one_pc <- function(pc) {
  ord <- measure_order_by_pc[[pc]]
  df <- sim_long %>%
    dplyr::filter(gPC == pc) %>%
    mutate(measure_1 = factor(measure_1, levels = ord),
           measure_2 = factor(measure_2, levels = rev(ord)))
  ggplot(df, aes(x = measure_1, y = measure_2, fill = abs_rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", abs_rho)), size = 2.2) +
    scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", limits = c(0, 1), name = expression("|" * rho * "|")) +
    labs(x = NULL, y = NULL, title = pc) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
          axis.text.y = element_text(size = 7.5))
}

panels <- map(pc_names, plot_one_pc)

# shared legend, extracted from one panel rebuilt with a legend
legend_src <- plot_one_pc(pc_names[1]) +
  theme(legend.position = "right", legend.title = element_text(size = 9), legend.text = element_text(size = 8))
shared_legend <- get_legend(legend_src)

grid_plot <- plot_grid(plotlist = panels, nrow = 1, ncol = 3, align = "hv")
fig9 <- plot_grid(grid_plot, shared_legend, nrow = 1, rel_widths = c(1, 0.06))

fig9 <- plot_grid(
  ggdraw() + draw_label(
    "Cross-measure correspondence within gPC1, gPC2, gPC3\n(12x12 |rho| between measures' loadings on that component; each panel independently clustered via hclust)",
    fontface = "bold", size = 11, x = 0.5, hjust = 0.5),
  fig9, ncol = 1, rel_heights = c(0.08, 1)
)

ggsave(file.path(out_dir, "SuppFig9_gPC_clustered_similarity.pdf"), fig9, width = 16, height = 6.5)
ggsave(file.path(out_dir, "SuppFig9_gPC_clustered_similarity.png"), fig9, width = 16, height = 6.5, dpi = 300)

ts("Saved:", file.path(out_dir, "SuppFig9_gPC_clustered_similarity.pdf"))
cat("\n=== ALL DONE ===\n")
