# ===========================================================================
# Final_leiden_searchscript.R
# Manuscript item: Figure 2 (Leiden seed search)
#
# What this does:
#   Runs Leiden clustering over a range of random seeds and writes the
#   cluster assignments side by side, so the seed used in Figure 2 can be
#   chosen and documented.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# leiden_seed_search.R
# Runs Leiden clustering across a range of seeds and exports a
# CSV table with one column per seed so you can inspect cluster
# compositions and pick the seed that best matches the expected
# biological groupings (e.g. CT isolated in its own cluster).
#
# OUTPUT: leiden_seed_comparison.csv  (same folder as this script,
#         or change `out_dir` below)
#
# Once you identify the right seed, set it in Fig2_script_updated.R
# at the line:   set.seed(<your_seed>)   above cluster_leiden().

source("config.R")

library(tidyverse)
library(igraph)

# ---------------------
# Paths - adjust as needed
# ---------------------
rg_path <- file.path(proj_dir, "Data", "rg_matrix.rds")
out_dir  <- file.path(proj_dir, "Outputs/Figure_2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------
# Load & clean matrix
# ---------------------
rg_matrix <- readRDS(rg_path)
clean_names <- function(x) gsub("_as_pheno.*", "", gsub("GWAS_", "", x))
rownames(rg_matrix) <- clean_names(rownames(rg_matrix))
colnames(rg_matrix) <- clean_names(colnames(rg_matrix))

# Measure type extraction (matches main script)
get_measure_type <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  first <- parts[1]
  if (first %in% c("global", "regional", "tract") && length(parts) >= 2) {
    if (parts[2] == "GMV" || first == "SubCorVol") return("GMV_SubCorVol")
    else return(parts[2])
  } else if (first == "SubCorVol") return("GMV_SubCorVol")
  return(first)
}
measure_types <- sapply(colnames(rg_matrix), get_measure_type)

# ---------------------
# Build graph (identical to main script)
# ---------------------
diag(rg_matrix) <- 0
adj_matrix <- (rg_matrix + 1) / 2
adj_matrix[adj_matrix < 0.3] <- 0
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected",
                                 weighted = TRUE, diag = FALSE)

# ---------------------
# Seeds to try
# Extend this vector as needed - add more seeds if none of these
# reproduce the expected clustering.
# ---------------------
seeds_to_try <- c(1, 7, 42, 99, 123, 314, 500, 1234, 2024,
                  2025, 9999, 42042, 11, 22, 33, 55, 77, 88,
                  100, 200, 300, 400, 999, 1111, 5555)

# ---------------------
# Run Leiden for each seed
# ---------------------
run_leiden <- function(seed) {
  set.seed(seed)
  lr  <- cluster_leiden(g, weights = E(g)$weight,
                        resolution_parameter = 0.56)
  raw <- as.character(membership(lr))
  names(raw) <- names(membership(lr))
  
  # Collapse small clusters (<10 members) to N/A
  sz  <- table(raw)
  raw[raw %in% names(sz[sz < 10])] <- "N/A"
  
  # Renumber non-N/A clusters 1, 2, 3, ... in ascending order
  non_na <- sort(as.numeric(unique(raw[raw != "N/A"])))
  rmap   <- setNames(as.character(seq_along(non_na)), as.character(non_na))
  raw    <- ifelse(raw == "N/A", "N/A", rmap[raw])
  raw
}

cat("Running Leiden for", length(seeds_to_try), "seeds...\n")
seed_results <- setNames(
  lapply(seeds_to_try, run_leiden),
  paste0("seed_", seeds_to_try)
)

# ---------------------
# Build comparison table
# ---------------------
trait_order <- names(seed_results[[1]])

seed_df <- as.data.frame(seed_results, row.names = trait_order)

# Prepend phenotype column
seed_df <- cbind(
  phenotype = measure_types[match(trait_order, names(measure_types))],
  seed_df
)

# Prepend summary rows
n_clusters_row  <- c("n_clusters",
                     sapply(seed_results, function(x) length(unique(x[x != "N/A"]))))
traits_in_na_row <- c("n_traits_NA",
                      sapply(seed_results, function(x) sum(x == "N/A")))

seed_df <- rbind(
  setNames(n_clusters_row,   names(seed_df)),
  setNames(traits_in_na_row, names(seed_df)),
  seed_df
)

# ---------------------
# Save
# ---------------------
out_path <- file.path(out_dir, "leiden_seed_comparison.csv")
write.csv(seed_df, out_path, row.names = TRUE)
cat("Saved:", out_path, "\n")

# ---------------------
# Quick console summary
# ---------------------
cat("\n=== Cluster count per seed ===\n")
summary_tbl <- data.frame(
  seed       = seeds_to_try,
  n_clusters = sapply(seed_results, function(x) length(unique(x[x != "N/A"]))),
  n_NA       = sapply(seed_results, function(x) sum(x == "N/A"))
)
print(summary_tbl)

# Flag seeds where CT traits are all in a single exclusive cluster
cat("\n=== Seeds where CT is in its own isolated cluster ===\n")
for (s in names(seed_results)) {
  cl   <- seed_results[[s]]
  ct_clusters  <- unique(cl[measure_types[names(cl)] == "CT"])
  ct_clusters  <- ct_clusters[ct_clusters != "N/A"]
  if (length(ct_clusters) == 1) {
    cluster_members <- measure_types[names(cl)[cl == ct_clusters]]
    if (all(cluster_members == "CT")) {
      cat(sprintf("  %s : CT is alone in cluster %s (%d traits)\n",
                  s, ct_clusters, sum(cl == ct_clusters, na.rm = TRUE)))
    }
  }
}
cat("(No output above = no seed produced a CT-exclusive cluster;\n",
    " consider widening seeds_to_try or adjusting resolution_parameter)\n")