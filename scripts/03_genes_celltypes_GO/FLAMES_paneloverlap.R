# ===========================================================================
# FLAMES_paneloverlap.R
# Manuscript item: Results: clinical gene panels
#
# What this does:
#   Tests whether FLAMES genes overlap clinical cortical malformation /
#   microcephaly gene panels (ACMG, PanelApp), overall and per family,
#   against a protein-coding background.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# FLAMES x clinical gene panel overlap (ACMG + Genomics England
# PanelApp: cortical malformations + severe microcephaly)
# Combines all three gene panels into one list (they
# are all broadly related cephalic/cortical malformation panels),
# then tests overlap with FLAMES-prioritised genes: once overall,
# and once per phenotype family, using a hypergeometric/Fisher's
# exact test against the MAGMA protein-coding background.

source("config.R")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

# ---- Paths ----
base_dir <- proj_dir

acmg_path        <- file.path(base_dir, "Data", "ACMG_brainmalformations.csv")
cortical_path    <- file.path(base_dir, "Data", "Malformations of cortical development.tsv")
microcephaly_path <- file.path(base_dir, "Data", "Severe microcephaly.tsv")

flames_path   <- file.path(proj_dir, "Data", "FLAMES_gene_prioritisation_results.txt")
gene_loc_path <- file.path(base_dir, "Data", "NCBI37.3.gene.loc")

out_dir <- file.path(base_dir, "Outputs", "gene_panel_overlap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Family map (matches the convention used throughout this project) ----
family_map <- c(
  SA = "Brain size", FI = "Brain size", GMV = "Brain size", IC = "Brain size",
  SubCorVol = "Brain size",
  CT = "Cortical thickness",
  GC = "Curvature", MC = "Curvature",
  FA = "Microstructural coherence", ICVF = "Microstructural coherence",
  MD = "Microstructural diffusivity", ISOVF = "Microstructural diffusivity",
  OD = "Orientation dispersion"
)

# ================================================================
# SECTION 1: Load and combine the three gene panels
# ================================================================
cat("Loading gene panels...\n")

acmg_raw <- read_csv(acmg_path, show_col_types = FALSE)

# Gene column is formatted "SYMBOLHGNC:12345" with no separator --
# strip everything from "HGNC" onward to recover the bare symbol.
# Restricted to Classification == "Definitive" per request.
acmg_genes <- acmg_raw %>%
  dplyr::filter(Classification == "Definitive") %>%
  mutate(symbol = str_remove(Gene, "HGNC.*$") %>% str_trim() %>% toupper()) %>%
  pull(symbol)
cat(sprintf("  ACMG_brainmalformations.csv (Definitive only): %d / %d rows -> %d unique genes\n",
            sum(acmg_raw$Classification == "Definitive"), nrow(acmg_raw), n_distinct(acmg_genes)))

read_panelapp_tsv <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)
  if ("Entity type" %in% colnames(df)) {
    df <- df %>% dplyr::filter(`Entity type` == "gene")
  }
  df %>% pull(`Gene Symbol`) %>% str_trim() %>% toupper()
}

cortical_genes     <- read_panelapp_tsv(cortical_path)
microcephaly_genes <- read_panelapp_tsv(microcephaly_path)
cat(sprintf("  Malformations of cortical development.tsv: %d genes\n", n_distinct(cortical_genes)))
cat(sprintf("  Severe microcephaly.tsv: %d genes\n", n_distinct(microcephaly_genes)))

panel_genes <- unique(c(acmg_genes, cortical_genes, microcephaly_genes))
panel_genes <- panel_genes[!is.na(panel_genes) & panel_genes != ""]
cat(sprintf("  Combined (union, deduplicated): %d unique genes\n", length(panel_genes)))

write.csv(tibble(symbol = sort(panel_genes)), file.path(out_dir, "combined_gene_panel.csv"),
          row.names = FALSE)

# ================================================================
# SECTION 2: Load FLAMES genes and background
# ================================================================
cat("\nLoading FLAMES genes and MAGMA background...\n")

gene_loc <- read.table(gene_loc_path, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("entrez_id", "chr", "start", "end", "strand", "symbol")) %>%
  dplyr::filter(chr %in% c(1:22, "X", "Y")) %>%
  distinct(symbol, .keep_all = TRUE)
background_genes <- toupper(unique(gene_loc$symbol))
cat(sprintf("  MAGMA background: %d genes\n", length(background_genes)))

flames_raw <- read.table(flames_path, header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::filter(grepl("^regional_|^global_|^tract_|^SubCor_Vol", phenotype)) %>%
  mutate(
    phenotype_clean = str_remove(phenotype, "^(regional_|global_|tract_|SubCor_Vol_)"),
    phenotype_type = case_when(
      str_detect(phenotype, "^SubCor_Vol") ~ "SubCorVol",
      TRUE ~ str_split(phenotype_clean, "_") %>% map_chr(1)
    ),
    family = family_map[phenotype_type],
    symbol = toupper(symbol)
  ) %>%
  dplyr::filter(!is.na(family))

flames_genes_all <- unique(flames_raw$symbol)
cat(sprintf("  FLAMES: %d unique genes (all families)\n", length(flames_genes_all)))

# ================================================================
# SECTION 3: Overlap test -- fisher's exact / hypergeometric
# ================================================================
hypergeom_test <- function(query, gene_set, background, label) {
  query    <- intersect(query, background)
  gene_set <- intersect(gene_set, background)
  n_bg <- length(background)
  a <- length(intersect(query, gene_set))
  b <- length(setdiff(query, gene_set))
  c <- length(setdiff(gene_set, query))
  d <- n_bg - a - b - c
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  tibble(group = label, n_query = length(query), n_panel = length(gene_set),
         n_background = n_bg, n_overlap = a, OR = ft$estimate, p_value = ft$p.value,
         overlap_genes = paste(sort(intersect(query, gene_set)), collapse = ";"))
}

cat("\nRunning overlap tests...\n")

# Overall (headline test, reported standalone -- not FDR-corrected
# alongside the per-family tests below)
overall_result <- hypergeom_test(flames_genes_all, panel_genes, background_genes, "All families")
cat(sprintf("  Overall: n_overlap=%d, OR=%.2f, p=%.3e\n",
            overall_result$n_overlap, overall_result$OR, overall_result$p_value))

# Per-family
family_results <- map_dfr(unique(family_map), function(fam) {
  fam_genes <- unique(flames_raw$symbol[flames_raw$family == fam])
  if (length(fam_genes) < 3) return(NULL)
  hypergeom_test(fam_genes, panel_genes, background_genes, fam)
}) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

cat("\nPer-family results:\n")
print(family_results %>% dplyr::select(group, n_query, n_overlap, OR, p_value, fdr))

results_all <- bind_rows(overall_result %>% mutate(fdr = NA_real_), family_results)
write.csv(results_all %>% dplyr::select(-overlap_genes),
          file.path(out_dir, "flames_panel_overlap_results.csv"), row.names = FALSE)
write.csv(results_all %>% dplyr::select(group, overlap_genes),
          file.path(out_dir, "flames_panel_overlap_genes.csv"), row.names = FALSE)

cat("\n=== DONE ===\n")
cat("Outputs in:", out_dir, "\n")
cat("  combined_gene_panel.csv\n")
cat("  flames_panel_overlap_results.csv\n")
cat("  flames_panel_overlap_genes.csv\n")