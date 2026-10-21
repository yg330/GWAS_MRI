# ===========================================================================
# Final_Figure6_arealizationscript.R
# Manuscript item: Figure 5 and Supplementary Tables 10-13 (script is named Figure 6)
#
# What this does:
#   Full arealisation pipeline: spatial clustering of FLAMES hits across
#   parcels, quality filter, enrichment against the Qian and Bhaduri gene
#   sets, cell-type and eRegulon enrichment, developmental-stage
#   specificity, anterior-posterior partition and pathway enrichment.
#   Absorbs the earlier Arealisation_genesscript.R and
#   Arealisation_enrichmentscript.R.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Figure 6: Arealisation genes for cortical measures
# MASTER PIPELINE -- single source script
#
# Consolidates all three prior scripts into one file:
#   PART 1  Spatial clustering of FLAMES gene hits + arealisation
#           candidate identification      (was: Arealisation_genesscript.R)
#   PART 2  Quality filter + DEG enrichment (Qian/Bhaduri)
#   PART 3  Cell type enrichment (Wang 2025 pseudobulk)
#   PART 4  eRegulon target-gene enrichment (pooled)
#   PART 5  Developmental stage specificity
#                                          (was: Arealisation_enrichmentscript.R,
#                                           Sections 1-4/5/6/7 respectively)
#   PART 6  AP spatial partition (frontal/posterior/intermediate)
#           + pathway (GO:BP) enrichment
#   PART 7  AP-specific TF eRegulon enrichment
#                                          (was: AP pathway/TF script)
#   PART 8  Figure 6 assembly (panels A-E) + table exports
#
# All intermediate/heavy computations (spatial permutation tests,
# Fisher tests, GAM/glm enrichment models) are IDENTICAL in logic
# to the three source scripts -- only paths, shared-variable
# naming, and output organisation have been consolidated. Nothing
# in Parts 1-7 has been re-derived from scratch.
#
# OUTPUTS
#   <out_dir>/tables/    every result table (CSV) from all 7 parts
#   <out_dir>/figures/   individual panel PDFs
#   <out_dir>/Figure6_arealisation.pdf / .png   assembled figure
#
# NOTE: Part 1's spatial permutation test (N_PERM=1000 per
# gene x phenotype combination) is the single most expensive step
# in this script, exactly as in the original -- consolidating
# into one file hasn't changed its cost.

source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readxl)
  library(msigdbr)
  library(ggplot2)
  library(cowplot)
  library(ggseg)
  library(ggsegGlasser)
  library(sf)
  library(igraph)
  library(scales)
})

ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# ---- Paths (shared across all parts) ----
base_dir      <- proj_dir
flames_path   <- file.path(proj_dir, "Data", "FLAMES_gene_prioritisation_results.txt")
gene_loc_path <- file.path(base_dir, "Data", "NCBI37.3.gene.loc")
magma_loc     <- gene_loc_path  # same file, used under both names in the source scripts
qian_path     <- file.path(base_dir, "Data", "Supplementary_Table_6_Anterior-posterior DEGs.xlsx")
bhaduri_path  <- file.path(base_dir, "Data", "Bhaduri_arealisation_DEG.csv")
eregulon_xlsx <- file.path(base_dir, "Data", "Wang_eregulonlist.xlsx")

cache_dir     <- file.path(base_dir, "Outputs", "cache")
out_dir     <- file.path(base_dir, "Outputs", "Figure6_arealisation")
tables_dir  <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Global typography (matches Figures 4/5 convention) ----
BASE_SIZE  <- 10
TITLE_SIZE <- 10
LABEL_SIZE <- 9

# ---- Shared theme (matches the theme_bw() convention used for ----
# stats/dotplot panels in Figures 4/5, wired to the same BASE_SIZE/
# TITLE_SIZE/LABEL_SIZE constants above -- mirrors Figure 2's
# theme_fig2() pattern so this figure's panels stop each picking a
# different base theme (theme_void/theme_bw/theme_minimal) ad hoc).
# base_size can be overridden per-panel (e.g. smaller for compact
# panels) while axis/legend text still scales off LABEL_SIZE.
theme_fig6 <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(
      axis.text    = element_text(size = LABEL_SIZE),
      axis.title   = element_text(size = TITLE_SIZE),
      legend.text  = element_text(size = LABEL_SIZE),
      legend.title = element_text(size = LABEL_SIZE),
      strip.text   = element_text(size = LABEL_SIZE)
    )
}

# ---- Family map + colours (matches Figures 4/5) ----
family_map <- c(
  SA = "Brain size", FI = "Brain size", GMV = "Brain size", IC = "Brain size",
  SubCorVol = "Brain size",
  CT = "Cortical thickness",
  GC = "Curvature", MC = "Curvature",
  FA = "Microstructural coherence", ICVF = "Microstructural coherence",
  MD = "Microstructural diffusivity", ISOVF = "Microstructural diffusivity",
  OD = "Orientation dispersion"
)

family_colours <- c(
  "Brain size"                  = "#00C28A",
  "Cortical thickness"          = "#0072B2",
  "Curvature"                   = "#E69F00",
  "Microstructural coherence"   = "#CC79A7",
  "Microstructural diffusivity" = "#56B4E9",
  "Orientation dispersion"      = "#D55E00"
)

# Per-measure colours (matches the phenotype_colours convention used
# in the Figure 6 TF/eRegulon script) -- used for Panel A instead of
# the outlier-group colouring, which has been dropped from the
# manuscript.
phenotype_colours <- c(
  CT    = "#0072B2", FA    = "#C0006E", FI    = "#7DDCB8",
  GC    = "#6B4800", GMV   = "#003D28", IC    = "#2DC48E",
  ICVF  = "#E8A8CB", ISOVF = "#A8D8F0", MC    = "#F5C842",
  MD    = "#0D3F5F", OD    = "#D55E00", SA    = "#007A52"
)

group_order <- c("First_trimester", "Second_trimester",
                 "Third_trimester", "Infancy", "Adolescence")

WINDOW_BP   <- 500000
MIN_PARCELS <- 9

# ---- Outlier parcel groups (from Figure 3) ----
outlier_group1 <- c(
  "RSC","PHA1","PHA2","PHA3","PeEc","EC","PreS","H",
  "TF","TE1a","TE1m","TE1p","TE2a","TGv","TGd",
  "STSvp","STSva","STSdp","STSda","TE2p","PHT",
  "TPOJ1","STV","PSL","RI","52","43","OP4","OP1",
  "OP2-3","PFcm","PoI1","PoI2","Ig","FOP4","MI",
  "FOP1","NPC","AVI","AAIC","PI","Pir","OFC") %>% unique()

outlier_group2 <- c(
  "FOP5","FOP3","FOP2","PFop","AIP","VIP","MIP",
  "LIPd","LIPv","VVC","PHA3","IP0","IP1","IP2",
  "IFJp","IFJa","IFSp","IFSa","p9-46v","46","a9-46v",
  "9-46d","p47r","a47r","13l","47m","47l","11l",
  "OFC","3b","3a","2","1","4","6mp","6d","23c",
  "5mv","5m","PEF","8C") %>% unique()

# ================================================================
# PART 1: Spatial clustering of FLAMES gene hits
# ================================================================
ts("=== PART 1: Spatial clustering ===")

cat("Building Glasser parcel adjacency matrix...\n")
glasser_atlas <- ggsegGlasser::glasser()  # explicit :: -- ggseg masks this with its own glasser()

lh_sf <- glasser_atlas$data$sf %>%
  dplyr::filter(grepl("^lh_L_", label)) %>%
  mutate(parcel = gsub("^lh_L_(.*)_ROI$", "\\1", label))

all_parcels <- unique(lh_sf$parcel)
cat(sprintf("  %d left-hemisphere Glasser parcels found\n", length(all_parcels)))

suppressWarnings({
  adj_mat <- sf::st_touches(lh_sf$geometry, sparse = FALSE)
})
rownames(adj_mat) <- lh_sf$parcel
colnames(adj_mat) <- lh_sf$parcel

lh_sf_lat <- glasser_atlas$data$sf %>%
  dplyr::filter(grepl("^lh_L_", label), view == "lateral") %>%
  mutate(parcel = gsub("^lh_L_(.*)_ROI$", "\\1", label))
lh_sf_med <- glasser_atlas$data$sf %>%
  dplyr::filter(grepl("^lh_L_", label), view == "medial") %>%
  mutate(parcel = gsub("^lh_L_(.*)_ROI$", "\\1", label))

# Bilateral versions for plotting -- the Glasser parcellation here is
# bilaterally averaged, so a gene's significant parcels apply
# identically to both hemispheres. Matches on the hemisphere-stripped
# parcel name regardless of "lh_L_"/"rh_R_" prefix, so both sides
# highlight together when a parcel is in hit_p.
full_sf_lat <- glasser_atlas$data$sf %>%
  dplyr::filter(view == "lateral", grepl("^(lh|rh)_[LR]_.*_ROI$", label)) %>%
  mutate(parcel = gsub("^(lh|rh)_[LR]_(.*)_ROI$", "\\2", label))
full_sf_med <- glasser_atlas$data$sf %>%
  dplyr::filter(view == "medial", grepl("^(lh|rh)_[LR]_.*_ROI$", label)) %>%
  mutate(parcel = gsub("^(lh|rh)_[LR]_(.*)_ROI$", "\\2", label))

cat("Loading FLAMES gene prioritisation results...\n")
flames_raw <- read.table(flames_path, header = TRUE, stringsAsFactors = FALSE)

flames <- flames_raw %>%
  dplyr::filter(grepl("^regional_", phenotype)) %>%
  dplyr::filter(!grepl("SubCorVol|SubCor", phenotype)) %>%
  mutate(
    phenotype_clean = str_remove(phenotype, "^regional_"),
    phenotype_type  = str_split(phenotype_clean, "_") %>% map_chr(1),
    parcel          = str_replace(phenotype_clean, paste0("^", phenotype_type, "_"), "")
  ) %>%
  dplyr::filter(parcel %in% all_parcels)

cat(sprintf("  %d FLAMES hits across %d regional phenotypes\n",
            nrow(flames), n_distinct(flames$phenotype)))

cat("Running spatial clustering tests...\n")
set.seed(42)
N_PERM <- 1000

gene_pheno_hits <- flames %>%
  group_by(symbol, phenotype_type) %>%
  summarise(hit_parcels = list(unique(parcel)), .groups = "drop") %>%
  mutate(n_hits = map_int(hit_parcels, length)) %>%
  dplyr::filter(n_hits >= MIN_PARCELS)

cat(sprintf("  %d gene x phenotype combinations with >= %d parcels\n",
            nrow(gene_pheno_hits), MIN_PARCELS))

count_adjacent_pairs <- function(parcels, adj) {
  n <- length(parcels)
  if (n < 2) return(0L)
  p <- parcels[parcels %in% rownames(adj)]
  if (length(p) < 2) return(0L)
  sub <- adj[p, p, drop = FALSE]
  sum(sub[upper.tri(sub)])
}

results <- gene_pheno_hits %>%
  mutate(
    obs_adj = map_int(hit_parcels, count_adjacent_pairs, adj = adj_mat),
    perm_dist = map2(hit_parcels, n_hits, function(hp, nh) {
      replicate(N_PERM, {
        rand_p <- sample(all_parcels, nh, replace = FALSE)
        count_adjacent_pairs(rand_p, adj_mat)
      })
    }),
    perm_mean = map_dbl(perm_dist, mean),
    perm_sd   = map_dbl(perm_dist, sd),
    spatial_score = obs_adj / pmax(perm_mean, 0.5),
    perm_p = map2_dbl(obs_adj, perm_dist, function(obs, pd) mean(pd >= obs))
  ) %>%
  dplyr::select(-perm_dist)

results$fdr <- p.adjust(results$perm_p, method = "BH")

sig_clusters <- results %>%
  dplyr::filter(fdr < 0.05) %>%
  arrange(fdr, desc(spatial_score)) %>%
  mutate(
    n_g1 = map_int(hit_parcels, ~ sum(.x %in% outlier_group1)),
    n_g2 = map_int(hit_parcels, ~ sum(.x %in% outlier_group2)),
    outlier_overlap = case_when(
      n_g1 > 0 & n_g2 > 0 ~ "Both groups",
      n_g1 > 0             ~ "Paralimbic only",
      n_g2 > 0             ~ "Insula/SM only",
      TRUE                 ~ "None"
    )
  )

cat(sprintf("  %d FDR-significant spatial clusters (FDR < 0.05)\n", nrow(sig_clusters)))

compute_coherence <- function(hit_p, adj, obs_adj_val) {
  p <- hit_p[hit_p %in% rownames(adj)]
  n <- length(p)
  if (n < 2) return(list(n_components = 1L, LCC_fraction = 1, compactness = NA_real_))
  sub <- adj[p, p, drop = FALSE]
  g   <- igraph::graph_from_adjacency_matrix(sub, mode = "undirected", diag = FALSE)
  comps <- igraph::components(g)
  n_comp <- comps$no
  lcc_frac <- max(comps$csize) / n
  min_edges <- n - n_comp
  compactness <- if (min_edges > 0) obs_adj_val / min_edges else NA_real_
  list(n_components = n_comp, LCC_fraction = lcc_frac, compactness = compactness)
}

sig_clusters <- sig_clusters %>%
  mutate(
    coherence    = map2(hit_parcels, obs_adj, compute_coherence, adj = adj_mat),
    n_components = map_int(coherence, "n_components"),
    LCC_fraction = map_dbl(coherence, "LCC_fraction"),
    compactness  = map_dbl(coherence, "compactness")
  ) %>%
  dplyr::select(-coherence)

write.csv(results %>% mutate(hit_parcels = map_chr(hit_parcels, paste, collapse = ";")),
          file.path(tables_dir, "spatial_clustering_results.csv"), row.names = FALSE)

# ---- Cross-phenotype sanity-check annotation ----
flames_hits <- flames %>% dplyr::select(symbol, phenotype_type, parcel) %>% distinct()

sig_lookup <- sig_clusters %>%
  dplyr::select(symbol, phenotype_type, hit_parcels) %>%
  mutate(family = family_map[phenotype_type])

jaccard <- function(a, b) {
  a <- unlist(a); b <- unlist(b)
  if (length(union(a, b)) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

sig_clusters_annotated <- sig_clusters %>%
  mutate(family = family_map[phenotype_type]) %>%
  mutate(
    same_fam_entries = pmap(list(symbol, phenotype_type, family), function(sym, phen, fam) {
      sig_lookup %>% dplyr::filter(symbol == sym, family == fam, phenotype_type != phen)
    }),
    n_same_family_sig = map_int(same_fam_entries, nrow),
    same_family_sig_phenos = map_chr(same_fam_entries, function(e)
      if (nrow(e) == 0) "" else paste(sort(unique(e$phenotype_type)), collapse = ";")),
    same_family_jaccard = map2_chr(hit_parcels, same_fam_entries, function(hp, entries) {
      if (nrow(entries) == 0) return("")
      entries %>% rowwise() %>% mutate(jac = jaccard(hp, hit_parcels)) %>% ungroup() %>%
        mutate(entry = sprintf("%s:%.2f", phenotype_type, jac)) %>%
        pull(entry) %>% paste(collapse = ";")
    }),
    mean_same_family_jaccard = map2_dbl(hit_parcels, same_fam_entries, function(hp, entries) {
      if (nrow(entries) == 0) return(NA_real_)
      mean(map_dbl(entries$hit_parcels, ~ jaccard(hp, .x)), na.rm = TRUE)
    }),
    other_fam_entries = pmap(list(symbol, family), function(sym, fam) {
      sig_lookup %>% dplyr::filter(symbol == sym, family != fam)
    }),
    n_other_family_sig = map_int(other_fam_entries, ~ n_distinct(.x$family)),
    other_family_sig = map_chr(other_fam_entries, function(entries)
      if (nrow(entries) == 0) "" else paste(sort(unique(entries$family)), collapse = ";")),
    other_family_jaccard = map2_chr(hit_parcels, other_fam_entries, function(hp, entries) {
      if (nrow(entries) == 0) return("")
      entries %>% group_by(family) %>%
        summarise(mean_jac = mean(map_dbl(hit_parcels, ~ jaccard(hp, .x)), na.rm = TRUE), .groups = "drop") %>%
        mutate(entry = sprintf("%s:%.2f", family, mean_jac)) %>%
        pull(entry) %>% paste(collapse = ";")
    }),
    overlapping_pheno_parcels = map2_chr(symbol, hit_parcels, function(sym, hp) {
      hp <- unlist(hp)
      other_hits <- flames_hits %>% dplyr::filter(symbol == sym, parcel %in% hp)
      if (nrow(other_hits) == 0) return("")
      other_hits %>% group_by(phenotype_type) %>%
        summarise(parcels = paste(sort(parcel), collapse = ","), .groups = "drop") %>%
        mutate(entry = paste0(phenotype_type, ":", parcels)) %>%
        pull(entry) %>% paste(collapse = ";")
    })
  ) %>%
  dplyr::select(-same_fam_entries, -other_fam_entries)

write.csv(
  sig_clusters_annotated %>% mutate(hit_parcels = map_chr(hit_parcels, paste, collapse = ";")),
  file.path(tables_dir, "spatial_clustering_significant.csv"), row.names = FALSE)

# ---- Arealisation lookup (hg19 coordinates + windows) ----
gene_loc <- read.table(gene_loc_path, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("entrez_id", "chr", "start_hg19", "end_hg19", "strand", "symbol")) %>%
  dplyr::filter(chr %in% c(1:22, "X", "Y")) %>%
  mutate(chr = as.character(chr))

sig_annotated <- sig_clusters_annotated %>%
  mutate(hit_parcels_str = map_chr(hit_parcels, paste, collapse = ";"),
         n_hit_parcels   = map_int(hit_parcels, length)) %>%
  dplyr::select(symbol, phenotype_type, family, n_hit_parcels, hit_parcels_str,
                spatial_score, fdr, outlier_overlap, n_components, LCC_fraction, compactness,
                n_same_family_sig, same_family_sig_phenos, same_family_jaccard, mean_same_family_jaccard,
                n_other_family_sig, other_family_sig, other_family_jaccard, overlapping_pheno_parcels) %>%
  dplyr::rename(hit_parcels = hit_parcels_str) %>%
  left_join(gene_loc %>% dplyr::select(symbol, entrez_id, chr, start_hg19, end_hg19),
            by = "symbol", relationship = "many-to-many") %>%
  mutate(start_window = pmax(1L, start_hg19 - WINDOW_BP), end_window = end_hg19 + WINDOW_BP)

arealisation_lookup <- sig_annotated %>%
  dplyr::filter(!is.na(chr)) %>%
  dplyr::select(symbol, entrez_id, chr, start_hg19, end_hg19, start_window, end_window,
                phenotype_type, family, n_hit_parcels, hit_parcels, spatial_score, fdr, outlier_overlap,
                n_components, LCC_fraction, compactness, n_same_family_sig, same_family_sig_phenos,
                same_family_jaccard, mean_same_family_jaccard, n_other_family_sig, other_family_sig,
                other_family_jaccard, overlapping_pheno_parcels) %>%
  arrange(chr, start_hg19, symbol, phenotype_type)

write.csv(arealisation_lookup, file.path(tables_dir, "arealisation_lookup.csv"), row.names = FALSE)

ts("=== PART 1 done ===")

# ================================================================
# PART 2: Quality filter + DEG enrichment (Qian/Bhaduri)
# ================================================================
ts("=== PART 2: Quality filter + DEG enrichment ===")

sig <- sig_clusters_annotated %>%
  mutate(hit_parcels = map_chr(hit_parcels, paste, collapse = ";"),
         n_hits       = map_int(str_split(hit_parcels, ";"), length),
         LCC_fraction = as.numeric(LCC_fraction),
         n_components = as.integer(n_components),
         compactness  = as.numeric(compactness))

sig_filtered <- sig %>%
  dplyr::filter(n_hits <= 60 & LCC_fraction >= 0.30 & n_components <= 8)

# ---- Beta-level QC exclusion ----
# WNT3/SA, RUNX2/MD, KANSL1/GMV are spatially-clustered FLAMES hits
# that passed the filters above but did NOT hold up under lead-SNP
# beta validation (Arealisation_beta_hitparcel_validation.R,
# PC1_vs_arealisation_beta_correlation.R): |beta| inside the gene's
# own hit_parcels was not distinguishable from outside, and/or no
# distance-decay away from the hit set -- i.e. not bona fide
# arealisation signal, just an artefact of the spatial-clustering
# test. Excluded here (rather than only post-hoc in the CSV) so the
# exclusion propagates through every downstream part of this
# pipeline (DEG enrichment, cell-type enrichment, eRegulon
# enrichment, AP partition/pathway enrichment, Figure 6 panels).
beta_qc_exclude <- tibble::tribble(
  ~symbol,   ~phenotype_type,
  "WNT3",    "SA",
  "RUNX2",   "MD",
  "KANSL1",  "GMV"
)

n_before_betaqc <- nrow(sig_filtered)
sig_filtered <- sig_filtered %>%
  dplyr::anti_join(beta_qc_exclude, by = c("symbol", "phenotype_type"))

cat(sprintf("  Before filter: %d entries, %d unique genes\n", nrow(sig), n_distinct(sig$symbol)))
cat(sprintf("  After spatial-clustering filter: %d entries, %d unique genes\n",
            n_before_betaqc, n_distinct(sig$symbol[sig$n_hits <= 60 & sig$LCC_fraction >= 0.30 & sig$n_components <= 8])))
cat(sprintf("  After beta-QC exclusion (WNT3/SA, RUNX2/MD, KANSL1/GMV removed): %d entries, %d unique genes\n",
            nrow(sig_filtered), n_distinct(sig_filtered$symbol)))

write.csv(sig_filtered, file.path(tables_dir, "arealisation_filtered.csv"), row.names = FALSE)

arealisation_genes_raw <- unique(sig_filtered$symbol)
arealisation_by_family_raw <- purrr::map(unique(family_map), function(fam) {
  phenos <- names(family_map)[family_map == fam]
  unique(sig_filtered$symbol[sig_filtered$phenotype_type %in% phenos])
}) %>% setNames(unique(family_map))

background_genes   <- unique(gene_loc$symbol)
arealisation_genes <- intersect(arealisation_genes_raw, background_genes)

cat(sprintf("  MAGMA background: %d genes | Arealisation genes in background: %d\n",
            length(background_genes), length(arealisation_genes)))

qian_sheets <- readxl::excel_sheets(qian_path)
qian_list <- purrr::map(qian_sheets, function(sh) {
  readxl::read_excel(qian_path, sheet = sh) %>%
    dplyr::filter(pvals_adj < 0.05) %>%
    mutate(direction = ifelse(log_fold_change > 0, "anterior", "posterior"))
}) %>% setNames(qian_sheets)

qian_sets <- purrr::map(qian_sheets, function(sh) {
  list(anterior = dplyr::filter(qian_list[[sh]], direction == "anterior")$gene,
       posterior = dplyr::filter(qian_list[[sh]], direction == "posterior")$gene,
       all = qian_list[[sh]]$gene)
}) %>% setNames(qian_sheets)

bhaduri_ap <- c(pfc = "anterior", motor = "middle_anterior", somatosensory = "middle",
                parietal = "middle", temporal = "middle_posterior", v1 = "posterior")

bhaduri <- read.csv(bhaduri_path, stringsAsFactors = FALSE) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  mutate(ap_label = bhaduri_ap[area])

bhaduri_sets <- list(
  by_area = split(bhaduri$gene, bhaduri$area),
  anterior = unique(bhaduri$gene[bhaduri$ap_label == "anterior"]),
  posterior = unique(bhaduri$gene[bhaduri$ap_label == "posterior"]),
  all = unique(bhaduri$gene)
)

qian_all_genes      <- unique(unlist(purrr::map(qian_sets, ~ .x$all)))
pooled_arealisation <- unique(c(qian_all_genes, bhaduri_sets$all))

cat(sprintf("  Pooled: %d Qian + %d Bhaduri = %d unique\n",
            length(qian_all_genes), length(bhaduri_sets$all), length(pooled_arealisation)))

fisher_test <- function(query, gene_set, background) {
  query    <- intersect(query,    background)
  gene_set <- intersect(gene_set, background)
  n_bg <- length(background)
  a <- length(intersect(query, gene_set)); b <- length(setdiff(query, gene_set))
  c <- length(setdiff(gene_set, query));   d <- n_bg - a - b - c
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  data.frame(n_query = length(query), n_gene_set = length(gene_set), n_background = n_bg,
             n_overlap = a, overlap_genes = paste(sort(intersect(query, gene_set)), collapse = ";"),
             OR = ft$estimate, p_value = ft$p.value, stringsAsFactors = FALSE)
}

enrich_results <- list()

res <- fisher_test(arealisation_genes, pooled_arealisation, background_genes)
enrich_results[[1]] <- data.frame(level = "primary", source = "Qian+Bhaduri", gene_set = "pooled_all",
                                  cell_type = "all", direction = "all", family = "all", res, stringsAsFactors = FALSE)

for (fam in unique(family_map)) {
  phenos <- names(family_map)[family_map == fam]
  fam_genes <- intersect(unique(sig_filtered$symbol[sig_filtered$phenotype_type %in% phenos]), background_genes)
  if (length(fam_genes) < 3) next
  res <- fisher_test(fam_genes, pooled_arealisation, background_genes)
  enrich_results[[length(enrich_results) + 1]] <- data.frame(
    level = "family", source = "Qian+Bhaduri", gene_set = "pooled_all", cell_type = "all",
    direction = "all", family = fam, res, stringsAsFactors = FALSE)
}

for (sh in qian_sheets) {
  gs <- intersect(qian_sets[[sh]]$all, background_genes)
  if (length(gs) < 5) next
  res <- fisher_test(arealisation_genes, gs, background_genes)
  enrich_results[[length(enrich_results) + 1]] <- data.frame(
    level = "cell_type", source = "Qian", gene_set = sh, cell_type = sh,
    direction = "all", family = "all", res, stringsAsFactors = FALSE)
}

for (area_i in names(bhaduri_sets$by_area)) {
  gs <- intersect(bhaduri_sets$by_area[[area_i]], background_genes)
  if (length(gs) < 5) next
  res <- fisher_test(arealisation_genes, gs, background_genes)
  enrich_results[[length(enrich_results) + 1]] <- data.frame(
    level = "area", source = "Bhaduri", gene_set = paste0("area_", area_i), cell_type = "all",
    direction = bhaduri_ap[area_i], family = "all", res, stringsAsFactors = FALSE)
}

# -- Overall (pooled) Qian and Bhaduri tests, INDEPENDENTLY --
# (new for Panel C: single overall test per source, not split by
# subtype/area, and not pooled together the way "primary" is)
res_qian_overall    <- fisher_test(arealisation_genes, qian_all_genes, background_genes)
res_bhaduri_overall <- fisher_test(arealisation_genes, bhaduri_sets$all, background_genes)

enrich_results[[length(enrich_results) + 1]] <- data.frame(
  level = "overall", source = "Qian", gene_set = "qian_overall", cell_type = "all",
  direction = "all", family = "all", res_qian_overall, stringsAsFactors = FALSE)
enrich_results[[length(enrich_results) + 1]] <- data.frame(
  level = "overall", source = "Bhaduri", gene_set = "bhaduri_overall", cell_type = "all",
  direction = "all", family = "all", res_bhaduri_overall, stringsAsFactors = FALSE)

enrich_df <- bind_rows(enrich_results) %>%
  group_by(level) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(level, p_value)

write.csv(enrich_df, file.path(tables_dir, "enrichment_results.csv"), row.names = FALSE)
cat(sprintf("  %d tests; %d FDR < 0.05\n", nrow(enrich_df), sum(enrich_df$fdr < 0.05, na.rm = TRUE)))

ts("=== PART 2 done ===")

# ================================================================
# PART 3: Cell type enrichment (Wang 2025 pseudobulk)
# ================================================================
ts("=== PART 3: Cell type enrichment ===")

pseudobulk_pc   <- readRDS(file.path(cache_dir, "pseudobulk_protein_coding.rds"))
pc_pool         <- rownames(pseudobulk_pc)
gene_lengths_df <- readRDS(file.path(cache_dir, "gene_lengths_hg19.rds"))
mean_expr_all   <- rowMeans(pseudobulk_pc, na.rm = TRUE)

ar_in_pb <- intersect(arealisation_genes_raw, pc_pool)
arealisation_by_family <- purrr::map(arealisation_by_family_raw, ~ intersect(.x, pc_pool)) %>%
  setNames(names(arealisation_by_family_raw))

gene_covars_ar <- data.frame(symbol = pc_pool, stringsAsFactors = FALSE) %>%
  left_join(gene_lengths_df %>% dplyr::select(symbol, gene_length), by = "symbol") %>%
  mutate(log_length = log10(gene_length + 1))

bin_meta <- data.frame(bin = colnames(pseudobulk_pc), stringsAsFactors = FALSE) %>%
  mutate(cell_type = gsub("__.*", "", bin),
         stage     = factor(gsub(".*__", "", bin), levels = group_order))

logreg_enrich <- function(query_genes, predictor_col, data_df) {
  df <- data_df %>% dplyr::filter(!is.na(log_length), !is.na(.data[[predictor_col]])) %>%
    mutate(is_query = as.integer(symbol %in% query_genes))
  if (sum(df$is_query) < 3) return(NULL)
  if (var(df[[predictor_col]]) == 0) return(NULL)
  tryCatch({
    fm  <- as.formula(paste("is_query ~ log_length +", predictor_col))
    fit <- glm(fm, data = df, family = binomial())
    co  <- summary(fit)$coefficients
    if (!predictor_col %in% rownames(co)) return(NULL)
    row <- co[predictor_col, , drop = FALSE]
    data.frame(OR = exp(row[, "Estimate"]), p_value = row[, "Pr(>|z|)"],
               n_overlap = sum(df$is_query == 1 & df[[predictor_col]] > 0),
               n_query = sum(df$is_query), stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}

run_ct_enrichment <- function(query_genes, label) {
  query_genes <- intersect(query_genes, pc_pool)
  if (length(query_genes) < 5) return(NULL)
  map_dfr(colnames(pseudobulk_pc), function(bn) {
    expr_bin <- pseudobulk_pc[, bn, drop = TRUE]
    log_rel_expr <- log10(expr_bin / (mean_expr_all + 0.01) + 0.01)
    df <- gene_covars_ar %>% mutate(log_rel_expr = log_rel_expr[symbol])
    res <- logreg_enrich(query_genes, "log_rel_expr", df)
    if (is.null(res)) return(NULL)
    data.frame(gene_set = label, bin = bn, cell_type = bin_meta$cell_type[bin_meta$bin == bn],
               stage = as.character(bin_meta$stage[bin_meta$bin == bn]), res, stringsAsFactors = FALSE)
  })
}

ct_results <- list()
ct_results[["All arealisation"]] <- run_ct_enrichment(ar_in_pb, "All arealisation")
for (fam in names(arealisation_by_family)) {
  gs <- arealisation_by_family[[fam]]
  if (length(gs) < 5) next
  ct_results[[fam]] <- run_ct_enrichment(gs, fam)
}

ct_df <- bind_rows(ct_results) %>%
  group_by(gene_set) %>% mutate(fdr = p.adjust(p_value, method = "BH")) %>% ungroup() %>%
  mutate(stage = factor(stage, levels = group_order)) %>% arrange(gene_set, p_value)

write.csv(ct_df, file.path(tables_dir, "arealisation_celltype_enrichment.csv"), row.names = FALSE)
cat(sprintf("  %d tests; %d FDR < 0.05\n", nrow(ct_df), sum(ct_df$fdr < 0.05, na.rm = TRUE)))

ts("=== PART 3 done ===")

# ================================================================
# PART 4: eRegulon target-gene enrichment (pooled)
# ================================================================
ts("=== PART 4: eRegulon enrichment (pooled) ===")

ereg_df <- NULL

if (!file.exists(eregulon_xlsx)) {
  cat("  WARNING: Wang_eregulonlist.xlsx not found -- skipping Part 4\n")
} else {
  sheets_er <- readxl::excel_sheets(eregulon_xlsx)
  eregulon_raw <- map_dfr(sheets_er, function(sh) {
    df <- readxl::read_excel(eregulon_xlsx, sheet = sh,
                             col_names = c("eregulon", "eregulon_n", "target_genes"), skip = 1)
    df %>% dplyr::filter(!is.na(eregulon), !is.na(target_genes)) %>%
      mutate(tf = str_remove(eregulon, "_extended.*|_\\+.*|_-.*"),
             targets = str_split(target_genes, ",")) %>%
      dplyr::select(eregulon, tf, targets)
  })
  
  eregulon_list <- eregulon_raw %>%
    unnest(targets) %>% mutate(target = str_trim(targets)) %>%
    group_by(tf) %>% summarise(target_vec = list(unique(target)), .groups = "drop") %>%
    { setNames(.$target_vec, .$tf) }
  
  if (length(eregulon_list) > 0) {
    gene_covars_er <- data.frame(symbol = pc_pool, stringsAsFactors = FALSE) %>%
      left_join(gene_lengths_df %>% dplyr::select(symbol, gene_length), by = "symbol") %>%
      mutate(log_length = log10(gene_length + 1))
    
    run_ereg_enrichment <- function(query_genes, label) {
      query_genes <- intersect(query_genes, pc_pool)
      if (length(query_genes) < 5) return(NULL)
      map_dfr(names(eregulon_list), function(tf) {
        targets <- intersect(eregulon_list[[tf]], pc_pool)
        if (length(targets) < 5) return(NULL)
        df <- gene_covars_er %>% mutate(is_target = as.integer(symbol %in% targets))
        res <- logreg_enrich(query_genes, "is_target", df)
        if (is.null(res)) return(NULL)
        data.frame(gene_set = label, tf = tf, n_targets = length(targets), res, stringsAsFactors = FALSE)
      })
    }
    
    ereg_results <- list()
    ereg_results[["All arealisation"]] <- run_ereg_enrichment(ar_in_pb, "All arealisation")
    for (fam in names(arealisation_by_family)) {
      gs <- arealisation_by_family[[fam]]
      if (length(gs) < 5) next
      ereg_results[[fam]] <- run_ereg_enrichment(gs, fam)
    }
    
    ereg_df <- bind_rows(ereg_results) %>%
      group_by(gene_set) %>% mutate(fdr = p.adjust(p_value, method = "BH")) %>% ungroup() %>%
      mutate(tf_is_arealisation = tf %in% arealisation_genes_raw) %>%
      arrange(gene_set, p_value)
    
    write.csv(ereg_df, file.path(tables_dir, "arealisation_eregulon_enrichment.csv"), row.names = FALSE)
    cat(sprintf("  %d tests; %d FDR < 0.05\n", nrow(ereg_df), sum(ereg_df$fdr < 0.05, na.rm = TRUE)))
  }
}

ts("=== PART 4 done ===")

# ================================================================
# PART 5: Developmental stage specificity
# ================================================================
ts("=== PART 5: Developmental stage specificity ===")

set.seed(42)
N_PERM_STAGE <- 1000
n_ar <- length(ar_in_pb)

if (exists("ct_df") && any(ct_df$fdr < 0.05 & ct_df$OR > 1, na.rm = TRUE)) {
  top_cts <- ct_df %>% dplyr::filter(gene_set == "All arealisation", fdr < 0.05, OR > 1) %>%
    count(cell_type, sort = TRUE) %>% slice_head(n = 5) %>% pull(cell_type)
} else {
  top_cts <- unique(bin_meta$cell_type)[1:5]
}

stage_results <- map_dfr(top_cts, function(ct) {
  bins_ct <- bin_meta$bin[bin_meta$cell_type == ct]
  bins_ct <- bins_ct[bins_ct %in% colnames(pseudobulk_pc)]
  if (length(bins_ct) == 0) return(NULL)
  map_dfr(bins_ct, function(bn) {
    expr_bin <- pseudobulk_pc[, bn, drop = TRUE]
    log_rel_expr <- log10(expr_bin / (mean_expr_all + 0.01) + 0.01)
    obs_mean <- mean(log_rel_expr[ar_in_pb], na.rm = TRUE)
    perm_means <- replicate(N_PERM_STAGE, mean(log_rel_expr[sample(pc_pool, n_ar, replace = FALSE)], na.rm = TRUE))
    data.frame(cell_type = ct, bin = bn, stage = as.character(bin_meta$stage[bin_meta$bin == bn]),
               obs_mean = obs_mean, perm_mean = mean(perm_means), perm_sd = sd(perm_means),
               perm_p = mean(perm_means >= obs_mean), stringsAsFactors = FALSE)
  })
})

stage_results <- stage_results %>%
  mutate(stage = factor(stage, levels = group_order), fdr = p.adjust(perm_p, method = "BH"),
         z = (obs_mean - perm_mean) / pmax(perm_sd, 1e-6))

write.csv(stage_results, file.path(tables_dir, "arealisation_stage_specificity.csv"), row.names = FALSE)
cat(sprintf("  %d bin tests; %d FDR < 0.05\n", nrow(stage_results), sum(stage_results$fdr < 0.05, na.rm = TRUE)))

ts("=== PART 5 done ===")

# ================================================================
# PART 6: AP spatial partition + pathway (GO:BP) enrichment
# ================================================================
ts("=== PART 6: AP pathway enrichment ===")

frontal_parcels <- c(
  "8Ad","8Av","8BL","8BM","8C","9-46d","p9-46v","a9-46v","46",
  "9a","9m","9p","10r","10v","10pp","10d","IFSa","IFSp","IFJa","IFJp",
  "44","45","47l","47m","47s","a47r","p47r","FEF","55b","6r","6v",
  "SCEF","SFL","a32pr","d32","p32","p32pr","24dd","24dv","s32","6ma","6mp")

posterior_parcels <- c(
  "V1","V2","V3","V4","V3A","V3B","V3CD","V4t","V6","V6A","V7","V8",
  "MT","MST","LO1","LO2","LO3","PIT","FFC","VMV1","VMV2","VMV3",
  "VVC","FST","PH","PHA1","PHA2","PHA3","RSC","POS1","POS2",
  "IP0","IP1","IP2","AIP","VIP","MIP","LIPd","LIPv","7Am","7AL","7m","7PC","7Pm")

exclude_ap <- c("AC021218.2", "NPM1", "C16orf95")

gene_parcels <- sig_filtered %>%
  dplyr::filter(!symbol %in% exclude_ap) %>%
  mutate(parcel_list = strsplit(hit_parcels, ";")) %>%
  tidyr::unnest(parcel_list) %>%
  group_by(symbol) %>%
  summarise(all_parcels = list(unique(parcel_list)), n_parcels = n_distinct(parcel_list),
            n_frontal = sum(unique(parcel_list) %in% frontal_parcels),
            frontal_frac = n_frontal / n_distinct(parcel_list), .groups = "drop") %>%
  mutate(n_posterior = map_int(all_parcels, ~ sum(.x %in% posterior_parcels)),
         posterior_frac = n_posterior / n_parcels)

go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

specific_terms <- c(
  "GOBP_FIBROBLAST_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_EPIDERMAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_ERK1_AND_ERK2_CASCADE",
  "GOBP_RETINOIC_ACID_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_RETINOIC_ACID_METABOLIC_PROCESS",
  "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
  "GOBP_BMP_SIGNALING_PATHWAY",
  "GOBP_NEGATIVE_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY",
  "GOBP_NOTCH_SIGNALING_PATHWAY")

found_terms   <- intersect(specific_terms, unique(go_bp$gs_name))
missing_terms <- setdiff(specific_terms, found_terms)
if (length(missing_terms) > 0) {
  cat("  WARNING: requested GO:BP terms not found:\n")
  cat(paste(" -", missing_terms), sep = "\n")
}

candidate_sets <- go_bp %>% dplyr::filter(gs_name %in% found_terms) %>%
  split(x = .$gene_symbol, f = .$gs_name)

if (length(candidate_sets) == 0) {
  stop("No candidate GO:BP gene sets matched -- check the WARNING block above.")
}

fisher_binary <- function(query_genes, gene_set, background) {
  gene_set <- intersect(gene_set, background)
  if (length(gene_set) < 2) return(NULL)
  a <- sum(query_genes %in% gene_set); b <- length(query_genes) - a
  c <- length(gene_set) - a; d <- length(background) - a - b - c
  if (a == 0) return(NULL)
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  data.frame(n_query = length(query_genes), n_background = length(background), n_set = length(gene_set),
             n_overlap = a, overlap_genes = paste(sort(query_genes[query_genes %in% gene_set]), collapse = ";"),
             OR = ft$estimate, p_value = ft$p.value, stringsAsFactors = FALSE)
}

thresholds  <- c(0.30, 0.40, 0.50)
gene_groups <- c("frontal", "posterior", "intermediate")
all_results_ap <- list()

for (thresh in thresholds) {
  frontal      <- gene_parcels$symbol[gene_parcels$frontal_frac   >  thresh]
  posterior    <- gene_parcels$symbol[gene_parcels$posterior_frac >  thresh]
  intermediate <- gene_parcels$symbol[gene_parcels$frontal_frac <= thresh & gene_parcels$posterior_frac <= thresh]
  
  for (grp in gene_groups) {
    query_genes <- switch(grp, frontal = frontal, posterior = posterior, intermediate = intermediate)
    if (length(query_genes) < 3) next
    res <- map_dfr(names(candidate_sets), function(nm) {
      r <- fisher_binary(query_genes, candidate_sets[[nm]], background_genes)
      if (is.null(r)) return(NULL)
      data.frame(pathway = nm, group = grp, threshold = thresh, r, stringsAsFactors = FALSE)
    })
    if (nrow(res) == 0) next
    res <- res %>% mutate(fdr = p.adjust(p_value, method = "BH")) %>% arrange(p_value)
    all_results_ap[[paste0(grp, "_", thresh)]] <- res
  }
}

combined <- bind_rows(all_results_ap) %>%
  mutate(threshold_label = paste0(">", threshold * 100, "%"),
         group = factor(group, levels = c("frontal", "posterior", "intermediate")))

write.csv(combined, file.path(tables_dir, "arealisation_binary_gsea.csv"), row.names = FALSE)
cat(sprintf("  %d tests; %d p<0.05\n", nrow(combined), sum(combined$p_value < 0.05, na.rm = TRUE)))

ts("=== PART 6 done ===")

# ================================================================
# PART 7: AP-specific TF eRegulon enrichment
# ================================================================
ts("=== PART 7: AP TF eRegulon enrichment ===")

tfs_of_interest <- c(
  "PAX6", "SP8", "EMX2", "NR2F1", "MEIS2",
  "TCF7L2", "MEIS1", "BCL11B", "PRDM16",
  "RORA", "ZIC1", "ZEB1", "SOX5", "LHX2", "FOXG1")

sheets_er2 <- readxl::excel_sheets(eregulon_xlsx)
eregulon_raw2 <- map_dfr(sheets_er2, function(sh) {
  df <- readxl::read_excel(eregulon_xlsx, sheet = sh,
                           col_names = c("eregulon", "eregulon_n", "target_genes"), skip = 1)
  df %>% dplyr::filter(!is.na(eregulon), !is.na(target_genes)) %>%
    mutate(tf = str_remove(eregulon, "_extended.*|_\\+.*|_-.*"),
           mode = case_when(
             str_detect(eregulon, "_\\+_[+-]$") ~ "activator",
             str_detect(eregulon, "_-_[+-]$")   ~ "repressor",
             TRUE                                ~ NA_character_
           ),
           targets = str_split(target_genes, ",")) %>%
    dplyr::select(eregulon, tf, mode, targets)
})

tf_targets <- eregulon_raw2 %>%
  dplyr::filter(tf %in% tfs_of_interest, !is.na(mode)) %>%
  unnest(targets) %>% mutate(target = str_trim(targets)) %>%
  group_by(tf, mode) %>%
  summarise(targets = list(unique(target)), n_targets = n_distinct(target), .groups = "drop") %>%
  arrange(tf, mode)

tf_all_results <- list()

if (nrow(tf_targets) > 0) {
  for (thresh in thresholds) {
    frontal      <- gene_parcels$symbol[gene_parcels$frontal_frac   >  thresh]
    posterior    <- gene_parcels$symbol[gene_parcels$posterior_frac >  thresh]
    intermediate <- gene_parcels$symbol[gene_parcels$frontal_frac <= thresh & gene_parcels$posterior_frac <= thresh]
    
    for (grp in gene_groups) {
      query_genes <- switch(grp, frontal = frontal, posterior = posterior, intermediate = intermediate)
      if (length(query_genes) < 3) next
      res <- map_dfr(seq_len(nrow(tf_targets)), function(i) {
        r <- fisher_binary(query_genes, tf_targets$targets[[i]], background_genes)
        if (is.null(r)) return(NULL)
        data.frame(tf = tf_targets$tf[i], mode = tf_targets$mode[i], group = grp, threshold = thresh,
                   r, stringsAsFactors = FALSE)
      })
      if (nrow(res) == 0) next
      res <- res %>% mutate(fdr = p.adjust(p_value, method = "BH")) %>% arrange(p_value)
      tf_all_results[[paste0(grp, "_", thresh)]] <- res
    }
  }
}

tf_combined <- bind_rows(tf_all_results) %>%
  mutate(threshold_label = paste0(">", threshold * 100, "%"),
         group = factor(group, levels = c("frontal", "posterior", "intermediate")))

write.csv(tf_combined, file.path(tables_dir, "arealisation_TF_enrichment_AP.csv"), row.names = FALSE)
cat(sprintf("  %d tests; %d p<0.05\n", nrow(tf_combined), sum(tf_combined$p_value < 0.05, na.rm = TRUE)))

ts("=== PART 7 done ===")

# ================================================================
# PART 8: FIGURE 6 -- panels A-E
# ================================================================
ts("=== PART 8: Figure assembly ===")

# ---- Panel A: exemplar genes, BOTH views, coloured by measure --
exemplar_genes <- c("FOXP1", "BCL11B", "MEIS1", "THBS1", "SPRY2", "TCF7L2")
exemplar_present <- intersect(exemplar_genes, unique(sig_clusters_annotated$symbol))

exemplar_rows <- sig_clusters_annotated %>%
  mutate(hit_parcels = map_chr(hit_parcels, paste, collapse = ";")) %>%
  dplyr::filter(symbol %in% exemplar_present) %>%
  group_by(symbol) %>%
  slice_min(fdr, n = 1, with_ties = FALSE) %>%
  ungroup()

# Coloured by the phenotype's own measure colour (phenotype_colours)
# rather than outlier-group membership -- the outlier framing has
# been dropped from the manuscript, so this panel no longer needs
# to distinguish paralimbic/insula-SM outlier parcels at all.
make_gene_brain <- function(hit_p, sf_df, title_str, measure_col) {
  sf_df %>%
    mutate(fill_col = ifelse(parcel %in% hit_p, measure_col, "#EEEEEE")) %>%
    ggplot() +
    geom_sf(aes(geometry = geometry, fill = fill_col), colour = "white", linewidth = 0.1) +
    coord_sf(expand = FALSE) +
    scale_fill_identity() +
    theme_void(base_size = BASE_SIZE) +
    labs(title = title_str) +
    theme(plot.title = element_text(size = LABEL_SIZE, hjust = 0.5, face = "bold",
                                    margin = margin(b = 1)),
          plot.margin = margin(1, 1, 1, 1))
}

# Reverting to the same structure as the original (4-gene) version of
# this script, which never had a sizing problem: ONE hemisphere (left
# -- valid since the parcellation is bilaterally averaged, so the left
# side fully represents the result), shown as a lateral+medial PAIR
# per gene, at a fixed practical size. The bilateral "both hemispheres
# in one image" approach produced a very wide bounding box (two
# silhouettes with a gap between them); deriving the exact aspect
# ratio from that pushed the cell width past ggsave's 50-inch limit.
# Fixed sizing avoids that failure mode entirely, at the cost of only
# showing the left hemisphere (as before).
panelA_plots <- pmap(
  list(exemplar_rows$symbol, exemplar_rows$phenotype_type, exemplar_rows$hit_parcels),
  function(sym, phen, hp_str) {
    hit_p <- str_split(hp_str, ";")[[1]]
    mcol  <- phenotype_colours[[phen]]
    if (is.null(mcol) || is.na(mcol)) mcol <- "#B2182B"  # fallback if phen not in the lookup
    ttl   <- sprintf("%s (%s)", sym, phen)
    p_lat <- make_gene_brain(hit_p, lh_sf_lat, ttl, mcol)
    p_med <- make_gene_brain(hit_p, lh_sf_med, "", mcol)
    plot_grid(p_lat, p_med, nrow = 1, align = "h")
  })

# 6 genes, 2 rows x 3 cols (each cell = one gene's lat+med pair).
# Fixed size, same as the original 2x2/6x6 pattern scaled to 3 cols.
panelA_w <- 9
panelA_h <- 6
panel_A <- plot_grid(plotlist = panelA_plots, nrow = 2, ncol = 3, align = "hv")
ggsave(file.path(figures_dir, "PanelA_exemplars.pdf"), panel_A, width = panelA_w, height = panelA_h)

# ---- Panel B: within/between family Jaccard, permutation test --
# Same rationale as Figure 4's compute_family_similarity(): a
# Wilcoxon rank-sum test assumes independent samples, but pairwise
# Jaccard values are NOT independent here -- each gene and each
# phenotype_type appears in multiple pairs -- so its p-value would
# be anti-conservative. Instead, permute the exchangeable unit
# (the phenotype_type -> family label, n_perm shuffles) and
# recompute the within-minus-between mean-Jaccard difference under
# each shuffle, exactly mirroring Figure 4's Supp Fig 2 approach.
set.seed(42)
N_PERM_JAC <- 10000

# Rebuild ALL pairwise gene-level parcel-overlap (Jaccard) comparisons
# directly from sig_lookup (symbol, phenotype_type, family, hit_parcels)
# -- one row per pair of phenotype_type entries sharing the same gene --
# rather than reusing the pre-formatted same/other-family summary
# strings, so every pair's phenotype_type labels are available for
# permutation.
gene_pair_jaccard <- function(lookup_df) {
  lookup_df %>%
    group_by(symbol) %>%
    dplyr::filter(dplyr::n() >= 2) %>%
    group_modify(~ {
      idx <- combn(nrow(.x), 2)
      tibble(
        phen_a  = .x$phenotype_type[idx[1, ]],
        phen_b  = .x$phenotype_type[idx[2, ]],
        jaccard = map2_dbl(.x$hit_parcels[idx[1, ]], .x$hit_parcels[idx[2, ]], jaccard)
      )
    }) %>%
    ungroup()
}

pair_df_jac <- gene_pair_jaccard(sig_lookup) %>%
  mutate(fam_a = family_map[phen_a], fam_b = family_map[phen_b],
         within = fam_a == fam_b) %>%
  dplyr::filter(!is.na(jaccard))

within_vals  <- pair_df_jac$jaccard[pair_df_jac$within]
between_vals <- pair_df_jac$jaccard[!pair_df_jac$within]
obs_diff_jac <- mean(within_vals, na.rm = TRUE) - mean(between_vals, na.rm = TRUE)

# Null: shuffle the phenotype_type -> family mapping itself (the
# independent unit), then relabel every observed pair under that
# shuffle and recompute the within/between difference -- NOT
# shuffling the derived pairwise Jaccard values directly.
phenos_present_jac <- names(family_map)
null_diffs_jac <- replicate(N_PERM_JAC, {
  perm_fam <- setNames(sample(family_map), phenos_present_jac)
  w <- perm_fam[pair_df_jac$phen_a] == perm_fam[pair_df_jac$phen_b]
  mean(pair_df_jac$jaccard[w], na.rm = TRUE) - mean(pair_df_jac$jaccard[!w], na.rm = TRUE)
})
perm_p_jac <- mean(null_diffs_jac >= obs_diff_jac)

write.csv(pair_df_jac, file.path(tables_dir, "jaccard_within_between_pairs.csv"), row.names = FALSE)
cat(sprintf("  Jaccard within/between: obs_diff=%.4f, perm_p=%s (n_perm=%d)\n",
            obs_diff_jac, formatC(perm_p_jac, format = "e", digits = 2), N_PERM_JAC))

jaccard_df <- bind_rows(
  tibble(group = "Within",  jaccard = within_vals),
  tibble(group = "Between", jaccard = between_vals)
) %>% mutate(group = factor(group, levels = c("Within", "Between")))

p_B <- ggplot(jaccard_df, aes(x = group, y = jaccard, fill = group)) +
  geom_violin(alpha = 0.7, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("Within" = "#B2182B", "Between" = "#AAAAAA"), guide = "none") +
  labs(x = NULL, y = "Jaccard", tag = "B",
       subtitle = sprintf("Δ = %.3f, perm p = %s",
                          obs_diff_jac, formatC(perm_p_jac, format = "e", digits = 2))) +
  theme_fig6(BASE_SIZE - 2) +
  theme(axis.text = element_text(size = LABEL_SIZE - 2),
        plot.tag = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = LABEL_SIZE - 2))

# ---- Panel C: overall Qian and Bhaduri, independently -----------
panelC_df <- enrich_df %>%
  dplyr::filter(level == "overall") %>%
  mutate(log2OR = log2(pmax(OR, 0.01)),
         sig = ifelse(fdr < 0.05, "FDR < 0.05", "ns"),
         label = source)

p_C <- ggplot(panelC_df, aes(x = log2OR, y = label, colour = sig, size = n_overlap)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_point() +
  scale_colour_manual(values = c("FDR < 0.05" = "#B2182B", "ns" = "#AAAAAA"), guide = "none") +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(x = "log2 OR", y = NULL, tag = "C") +
  theme_fig6(BASE_SIZE - 2) +
  theme(axis.text = element_text(size = LABEL_SIZE - 2),
        plot.tag = element_text(face = "bold", size = 11))

# ---- Panel E: pathway enrichment, >50% threshold ----------------
pathway_df <- combined %>% dplyr::filter(threshold == 0.50)

panelE_df <- pathway_df %>%
  mutate(label = gsub("^(GOBP_|KEGG_|REACTOME_)", "", pathway),
         label = gsub("_", " ", label), label = tools::toTitleCase(tolower(label)),
         label = str_wrap(label, width = 22),   # spill to 2 lines -- E is now a narrower column
         log2OR = log2(pmax(OR, 0.01)),
         sig = case_when(fdr < 0.05 ~ "FDR < 0.05", p_value < 0.05 ~ "p < 0.05", TRUE ~ "ns"),
         group = factor(group, levels = c("frontal", "posterior", "intermediate")))

sig_pathways_e <- panelE_df %>% dplyr::filter(p_value < 0.05) %>% pull(pathway) %>% unique()
panelE_df <- panelE_df %>% dplyr::filter(pathway %in% sig_pathways_e)

p_E <- ggplot(panelE_df, aes(x = group, y = reorder(label, log2OR), size = OR, colour = sig)) +
  geom_point() +
  scale_colour_manual(values = c("FDR < 0.05" = "#B2182B", "p < 0.05" = "#F46D43", "ns" = "#AAAAAA"), name = NULL) +
  scale_size_continuous(range = c(2, 7), name = "Odds Ratio") +
  # drop=FALSE keeps "posterior" visible on the x-axis even where it
  # has no significant (or no tested) hits, rather than silently
  # disappearing because ggplot drops empty factor levels by default.
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, tag = "D", title = "Pathway enrichment\n(>50% territory)") +
  theme_fig6(BASE_SIZE - 1) +
  theme(axis.text = element_text(size = LABEL_SIZE - 1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        plot.title = element_text(size = TITLE_SIZE - 1),
        legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 11))

ggsave(file.path(figures_dir, "PanelB_jaccard.pdf"), p_B, width = 3, height = 3)
ggsave(file.path(figures_dir, "PanelC_DEG_overall.pdf"), p_C, width = 4, height = 2.5)
ggsave(file.path(figures_dir, "PanelE_pathway.pdf"), p_E, width = 4,
       height = max(3, 0.5 * n_distinct(panelE_df$label) + 1.5), limitsize = FALSE)

# ---- Assembly ----------------------------------------------------
# Top row:    Panel A only -- 6 exemplar gene brain plots (one
#             bilateral view each), 2 rows x 3 cols, sized to the
#             brain silhouette's true aspect ratio (see above).
# Bottom row: B (Jaccard), C (Qian/Bhaduri overall enrichment), and
#             D/pathway (GO:BP enrichment, >50% AP threshold) side
#             by side. Panel D-eRegulon has been dropped.
#
# Use panelA_w/panelA_h (computed above from the actual brain
# silhouette aspect ratio) to fix the top row's height in real
# inches, then give the bottom row a sensible fixed height rather
# than a relative fraction -- rel_heights on a plot_grid() forces
# both rows into whatever total figure height you pick, so if that
# total doesn't already match panel_A's true proportions, you're
# back to guessing and re-introducing letterboxing one level up.
bottom_row_h <- 4
fig6_w <- max(panelA_w, 14)
fig6_h <- panelA_h + bottom_row_h

bottom_row <- plot_grid(p_B, p_C, p_E, nrow = 1, rel_widths = c(1, 1, 1.3))
fig6 <- plot_grid(panel_A, bottom_row, ncol = 1,
                  rel_heights = c(panelA_h, bottom_row_h))

ggsave(file.path(out_dir, "Figure6_arealisation.pdf"), fig6, width = 12, height = 8)
ggsave(file.path(out_dir, "Figure6_arealisation.png"), fig6, width = 12, height = 8, dpi = 300)

cat("\n=== ALL DONE ===\n")
cat("Tables:  ", tables_dir, "\n")
cat("Panels:  ", figures_dir, "\n")
cat("Figure:  ", file.path(out_dir, "Figure6_arealisation.pdf/.png"), "\n")