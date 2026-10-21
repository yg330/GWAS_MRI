# ===========================================================================
# Final_Figure5_TFenrichement_script.R
# Manuscript item: Supplementary Figure 7 (TF / eRegulon analyses)
#
# What this does:
#   FLAMES gene and transcription-factor expression across developmental age
#   and pseudotime, and eRegulon (AUCell) activity by cell type for each
#   family.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Figure 5 (v3): Gene regulatory mechanisms underlying brain
#                structure
#
# PANELS
# ------
# A: Paired OR forest plot (unchanged from the original Fig 6
#    script).
#
# B: Three-row stack, chronological age x-axis:
#      1. Relative expression of each family's FULL FLAMES gene
#         set
#      2. Relative expression of each family's FLAMES TFs
#      3. Broad cell-type proportion (ONE shared banner, not
#         duplicated per row -- same underlying cells/bins)
#
# C: Identical three-row stack, but x-axis = average pseudotime
#    instead of chronological age. Same y-scale as B throughout
#    (computed jointly) so all four trajectory curves (age-gene,
#    age-TF, pseudotime-gene, pseudotime-TF) are directly
#    comparable.
#
# D: Family x broad-cell-type heatmap of mean AUCell z-score,
#    split Activating | Repressing, across ALL significant
#    FLAMES-overlapping eRegulons.
#
# RELATIVE EXPRESSION: at each point in time (age or pseudotime),
# z-score every gene against the full transcriptome measured at
# that same point -- genes binned, pseudobulked per bin, z-scored
# within each bin against the whole gene universe, then
# GAM-smoothed across bins.

source("config.R")

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(mgcv)
  library(biomaRt)
  library(patchwork)
  library(scales)
  library(grid)
})

select   <- dplyr::select;   filter    <- dplyr::filter
mutate   <- dplyr::mutate;   arrange   <- dplyr::arrange
rename   <- dplyr::rename;   summarise <- dplyr::summarise
group_by <- dplyr::group_by; ungroup   <- dplyr::ungroup

ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
set.seed(42)

# ---- Colour maps / orderings ----
family_colours <- c(
  "Brain size"                  = "#00C28A",
  "Cortical thickness"          = "#0072B2",
  "Curvature"                   = "#E69F00",
  "Microstructural coherence"   = "#CC79A7",
  "Microstructural diffusivity" = "#56B4E9",
  "Orientation dispersion"      = "#D55E00"
)
family_map <- c(
  SA = "Brain size", FI = "Brain size", GMV = "Brain size", IC = "Brain size",
  CT = "Cortical thickness",
  GC = "Curvature", MC = "Curvature",
  FA = "Microstructural coherence", ICVF = "Microstructural coherence",
  MD = "Microstructural diffusivity", ISOVF = "Microstructural diffusivity",
  OD = "Orientation dispersion"
)
family_order <- c(
  "Brain size", "Cortical thickness", "Curvature",
  "Microstructural coherence", "Microstructural diffusivity",
  "Orientation dispersion"
)
phenotype_order <- c("SA","FI","GMV","IC","CT",
                     "GC","MC","FA","ICVF","MD","ISOVF","OD")

broad_order <- c("Progenitors","Excitatory neurons","Inhibitory neurons",
                 "Cajal-Retzius","Astrocytes","OPC","Oligodendrocyte",
                 "Microglia")
# Progenitors = purple, Excitatory neurons = brown (matches the
# lineage colouring used in the fine-scale pseudotime strip plot
# below), everything else = shades of grey through to black.
grey_ramp <- colorRampPalette(c("#8C8C8C", "#000000"))(6)
broad_colours <- c(
  "Progenitors"         = "#6A51A3",
  "Excitatory neurons"  = "#8C510A",
  "Inhibitory neurons"  = grey_ramp[1],
  "Cajal-Retzius"       = grey_ramp[2],
  "Astrocytes"          = grey_ramp[3],
  "OPC"                 = grey_ramp[4],
  "Oligodendrocyte"     = grey_ramp[5],
  "Microglia"           = grey_ramp[6]
)

collapse_subclass <- function(ct) {
  dplyr::case_when(
    ct %in% c("Radial glia", "IPC-EN", "IPC-Glia") ~ "Progenitors",
    ct == "Glutamatergic neuron"                   ~ "Excitatory neurons",
    ct == "GABAergic neuron"                       ~ "Inhibitory neurons",
    ct == "Cajal-Retzius cell"                     ~ "Cajal-Retzius",
    ct == "Astrocyte"                              ~ "Astrocytes",
    ct == "OPC"                                    ~ "OPC",
    ct == "Oligodendrocyte"                        ~ "Oligodendrocyte",
    ct == "Microglia"                              ~ "Microglia",
    TRUE                                           ~ "Other"
  )
}

# ---- Global typography (matches Figure 4 / manuscript convention) ----
BASE_SIZE  <- 10
TITLE_SIZE <- 10
LABEL_SIZE <- 9

# ---- Extract a clean legend guide from a dummy plot ----
# theme_void() hides axis furniture but NOT the underlying geoms
# themselves -- embedding the raw plot object into a layout lets
# those geoms render at whatever size the layout slot gives them,
# which is what produced the giant colour-ribbon artifact. Pulling
# just the "guide-box" grob out avoids that entirely: only the
# actual legend keys + text get placed, never the plot panel.
extract_legend <- function(p) {
  built <- ggplot_gtable(ggplot_build(p))
  idx   <- which(sapply(built$grobs, function(x) x$name) == "guide-box")
  if (length(idx) == 0) return(nullGrob())
  built$grobs[[idx]]
}

# ---- Paths ----
base_dir  <- proj_dir
data_dir  <- file.path(base_dir, "Data")
cache_dir <- file.path(base_dir, "Outputs", "cache")
fig6_dir  <- file.path(base_dir, "Outputs", "Figure_6")
out_dir   <- file.path(base_dir, "Outputs", "Figure_5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

flames_path   <- file.path(proj_dir, "Data", "FLAMES_gene_prioritisation_results.txt")
eregulon_path <- file.path(data_dir, "Wang_eregulonlist.xlsx")
aucell_path   <- file.path(data_dir, "Wang_eregulon_aucell.xlsx")
wang_path     <- file.path(cache_dir, "wang_filtered.rds")
meta_path     <- file.path(data_dir, "Wang_celltypemetadata.xlsx")
pc_path       <- file.path(cache_dir, "pseudobulk_protein_coding.rds")
gl_cache_path <- file.path(cache_dir, "gene_lengths_hg19.rds")

enrich_cache_path <- file.path(fig6_dir, "Fig6_enrichment_phenotype_full.csv")

# ============================================================
# SECTION 1: Gene universe, FLAMES genes, gene lengths
# ============================================================
ts("Loading gene universe...")
pc_data <- readRDS(pc_path)
universe_genes <- rownames(pc_data)
universe_genes <- universe_genes[!is.na(universe_genes) & universe_genes != ""]
rm(pc_data); gc()

flames_raw <- read.table(flames_path, header = TRUE, stringsAsFactors = FALSE)
flames_genes <- flames_raw %>%
  select(symbol, phenotype) %>%
  distinct() %>%
  mutate(
    pheno_short = str_extract(phenotype, paste(names(family_map), collapse = "|")),
    family = family_map[pheno_short]
  ) %>%
  filter(!is.na(family))

phenotypes <- intersect(phenotype_order, unique(flames_genes$pheno_short))
ph_gene_sets <- lapply(phenotypes, function(ph) {
  flames_genes %>% filter(pheno_short == ph) %>% pull(symbol) %>% unique()
})
names(ph_gene_sets) <- phenotypes

if (file.exists(gl_cache_path)) {
  gene_lengths <- readRDS(gl_cache_path)
} else {
  ts("Fetching gene lengths...")
  mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
  gene_lengths <- getBM(
    attributes = c("hgnc_symbol","start_position","end_position"), mart = mart
  ) %>%
    filter(hgnc_symbol != "") %>%
    mutate(gene_length = end_position - start_position) %>%
    group_by(hgnc_symbol) %>%
    summarise(gene_length = max(gene_length), .groups = "drop") %>%
    rename(symbol = hgnc_symbol)
  saveRDS(gene_lengths, gl_cache_path)
}
universe_df <- tibble(symbol = universe_genes) %>%
  left_join(gene_lengths, by = "symbol") %>%
  mutate(
    gene_length = replace_na(gene_length, median(gene_lengths$gene_length, na.rm = TRUE)),
    log_gene_length = log(gene_length)
  )

genes_by_family <- flames_genes %>%
  distinct(symbol, family) %>%
  group_by(family) %>%
  summarise(genes = list(sort(unique(symbol))), .groups = "drop")

# ============================================================
# SECTION 2: eRegulon list + FLAMES TF identity
# ============================================================
ts("Loading eRegulon list...")
ereg_raw <- read_excel(eregulon_path)
colnames(ereg_raw) <- str_trim(colnames(ereg_raw))
ereg_list <- ereg_raw %>%
  rename(eregulon = 1, tf_full = 2, target_str = 3) %>%
  mutate(
    tf = str_remove(tf_full, "_extended.*|_[+-]_.*|_\\(.*"),
    target_genes = str_split(target_str, ",\\s*"),
    polarity_label = case_when(
      str_detect(eregulon, "_\\+_\\+") ~ "Activating",
      str_detect(eregulon, "_-_\\+")   ~ "Repressing",
      TRUE                              ~ "Activating"
    )
  )
all_tf_symbols   <- unique(ereg_list$tf)
all_target_genes <- unique(unlist(ereg_list$target_genes))

if (!file.exists(enrich_cache_path)) {
  stop("Missing cached enrichment results at:\n  ", enrich_cache_path,
       "\nRun the original Fig 6 script's Section 3b once to generate ",
       "Fig6_enrichment_phenotype_full.csv before running this script.")
}
enrich_results <- read.csv(enrich_cache_path, stringsAsFactors = FALSE) %>%
  mutate(phenotype = factor(phenotype, levels = phenotype_order),
         family    = factor(family,    levels = family_order))

flames_tfs <- enrich_results %>%
  # Previously also required filter(sig != "ns") here -- i.e. the TF's
  # eRegulon had to be significant in the separate, earlier per-eRegulon
  # enrichment analysis (Fig6_enrichment_phenotype_full.csv), on top of
  # the TF being a FLAMES gene. Removed per request: restricting to
  # eRegulon-significant TFs was dropping FLAMES TFs whose eRegulon
  # target-gene set didn't independently clear that separate test.
  # This affects Panels B and C, which track TF genes' own expression --
  # for that question, requiring the TF itself to be a FLAMES-prioritised
  # gene is the right criterion. Panel D (below) now uses a DIFFERENT,
  # independent criterion -- see ereg_family_polarity -- since it's
  # asking about eRegulon activity, a property of the regulatory unit,
# not of the TF gene's own locus.
rowwise() %>%
  mutate(tf_is_flames = tf %in% ph_gene_sets[[as.character(phenotype)]]) %>%
  ungroup() %>%
  filter(tf_is_flames)
ts("  FLAMES TFs (FLAMES-gene membership only):", length(unique(flames_tfs$tf)))

tfs_by_family <- flames_tfs %>%
  distinct(tf, family) %>%
  group_by(family) %>%
  summarise(genes = list(sort(unique(tf))), .groups = "drop")

# ============================================================
# SECTION 3: PANEL A -- paired OR forest plot (unchanged)
# ============================================================
ts("Building Panel A (paired OR forest)...")

run_paired_enrichment <- function(ph_genes, ph_label, universe_df,
                                  tf_symbols, target_genes) {
  df <- universe_df %>%
    mutate(is_flames = as.integer(symbol %in% ph_genes),
           is_tf     = as.integer(symbol %in% tf_symbols),
           is_target = as.integer(symbol %in% target_genes))
  run_glm <- function(predictor, label) {
    n_overlap <- sum(df[[predictor]] == 1 & df$is_flames == 1)
    if (n_overlap < 2) return(NULL)
    form <- as.formula(paste("is_flames ~", predictor, "+ log_gene_length"))
    fit <- tryCatch(glm(form, data = df, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(cf) || !predictor %in% rownames(cf)) return(NULL)
    tibble(phenotype = ph_label, family = family_map[ph_label], type = label,
           n_overlap = n_overlap, OR = exp(cf[predictor, "Estimate"]),
           CI_lo = exp(cf[predictor, "Estimate"] - 1.96 * cf[predictor, "Std. Error"]),
           CI_hi = exp(cf[predictor, "Estimate"] + 1.96 * cf[predictor, "Std. Error"]),
           p_value = cf[predictor, "Pr(>|z|)"])
  }
  bind_rows(run_glm("is_tf", "TF"), run_glm("is_target", "Target genes"))
}

paired_results <- map_dfr(phenotypes, function(ph) {
  run_paired_enrichment(ph_gene_sets[[ph]], ph, universe_df, all_tf_symbols, all_target_genes)
}) %>%
  group_by(type) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(sig = case_when(p_fdr < 0.001 ~ "***", p_fdr < 0.01 ~ "**",
                         p_fdr < 0.05 ~ "*", TRUE ~ "ns"),
         phenotype = factor(phenotype, levels = rev(phenotype_order)),
         family    = factor(family, levels = family_order),
         type      = factor(type, levels = c("TF", "Target genes")))

write.csv(paired_results, file.path(out_dir, "Fig5A_paired_enrichment.csv"), row.names = FALSE)

dodge_a <- position_dodge(width = 0.6)
sig_results <- paired_results %>% filter(sig != "ns")

# Family boundary lines: compute directly from the phenotypes
# ACTUALLY PRESENT in the significance-filtered data, in the same
# REVERSED order the y-axis factor uses (levels = rev(phenotype_order)
# above). The previous version computed these from the un-reversed
# phenotype_order while the plot itself uses the reversed order --
# that mismatch is what put the lines at the wrong y-positions.
phs_present   <- intersect(rev(phenotype_order), as.character(unique(sig_results$phenotype)))
fams_present  <- family_map[phs_present]
boundaries    <- which(fams_present != dplyr::lag(fams_present, default = ""))
ph_fam_df     <- data.frame(yintercept = boundaries - 0.5) %>% filter(yintercept > 0.5)

p_A <- ggplot(sig_results,
              aes(x = OR, y = phenotype, linetype = type, shape = type, group = type)) +
  geom_hline(data = ph_fam_df, aes(yintercept = yintercept),
             linetype = "dotted", colour = "grey60", linewidth = 0.4, inherit.aes = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.25, linewidth = 0.5,
                 colour = "black", position = dodge_a) +
  geom_point(size = 3, alpha = 0.9, colour = "black", position = dodge_a) +
  scale_shape_manual(values = c("TF" = 16, "Target genes" = 17), name = "Enrichment unit") +
  scale_linetype_manual(values = c("TF" = "solid", "Target genes" = "22"), name = "Enrichment unit") +
  scale_x_log10() +
  theme_bw(base_size = BASE_SIZE) +
  theme(axis.text.y = element_text(size = LABEL_SIZE), axis.text.x = element_text(size = LABEL_SIZE),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = LABEL_SIZE),
        legend.title = element_text(size = LABEL_SIZE)) +
  labs(x = "Odds ratio (log scale)", y = NULL, tag = "a")

ts("  Panel A done")

# ============================================================
# SECTION 4: Load Wang cells, restrict gene universe, join
#            BOTH age and pseudotime
# ============================================================
ts("Loading Wang SCE...")

wang       <- readRDS(wang_path)
keep_cells <- !colData(wang)$Class %in% c("Vascular", "Unknown")
wang       <- wang[, keep_cells]

genes_present <- intersect(universe_genes, rownames(wang))
wang_uni      <- wang[genes_present, ]
rm(wang); gc()
ts("  Universe genes present in Wang:", length(genes_present))

n_samp    <- min(15000, ncol(wang_uni))
samp_idx  <- sample(ncol(wang_uni), n_samp)
expr_samp <- assay(wang_uni[, samp_idx], "X")
cell_meta <- as.data.frame(colData(wang_uni[, samp_idx]))
cell_meta$cell_id <- colnames(wang_uni)[samp_idx]
rm(wang_uni); gc()

if (!"Subclass" %in% colnames(cell_meta)) {
  stop("colData(wang) has no 'Subclass' column. Check the object structure.")
}

wang_meta   <- read_excel(meta_path)
cell_id_col <- colnames(wang_meta)[1]
age_col     <- grep("postconceptional|age_in_days", colnames(wang_meta),
                    ignore.case = TRUE, value = TRUE)[1]
# NOTE: wang_meta has multiple per-lineage pseudotime columns
# (EN.L2_3.IT_pseudotime, RG_pseudotime, etc. -- each mostly NA
# outside its own lineage) alongside the one we actually want,
# Average_pseudotime. Match that column name exactly rather than
# a generic "pseudotime" pattern, which was grabbing whichever
# per-lineage column happened to appear first.
pt_col <- grep("^Average_pseudotime$", colnames(wang_meta),
               ignore.case = TRUE, value = TRUE)[1]
if (is.na(pt_col)) stop("Could not find an 'Average_pseudotime' column in wang_meta.")

cell_meta <- cell_meta %>%
  left_join(wang_meta %>% select(all_of(c(cell_id_col, age_col, pt_col))),
            by = setNames(cell_id_col, "cell_id"))
age_col_meta <- grep("postconceptional|age_in_days", colnames(cell_meta),
                     ignore.case = TRUE, value = TRUE)[1]
pt_col_meta  <- grep("^Average_pseudotime$", colnames(cell_meta),
                     ignore.case = TRUE, value = TRUE)[1]
cell_meta$pca_days   <- cell_meta[[age_col_meta]]
# Average_pseudotime is stored as character with literal "NA"
# strings mixed with real NA and numeric-looking values -- convert
# explicitly rather than relying on implicit coercion.
pt_raw <- cell_meta[[pt_col_meta]]
if (is.character(pt_raw)) pt_raw[pt_raw == "NA"] <- NA
cell_meta$pseudotime <- as.numeric(pt_raw)

ts("  Cells with valid age:", sum(!is.na(cell_meta$pca_days) & cell_meta$pca_days > 0),
   "/", nrow(cell_meta))
ts("  Cells with valid pseudotime:", sum(!is.na(cell_meta$pseudotime)), "/", nrow(cell_meta))

# ============================================================
# SECTION 5: Generic binning -> pseudobulk -> z-score, and
#            broad cell-type proportions, reused for both age
#            and pseudotime
# ============================================================
ts("Binning: expression z-scores + cell-type proportions...")

# x_vals: the continuous variable to bin on (already log-transformed
# for age upstream if desired -- this function itself is agnostic)
bin_and_zscore <- function(x_vals, expr_mat, n_bins = 24) {
  valid_idx  <- which(!is.na(x_vals))
  x_valid    <- x_vals[valid_idx]
  bin_breaks <- seq(min(x_valid), max(x_valid), length.out = n_bins + 1)
  bin_id     <- cut(x_valid, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE)
  bin_mid    <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2
  
  expr_valid  <- expr_mat[, valid_idx, drop = FALSE]
  gene_by_bin <- sapply(seq_len(n_bins), function(b) {
    cols <- which(bin_id == b)
    if (length(cols) < 10) return(rep(NA_real_, nrow(expr_valid)))
    Matrix::rowMeans(expr_valid[, cols, drop = FALSE], na.rm = TRUE)
  })
  rownames(gene_by_bin) <- rownames(expr_valid)
  log_gene_by_bin <- log10(gene_by_bin + 0.01)
  z_gene_by_bin <- apply(log_gene_by_bin, 2, function(col) {
    (col - mean(col, na.rm = TRUE)) / sd(col, na.rm = TRUE)
  })
  rownames(z_gene_by_bin) <- rownames(log_gene_by_bin)
  
  valid_bins <- colSums(!is.na(z_gene_by_bin)) > (nrow(z_gene_by_bin) * 0.5)
  list(z_mat = z_gene_by_bin[, valid_bins, drop = FALSE],
       bin_mid = bin_mid[valid_bins],
       valid_idx = valid_idx, bin_id = bin_id, valid_bins = valid_bins)
}

build_banner <- function(bin_result, subclass_vals, x_label_name) {
  broad_valid <- collapse_subclass(subclass_vals[bin_result$valid_idx])
  tibble(bin = bin_result$bin_id, broad = broad_valid) %>%
    filter(!is.na(bin), broad != "Other") %>%
    count(bin, broad, name = "n_cells") %>%
    group_by(bin) %>%
    mutate(proportion = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    mutate(x = bin_result$bin_mid[match(bin, which(bin_result$valid_bins))],
           broad = factor(broad, levels = broad_order)) %>%
    filter(!is.na(x)) %>%
    rename(!!x_label_name := x)
}

# -- Age (log-transformed) --
log_age_full <- ifelse(!is.na(cell_meta$pca_days) & cell_meta$pca_days > 0,
                       log(cell_meta$pca_days), NA_real_)
bin_age <- bin_and_zscore(log_age_full, expr_samp, n_bins = 24)
banner_age_df <- build_banner(bin_age, cell_meta$Subclass, "log_age")

other_check_age <- table(cell_meta$Subclass[bin_age$valid_idx][
  collapse_subclass(cell_meta$Subclass[bin_age$valid_idx]) == "Other"])
if (length(other_check_age) > 0) {
  cat("\n'Other' in Subclass, age binning (should be empty/near-empty):\n")
  print(other_check_age); cat("\n")
}

# -- Pseudotime (linear, no log transform) --
bin_pt <- bin_and_zscore(cell_meta$pseudotime, expr_samp, n_bins = 24)
banner_pt_df <- build_banner(bin_pt, cell_meta$Subclass, "pseudotime")

write.csv(banner_age_df, file.path(out_dir, "Fig5_celltype_banner_age.csv"), row.names = FALSE)
write.csv(banner_pt_df, file.path(out_dir, "Fig5_celltype_banner_pseudotime.csv"), row.names = FALSE)

# ============================================================
# SECTION 6: GAM-smooth family trajectories -- shared helper,
#            applied to both age and pseudotime bins, both gene
#            sets (FLAMES genes / FLAMES TFs)
# ============================================================
smooth_family_trajectory <- function(gene_set_df, z_mat, bin_x, x_label_name) {
  purrr::map_dfr(seq_len(nrow(gene_set_df)), function(i) {
    fam   <- gene_set_df$family[i]
    genes <- intersect(gene_set_df$genes[[i]], rownames(z_mat))
    if (length(genes) < 3) return(NULL)
    sub      <- z_mat[genes, , drop = FALSE]
    fam_mean <- colMeans(sub, na.rm = TRUE)
    d <- data.frame(x = bin_x, y = fam_mean)
    d <- d[!is.na(d$y), ]
    if (nrow(d) < 8) return(NULL)
    fit <- tryCatch(gam(y ~ s(x, k = 7, bs = "tp"), data = d), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    x_seq <- seq(min(d$x), max(d$x), length.out = 200)
    pred  <- predict(fit, newdata = data.frame(x = x_seq), se.fit = TRUE)
    out <- tibble(family = fam, x = x_seq, z = as.numeric(pred$fit),
                  ci_lo = as.numeric(pred$fit) - 1.96 * as.numeric(pred$se.fit),
                  ci_hi = as.numeric(pred$fit) + 1.96 * as.numeric(pred$se.fit),
                  n_genes = length(genes))
    rename(out, !!x_label_name := x)
  }) %>%
    mutate(family = factor(family, levels = family_order))
}

ts("Smoothing trajectories (age)...")
panelB_age_df <- smooth_family_trajectory(genes_by_family, bin_age$z_mat, bin_age$bin_mid, "log_age")
panelC_age_df <- smooth_family_trajectory(tfs_by_family,   bin_age$z_mat, bin_age$bin_mid, "log_age")

ts("Smoothing trajectories (pseudotime)...")
panelB_pt_df <- smooth_family_trajectory(genes_by_family, bin_pt$z_mat, bin_pt$bin_mid, "pseudotime")
panelC_pt_df <- smooth_family_trajectory(tfs_by_family,   bin_pt$z_mat, bin_pt$bin_mid, "pseudotime")

write.csv(panelB_age_df, file.path(out_dir, "Fig5B_flames_gene_age.csv"), row.names = FALSE)
write.csv(panelC_age_df, file.path(out_dir, "Fig5B_flames_tf_age.csv"), row.names = FALSE)
write.csv(panelB_pt_df, file.path(out_dir, "Fig5C_flames_gene_pseudotime.csv"), row.names = FALSE)
write.csv(panelC_pt_df, file.path(out_dir, "Fig5C_flames_tf_pseudotime.csv"), row.names = FALSE)

# Shared y-axis across ALL FOUR trajectory curves -- age and
# pseudotime versions stay directly comparable to each other
shared_ylim <- range(c(panelB_age_df$ci_lo, panelB_age_df$ci_hi,
                       panelC_age_df$ci_lo, panelC_age_df$ci_hi,
                       panelB_pt_df$ci_lo,  panelB_pt_df$ci_hi,
                       panelC_pt_df$ci_lo,  panelC_pt_df$ci_hi), na.rm = TRUE)

# ============================================================
# SECTION 7: Three-row stack (gene expr / TF expr / cell-type %)
#            -- shared helper, used for both age and pseudotime
# ============================================================
make_trio <- function(gene_df, tf_df, bottom_panel, x_col, ylim,
                      x_breaks, x_labels, tag, bottom_height = 1.3) {
  base_traj <- function(df, y_lab, show_tag = NULL) {
    ggplot(df, aes(x = .data[[x_col]], y = z, colour = family,
                   fill = family, group = family)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 1.1) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels, name = NULL) +
      scale_colour_manual(values = family_colours, guide = "none") +
      scale_fill_manual(values = family_colours, guide = "none") +
      scale_y_continuous(limits = ylim) +
      labs(y = y_lab, tag = show_tag) +
      theme_classic(base_size = BASE_SIZE) +
      theme(axis.text = element_text(size = LABEL_SIZE),
            axis.title = element_text(size = TITLE_SIZE),
            legend.position = "none", axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  p_gene <- base_traj(gene_df, "FLAMES genes:\nz-score", tag)
  p_tf   <- base_traj(tf_df,   "FLAMES TFs:\nz-score")
  
  (p_gene / p_tf / bottom_panel) + plot_layout(heights = c(2, 2, bottom_height))
}

# -- Panel B bottom: cell-type proportion banner (age; covers the
#    whole population, so a stacked-proportion area chart is valid) --
make_proportion_banner <- function(banner_df, x_col, x_breaks, x_labels, x_title) {
  ggplot(banner_df, aes(x = .data[[x_col]], y = proportion, fill = broad, group = broad)) +
    geom_area(position = "stack", alpha = 0.9, colour = "white", linewidth = 0.1) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, name = x_title) +
    scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = broad_colours, guide = "none") +
    labs(y = "Cell\ntype %") +
    theme_classic(base_size = BASE_SIZE - 2) +
    theme(axis.text = element_text(size = LABEL_SIZE - 2.5),
          axis.title = element_text(size = TITLE_SIZE - 2),
          legend.position = "none")
}

# -- Panel C bottom: per-cell strip plot (pseudotime; only defined
#    for the RG -> EN lineage, so a proportion chart would be
#    misleading -- shows individual cells per fine subtype instead) --
# Confirmed via colData(wang)$pseudotime coverage: IPC-Glia has ZERO
# non-NA pseudotime values despite being an IPC in the Type column,
# so it's excluded here.
fine_pt_types <- c("RG-vRG","RG-oRG","RG-tRG","IPC-EN",
                   "EN-Newborn","EN-IT-Immature","EN-Non-IT-Immature",
                   "EN-L2_3-IT","EN-L4-IT","EN-L5-IT","EN-L5_6-NP","EN-L5-ET",
                   "EN-L6-CT","EN-L6-IT","EN-L6b")

make_pseudotime_strip <- function(cell_meta, fine_types, x_breaks, x_labels, x_title) {
  strip_df <- cell_meta %>%
    filter(!is.na(pseudotime), Type %in% fine_types) %>%
    select(Type, pseudotime)
  type_order <- strip_df %>%
    group_by(Type) %>%
    summarise(m = mean(pseudotime), .groups = "drop") %>%
    arrange(m) %>%
    pull(Type)
  strip_df <- strip_df %>% mutate(Type = factor(Type, levels = type_order))
  
  # Colour by lineage, matching broad_colours' logic: RG- subtypes
  # (progenitors) get shades of purple, EN- subtypes (excitatory
  # neurons) get shades of brown, and IPC-EN -- the actual
  # transitional cell type between the two -- gets a genuine
  # purple-brown blend rather than an arbitrary third colour.
  type_chr    <- as.character(type_order)
  prog_types  <- type_order[str_detect(type_chr, "^RG-")]
  ex_types    <- type_order[str_detect(type_chr, "^EN-")]
  trans_types <- type_order[type_chr == "IPC-EN"]
  
  purple_ramp <- if (length(prog_types) > 0)
    setNames(colorRampPalette(c("#3F007D", "#B8A6D9"))(length(prog_types)), prog_types)
  else character(0)
  brown_ramp <- if (length(ex_types) > 0)
    setNames(colorRampPalette(c("#4A2C0B", "#D9B382"))(length(ex_types)), ex_types)
  else character(0)
  blend_colour <- colorRampPalette(c("#6A51A3", "#8C510A"))(3)[2]
  trans_ramp <- if (length(trans_types) > 0)
    setNames(rep(blend_colour, length(trans_types)), trans_types)
  else character(0)
  
  fine_colours <- c(purple_ramp, brown_ramp, trans_ramp)[type_chr]
  
  ggplot(strip_df, aes(x = pseudotime, y = Type, colour = Type)) +
    geom_jitter(height = 0.3, size = 0.4, alpha = 0.15) +
    scale_colour_manual(values = fine_colours, guide = "none") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, name = x_title) +
    labs(y = NULL) +
    theme_classic(base_size = BASE_SIZE - 2) +
    theme(axis.text = element_text(size = LABEL_SIZE - 2.5),
          axis.title = element_text(size = TITLE_SIZE - 2))
}


age_breaks <- log(c(100, 200, 500, 1000, 2000, 5000))
age_labels <- c("100","200","500","1k","2k","5k")
p_banner_age <- make_proportion_banner(banner_age_df, "log_age", age_breaks, age_labels,
                                       "Postconceptional age (days)")
p_B <- make_trio(panelB_age_df, panelC_age_df, p_banner_age, "log_age", shared_ylim,
                 age_breaks, age_labels, tag = "b", bottom_height = 1.3)

pt_range  <- range(c(panelB_pt_df$pseudotime, panelC_pt_df$pseudotime), na.rm = TRUE)
pt_breaks <- pretty(pt_range, n = 5)
p_strip_pt <- make_pseudotime_strip(cell_meta, fine_pt_types, pt_breaks,
                                    as.character(pt_breaks), "Average pseudotime")
p_C <- make_trio(panelB_pt_df, panelC_pt_df, p_strip_pt, "pseudotime", shared_ylim,
                 pt_breaks, as.character(pt_breaks), tag = "c", bottom_height = 2.5)

ggsave(file.path(out_dir, "Fig5B_age_trio.pdf"), p_B, width = 7, height = 7)
ggsave(file.path(out_dir, "Fig5C_pseudotime_trio.pdf"), p_C, width = 7, height = 8)
ts("  Panels B/C done")

# ============================================================
# SECTION 8: Panel D -- family x broad-cell-type activity,
#            ALL significant FLAMES-overlapping eRegulons
# ============================================================
ts("Building Panel D (family x broad cell type)...")

aucell_raw <- read_excel(aucell_path)
colnames(aucell_raw) <- str_trim(colnames(aucell_raw))
ct_cols <- colnames(aucell_raw)[-(1:2)]

collapse_fine_celltype <- function(ct) {
  dplyr::case_when(
    str_detect(ct, "^RG")                     ~ "Progenitors",
    str_detect(ct, "^(IPC|Tri-IPC|Tri_IPC)")   ~ "Progenitors",
    str_detect(ct, "^EN")                      ~ "Excitatory neurons",
    str_detect(ct, "^IN")                      ~ "Inhibitory neurons",
    str_detect(ct, "Cajal")                    ~ "Cajal-Retzius",
    str_detect(ct, "^Astrocyte")               ~ "Astrocytes",
    str_detect(ct, "^OPC")                     ~ "OPC",
    str_detect(ct, "^Oligodendrocyte")         ~ "Oligodendrocyte",
    str_detect(ct, "^Microglia")               ~ "Microglia",
    TRUE                                       ~ NA_character_
  )
}
ct_to_broad <- tibble(cell_type = ct_cols, broad = collapse_fine_celltype(ct_cols)) %>%
  filter(!is.na(broad))

aucell_long <- aucell_raw %>%
  rename(eregulon = 1, tf_raw = 2) %>%
  select(eregulon, all_of(ct_to_broad$cell_type)) %>%
  pivot_longer(-eregulon, names_to = "cell_type", values_to = "aucell") %>%
  left_join(ct_to_broad, by = "cell_type") %>%
  filter(!is.na(aucell))

aucell_z <- aucell_long %>%
  group_by(eregulon) %>%
  mutate(z = if (sd(aucell, na.rm = TRUE) > 0)
    (aucell - mean(aucell, na.rm = TRUE)) / sd(aucell, na.rm = TRUE) else 0) %>%
  ungroup()

ereg_family_polarity <- enrich_results %>%
  # Panel D asks about eRegulon ACTIVITY -- a property of the
  # regulatory unit (TF + its target genes) -- not about the TF
  # gene's own locus. The criterion that actually tests relevance to
  # a phenotype is whether the eRegulon's target-gene set is itself
  # significantly enriched among FLAMES genes for that phenotype;
  # whether the TF's own gene additionally clears FLAMES significance
  # is a separate, largely orthogonal fact (many real regulators sit
  # upstream of GWAS-significant target programs without being
  # GWAS/FLAMES hits themselves). So this deliberately does NOT use
  # flames_tfs / tf_is_flames -- restricted to eRegulon significance
  # only, regardless of the TF's own FLAMES status.
filter(sig != "ns") %>%
  distinct(eregulon, tf, family, polarity)

# ---- DIAGNOSTIC: does Progenitor vs Astrocyte/OPC activity differ, and ----
# does that difference actually vary by family, or does the heatmap
# just LOOK similar across families because of the shared colour
# scale? Three checks, broad to specific:
#   1. Overall: pooling all significant eRegulons, is Progenitor
#      activity different from Astrocyte/OPC activity at all?
#   2. Per-family: does that same contrast hold up within each
#      family separately?
#   3. Cross-family: is the SIZE of the Progenitor-vs-glia difference
#      itself distinguishable between families (Kruskal-Wallis on the
#      per-eRegulon differences, grouped by family)? This is the
#      direct test of "are families really different from each other
#      on this contrast" -- (1) and (2) alone can't answer that, since
#      every family could individually show the same-sized effect.
ts("Running Progenitor vs Astrocyte/OPC diagnostic tests...")

# Per-eRegulon, per-broad-cell-type mean z (one row per significant
# eRegulon x broad type, before collapsing across eRegulons into the
# family-level heatmap means) -- this is the correct unit for pairing,
# since z was originally normalised per eRegulon across cell types.
ereg_broad_z <- aucell_z %>%
  inner_join(ereg_family_polarity, by = "eregulon", relationship = "many-to-many") %>%
  group_by(eregulon, tf, family, polarity, broad) %>%
  summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop")

ereg_wide <- ereg_broad_z %>%
  pivot_wider(id_cols = c(eregulon, tf, family, polarity),
              names_from = broad, values_from = mean_z)

run_paired_check <- function(df, glia_col, label) {
  d <- df %>% dplyr::filter(!is.na(Progenitors), !is.na(.data[[glia_col]]))
  if (nrow(d) < 3) return(tibble(comparison = label, n = nrow(d), median_diff = NA, p_value = NA))
  diff <- d[[glia_col]] - d$Progenitors
  wt <- suppressWarnings(wilcox.test(d[[glia_col]], d$Progenitors, paired = TRUE))
  tibble(comparison = label, n = nrow(d), median_diff = median(diff, na.rm = TRUE), p_value = wt$p.value)
}

# 1. Overall (pooled across families)
overall_check <- bind_rows(
  run_paired_check(ereg_wide, "Astrocytes", "Progenitors vs Astrocytes (overall)"),
  run_paired_check(ereg_wide, "OPC",        "Progenitors vs OPC (overall)"))
cat("\n-- Overall Progenitor vs glia contrast --\n"); print(overall_check)

# 2. Per-family
per_family_check <- purrr::map_dfr(unique(ereg_wide$family), function(fam) {
  d <- ereg_wide %>% dplyr::filter(family == fam)
  bind_rows(
    run_paired_check(d, "Astrocytes", paste(fam, "- Progenitors vs Astrocytes")),
    run_paired_check(d, "OPC",        paste(fam, "- Progenitors vs OPC")))
}) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))
cat("\n-- Per-family Progenitor vs glia contrast (FDR-corrected) --\n"); print(per_family_check)

# 3. Does the SIZE of the Progenitor-vs-glia difference vary by family?
cross_family_check <- purrr::map_dfr(c("Astrocytes", "OPC"), function(glia_col) {
  d <- ereg_wide %>% dplyr::filter(!is.na(Progenitors), !is.na(.data[[glia_col]])) %>%
    mutate(diff = .data[[glia_col]] - Progenitors)
  if (n_distinct(d$family) < 2 || nrow(d) < 5) return(NULL)
  kt <- kruskal.test(diff ~ family, data = d)
  tibble(glia_type = glia_col, n = nrow(d), n_families = n_distinct(d$family),
         kruskal_p = kt$p.value)
})
cat("\n-- Does the Progenitor-vs-glia difference vary BY family? --\n"); print(cross_family_check)
cat("(A non-significant kruskal_p here means: no, the families are not\n")
cat(" statistically distinguishable on this contrast -- i.e. the concern\n")
cat(" that results look similar across families would be confirmed, not\n")
cat(" just a visual artefact of the shared heatmap colour scale.)\n\n")

write.csv(ereg_wide, file.path(out_dir, "SupplTable_progenitor_glia_per_eregulon.csv"), row.names = FALSE)
write.csv(bind_rows(overall_check, per_family_check),
          file.path(out_dir, "SupplTable_progenitor_glia_paired_tests.csv"), row.names = FALSE)
write.csv(cross_family_check, file.path(out_dir, "SupplTable_progenitor_glia_cross_family.csv"), row.names = FALSE)


panelD_df <- aucell_z %>%
  inner_join(ereg_family_polarity, by = "eregulon", relationship = "many-to-many") %>%
  group_by(family, polarity, broad) %>%
  summarise(mean_z = mean(z, na.rm = TRUE), n_eregulons = n_distinct(eregulon), .groups = "drop") %>%
  mutate(family = factor(family, levels = family_order),
         broad   = factor(broad,  levels = broad_order))

write.csv(panelD_df, file.path(out_dir, "Fig5D_family_celltype_activity.csv"), row.names = FALSE)

make_panel_d_heatmap <- function(pol_label, tag = NULL) {
  df <- panelD_df %>% filter(polarity == pol_label)
  ggplot(df, aes(x = broad, y = family, fill = mean_z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, name = "Mean AUCell\n(z-score)",
                         limits = c(-2, 2), oob = squish) +
    scale_x_discrete(position = "top", labels = function(x) str_wrap(x, 10)) +
    labs(x = NULL, y = NULL, title = pol_label, tag = tag) +
    theme_bw(base_size = BASE_SIZE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = LABEL_SIZE - 1),
          axis.text.y = element_text(size = LABEL_SIZE),
          plot.title  = element_text(face = "bold", size = TITLE_SIZE),
          panel.grid  = element_blank())
}

p_D_act <- make_panel_d_heatmap("Activating", tag = "d")
p_D_rep <- make_panel_d_heatmap("Repressing")
p_D <- p_D_act + p_D_rep + plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "Fig5D_family_celltype_activity.pdf"), p_D, width = 11, height = 4)
ts("  Panel D done")

# ============================================================
# SECTION 9: Assembly
# Layout: A | B (age trio) | C (pseudotime trio)  -- one row
#         legends (family colour, cell-type colour)
#         D across the bottom
# ============================================================
ts("Assembling Figure 5...")

leg_df <- tibble(family = factor(family_order, levels = family_order), x = 1, y = seq_along(family_order))
p_leg <- ggplot(leg_df, aes(x = x, y = y, colour = family)) +
  geom_point(size = 4) +
  scale_colour_manual(values = family_colours, name = "Phenotype family",
                      guide = guide_legend(nrow = 2, override.aes = list(size = 4))) +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9))

banner_leg_df <- tibble(broad = factor(broad_order, levels = broad_order), x = 1, y = seq_along(broad_order))
p_banner_leg <- ggplot(banner_leg_df, aes(x = x, y = y, fill = broad)) +
  geom_col() +
  scale_fill_manual(values = broad_colours, name = "Cell type",
                    guide = guide_legend(nrow = 2)) +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8))

# Extract ONLY the legend guides -- not the dummy plot panels
# (fixes the giant colour-ribbon artifact from embedding raw geoms)
leg_family_grob <- extract_legend(p_leg)
leg_broad_grob  <- extract_legend(p_banner_leg)
legends_row <- wrap_elements(full = leg_family_grob) + wrap_elements(full = leg_broad_grob) +
  plot_layout(widths = c(1, 1))

top_row <- p_A + wrap_elements(full = p_B) + wrap_elements(full = p_C) +
  plot_layout(widths = c(1, 2, 2))

fig5 <- (top_row / legends_row / wrap_elements(full = p_D)) +
  plot_layout(heights = c(6, 0.6, 3))

ggsave(file.path(out_dir, "Figure5_v3.pdf"), fig5, width = 15, height = 10)
ggsave(file.path(out_dir, "Figure5_v3.png"), fig5, width = 15, height = 10, dpi = 300)

ts("Done. Outputs in:", out_dir)