# ===========================================================================
# Siletti_FLAMES_enrichment.R
# Manuscript item: Supplementary Figure 6 (adult brain cell types)
#
# What this does:
#   Repeats the cell-type enrichment in the adult brain (Siletti atlas)
#   using the published expression specificity scores, with and without the
#   outlier phenotypes.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# siletti_flames_enrichment.R
#
# Cell type enrichment of FLAMES-prioritised genes in Siletti
# 2023 adult human brain atlas, using published expression
# specificity scores from Yao et al. 2025 (Nat Commun).
#
# Model per phenotype x supercluster:
#   is_FLAMES_gene ~ specificity + log_length
#
# where specificity = EP_g,c / sum_c(EP_g,c) from Yao et al.
# (expression proportion, equivalent to log_rel_expr in Wang
# pipeline but on a 0-1 scale rather than log-ratio scale).
#
# Background = protein-coding genes in Siletti pseudobulk
# (same 15k genes as specificity file).
#
# Output: enrich_siletti_full.rds and enrich_siletti_nooutlier.rds
# in cache/, same format as enrich_full/enrich_noout from Wang.
# Figure output: SupplementaryFigure2.pdf (+ _nooutlier sensitivity
# variant) in Outputs/Figure_4/.

source("config.R")

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)

ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# ---- Paths ----
data_dir    <- file.path(proj_dir, "Data")
base_dir    <- file.path(proj_dir, "Outputs")
cache_dir   <- file.path(base_dir, "cache")
flames_path <- file.path(proj_dir, "Data", "FLAMES_gene_prioritisation_results.txt")
out_dir    <- file.path(proj_dir, "Outputs/Figure_4")


spec_path   <- file.path(data_dir,
                         "Siletti_Supercluster_expression_specificity_TDEP_label.tsv")

dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)

# ---- Cell types and phenotypes to exclude ----
exclude_phenos <- c("SubCorVol", "SubCor", "SC_GMV")  # matches Wang pipeline
exclude_types  <- c("Vascular", "Fibroblast", "Bergmann_glia",
                    "Choroid_plexus", "Ependymal",
                    "Miscellaneous", "Splatter")

# ---- Phenotype family definitions (from Figure 4 script) ----
family_map <- list(
  "Brain size"                  = c("SA","FI","IC","GMV","SubCorVol"),
  "Cortical thickness"          = c("CT"),
  "Curvature"                   = c("GC","MC"),
  "Microstructural coherence"   = c("FA","ICVF"),
  "Microstructural diffusivity" = c("MD","ISOVF"),
  "Orientation dispersion"      = c("OD")
)
family_order <- c(
  "Brain size","Curvature","Cortical thickness",
  "Microstructural diffusivity",
  "Microstructural coherence","Orientation dispersion")
pheno_family_vec <- setNames(
  rep(names(family_map), sapply(family_map, length)),
  unlist(family_map))

outlier_traits <- c(
  "regional_FA_RSC","regional_FI_33pr","regional_FI_PeEc",
  "regional_FI_Pir","regional_FI_a24pr","regional_GC_33pr",
  "regional_GC_PHA1","regional_GC_PeEc","regional_GC_Pir",
  "regional_GC_PreS","regional_GC_RSC","regional_GC_a24pr",
  "regional_GC_p24","regional_GMV_EC","regional_GMV_PHA1",
  "regional_GMV_PHA2","regional_GMV_PeEc","regional_IC_Pir",
  "regional_IC_a24pr","regional_MC_33pr","regional_MC_PHA1",
  "regional_MC_PHA3","regional_MC_PeEc","regional_MC_Pir",
  "regional_MC_PreS","regional_MC_RSC","regional_MC_a24pr",
  "regional_MC_p24","regional_CT_PHA1","SubCorVol_CC_Central",
  "SubCorVol_CC_Mid_Anterior","SubCorVol_CC_Mid_Posterior",
  "regional_FI_FOP3","regional_FI_PoI1","regional_GC_4",
  "regional_GC_52","regional_GC_AAIC","regional_GC_AVI",
  "regional_GC_FOP2","regional_GC_FOP3","regional_GC_Ig",
  "regional_GC_LBelt","regional_GC_MBelt","regional_GC_MI",
  "regional_GC_PI","regional_GC_PoI1","regional_GC_PoI2",
  "regional_GC_6d","regional_IC_AVI","regional_IC_FOP3",
  "regional_IC_MI","regional_IC_PoI1","regional_MC_52",
  "regional_MC_AAIC","regional_MC_FOP3","regional_MC_Ig",
  "regional_MC_MBelt","regional_MC_MI","regional_MC_PI",
  "regional_MC_PoI1","regional_MC_PoI2")

# ================================================================
# 1. Load and prepare specificity data
# ================================================================
ts("Loading Siletti specificity scores...")
spec_raw <- read.table(spec_path, header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
cat(sprintf("Specificity file: %d rows, columns: %s\n",
            nrow(spec_raw),
            paste(colnames(spec_raw), collapse=", ")))

# Rename columns to standard names
spec_raw <- spec_raw %>%
  dplyr::rename(
    ensembl_id   = ENSGID,
    supercluster = Supercluster,
    specificity  = specificity)

# ---- Map Ensembl ID -> gene symbol via org.Hs.eg.db ----
# (biomaRt not used - Ensembl may be unavailable)
ts("Mapping Ensembl IDs to gene symbols via org.Hs.eg.db...")
id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys   = unique(spec_raw$ensembl_id),
  columns= c("ENSEMBL","SYMBOL"),
  keytype= "ENSEMBL") %>%
  dplyr::rename(ensembl_id  = ENSEMBL,
                gene_symbol = SYMBOL) %>%
  dplyr::filter(!is.na(gene_symbol)) %>%
  # Keep one symbol per Ensembl ID (first match)
  dplyr::distinct(ensembl_id, .keep_all=TRUE)

cat(sprintf("Mapped %d / %d Ensembl IDs to symbols\n",
            nrow(id_map), n_distinct(spec_raw$ensembl_id)))

spec_raw <- spec_raw %>%
  dplyr::left_join(id_map, by="ensembl_id") %>%
  dplyr::mutate(gene_symbol=as.character(gene_symbol)) %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "")
cat(sprintf("After symbol mapping: %d rows, %d genes\n",
            nrow(spec_raw),
            n_distinct(spec_raw$gene_symbol)))

# ---- Filter excluded cell types ----
spec_raw <- spec_raw %>%
  dplyr::filter(!supercluster %in% exclude_types)
cat(sprintf("Cell types after exclusion: %s\n",
            paste(sort(unique(spec_raw$supercluster)),
                  collapse=", ")))

# Background = all protein-coding genes in the specificity file
pc_pool      <- unique(spec_raw$gene_symbol)
cell_types   <- sort(unique(spec_raw$supercluster))
cat(sprintf("Background: %d genes, %d cell types\n",
            length(pc_pool), length(cell_types)))

# ---- Wide format: genes x cell types ----
# For duplicate gene symbols, keep highest specificity per cell type
spec_wide <- spec_raw %>%
  dplyr::group_by(gene_symbol, supercluster) %>%
  dplyr::slice_max(specificity, n=1, with_ties=FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene_symbol, supercluster, specificity) %>%
  tidyr::pivot_wider(names_from=supercluster,
                     values_from=specificity,
                     values_fill=0)

# Explicitly separate gene name column from cell type columns
gene_syms_wide <- as.character(spec_wide$gene_symbol)
ct_cols        <- setdiff(colnames(spec_wide), "gene_symbol")
spec_mat       <- as.matrix(spec_wide[, ct_cols])
rownames(spec_mat) <- gene_syms_wide
cat(sprintf("Specificity matrix: %d genes x %d cell types\n",
            nrow(spec_mat), ncol(spec_mat)))

# ================================================================
# 2. Gene lengths
# ================================================================
gene_length_cache <- file.path(cache_dir, "gene_lengths_hg19.rds")
if (file.exists(gene_length_cache)) {
  ts("Loading gene lengths from cache...")
  gene_lengths_df <- readRDS(gene_length_cache)
} else {
  ts("Querying Ensembl GRCh37 for gene lengths...")
  gene_lengths_df <- tryCatch({
    mart37 <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                               dataset="hsapiens_gene_ensembl",
                               host="https://grch37.ensembl.org")
    gene_lengths_raw <- biomaRt::getBM(
      attributes=c("hgnc_symbol","transcript_length"),
      filters="hgnc_symbol", values=pc_pool, mart=mart37)
    gene_lengths_raw %>%
      dplyr::filter(hgnc_symbol != "") %>%
      dplyr::group_by(hgnc_symbol) %>%
      dplyr::summarise(gene_length=max(transcript_length, na.rm=TRUE),
                       .groups="drop") %>%
      dplyr::rename(symbol=hgnc_symbol)
  }, error=function(e) {
    ts("biomaRt unavailable — skipping gene length covariate")
    data.frame(symbol=character(0), gene_length=integer(0))
  })
  if (nrow(gene_lengths_df) > 0)
    saveRDS(gene_lengths_df, gene_length_cache)
}

# Gene-level covariates
gene_covars <- data.frame(
  symbol=rownames(spec_mat),
  stringsAsFactors=FALSE) %>%
  dplyr::left_join(gene_lengths_df, by="symbol") %>%
  dplyr::mutate(log_length=log10(gene_length + 1))
use_length <- !all(is.na(gene_covars$log_length))
cat(sprintf("Gene lengths available: %d / %d genes\n",
            sum(!is.na(gene_covars$log_length)),
            nrow(gene_covars)))

# ================================================================
# 3. Load and prepare FLAMES data
# ================================================================
ts("Loading FLAMES data...")
df_raw <- read.table(flames_path, header=TRUE,
                     stringsAsFactors=FALSE)

prep_flames <- function(df_raw, exclude_outliers=FALSE) {
  df <- df_raw %>%
    dplyr::mutate(
      phenotype_clean = phenotype %>%
        str_remove("^global_") %>% str_remove("^regional_"),
      phenotype_type = dplyr::case_when(
        str_detect(phenotype_clean,"^tract_") ~
          str_extract(phenotype_clean,"(?<=^tract_)[A-Za-z0-9]+"),
        TRUE ~ str_split(phenotype_clean,"_") %>%
          purrr::map_chr(1))) %>%
    dplyr::filter(phenotype_type %in% names(pheno_family_vec)) %>%
    # BUGFIX: exclude_phenos (SubCorVol/SubCor/SC_GMV) used to be
    # applied only inside run_siletti_enrichment() -- AFTER the
    # cache-hit early return, so it silently never took effect once
    # enrich_siletti_full.rds / enrich_siletti_nooutlier.rds existed
    # on disk. Filtering here, at the shared data-prep stage, means
    # SubCorVol is removed from df_full/df_noout before either the
    # cache path or the fresh-compute path ever sees it, so both are
    # guaranteed consistent regardless of cache state.
    dplyr::filter(!phenotype_type %in% exclude_phenos)
  if (exclude_outliers) {
    df <- df %>% dplyr::filter(!phenotype %in% outlier_traits)
  }
  df
}

df_full  <- prep_flames(df_raw, exclude_outliers=FALSE)
df_noout <- prep_flames(df_raw, exclude_outliers=TRUE)
cat(sprintf("FLAMES full: %d rows | nooutlier: %d rows\n",
            nrow(df_full), nrow(df_noout)))

# ================================================================
# 4. Run enrichment
# ================================================================
run_siletti_enrichment <- function(df_in, suffix) {
  cache_file <- file.path(cache_dir,
                          sprintf("enrich_siletti_%s.rds", suffix))
  if (file.exists(cache_file)) {
    ts(sprintf("  Loading cached: %s", suffix))
    # BUGFIX: defensively re-apply exclude_phenos even to cached
    # results, in case the cache on disk was written by an older
    # version of this script (e.g. before df_in was pre-filtered in
    # prep_flames()). This is what was silently missing before --
    # a stale enrich_siletti_full.rds kept leaking SubCorVol rows
    # into every downstream table and figure.
    return(readRDS(cache_file) %>%
             dplyr::filter(!phenotype_type %in% exclude_phenos))
  }
  ts(sprintf("  Running Siletti enrichment: %s", suffix))

  # FLAMES gene sets per phenotype
  # (df_in is already filtered for exclude_phenos in prep_flames(),
  # this second filter is kept as a harmless defensive no-op)
  genes_by_pheno <- df_in %>%
    dplyr::filter(!phenotype_type %in% exclude_phenos) %>%
    dplyr::distinct(symbol, phenotype_type) %>%
    dplyr::group_by(phenotype_type) %>%
    dplyr::summarise(genes=list(sort(unique(symbol))),
                     .groups="drop") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pc_genes = list(intersect(genes, pc_pool)),
      n_pc     = length(pc_genes)) %>%
    dplyr::filter(n_pc >= 5) %>%
    dplyr::ungroup()
  
  cat(sprintf("  Phenotypes with >=5 FLAMES genes: %d\n",
              nrow(genes_by_pheno)))
  
  results <- lapply(seq_len(nrow(genes_by_pheno)), function(i) {
    ptype  <- genes_by_pheno$phenotype_type[i]
    pgenes <- genes_by_pheno$pc_genes[[i]]
    
    lapply(cell_types, function(ct) {
      spec_ct <- spec_mat[, ct]
      
      mdf <- gene_covars %>%
        dplyr::mutate(
          is_phenotype = as.integer(symbol %in% pgenes),
          specificity  = spec_ct[symbol])
      
      if (use_length) {
        mdf        <- mdf %>% dplyr::filter(!is.na(log_length))
        formula_lm <- is_phenotype ~ specificity + log_length
      } else {
        formula_lm <- is_phenotype ~ specificity
      }
      
      fit <- tryCatch(
        glm(formula_lm, data=mdf, family=binomial()),
        error=function(e) NULL)
      if (is.null(fit)) return(NULL)
      
      cf <- tryCatch(summary(fit)$coefficients,
                     error=function(e) NULL)
      if (is.null(cf) || !"specificity" %in% rownames(cf)) {
        return(NULL)
      }
      
      data.frame(
        phenotype_type = ptype,
        celltype       = ct,
        # No timepoint for adult Siletti - set to "Adult"
        timepoint      = "Adult",
        beta           = cf["specificity","Estimate"],
        se             = cf["specificity","Std. Error"],
        p_value        = cf["specificity","Pr(>|z|)"],
        stringsAsFactors=FALSE)
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  # FDR correction within phenotype
  results <- results %>%
    dplyr::group_by(phenotype_type) %>%
    dplyr::mutate(p_fdr = p.adjust(p_value, method="BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(significant = p_fdr < 0.05 & beta > 0)
  
  saveRDS(results, cache_file)
  ts(sprintf("  Done: %d rows, %d significant bins",
             nrow(results), sum(results$significant, na.rm=TRUE)))
  results
}

enrich_siletti_full  <- run_siletti_enrichment(df_full,  "full")
enrich_siletti_noout <- run_siletti_enrichment(df_noout, "nooutlier")

# ================================================================
# 5. Load MAGMA Siletti data
# ================================================================
ts("Loading MAGMA Siletti results...")
magma_siletti_raw <- read.csv(
  file.path(data_dir, "MAGMA_Siletti.csv"),
  header=TRUE, sep=",", stringsAsFactors=FALSE)

magma_siletti_raw <- magma_siletti_raw %>%
  dplyr::filter(!cell_type %in% exclude_types)
cat(sprintf("Cell types after exclusion: %s\n",
            paste(sort(unique(magma_siletti_raw$cell_type)),
                  collapse=", ")))

magma_siletti_raw <- magma_siletti_raw %>%
  dplyr::filter(!phenotype %in% exclude_phenos)
# BUGFIX: this message previously said "Cell types after exclusion"
# (copy-pasted from the block above) while actually printing the
# *phenotypes* remaining after the exclude_phenos filter -- fixed.
cat(sprintf("Phenotypes after exclusion: %s\n",
            paste(sort(unique(magma_siletti_raw$phenotype)),
                  collapse=", ")))


# Recompute FDR within phenotype, keep beta > 0 only
magma_siletti <- magma_siletti_raw %>%
  dplyr::rename(phenotype_type = phenotype,
                celltype       = cell_type) %>%
  dplyr::group_by(phenotype_type) %>%
  dplyr::mutate(fdr_magma = p.adjust(p, method="BH")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(beta > 0, fdr_magma < 0.05) %>%
  dplyr::distinct(phenotype_type, celltype)

cat(sprintf("MAGMA Siletti significant (beta>0, FDR<0.05): %d rows\n",
            nrow(magma_siletti)))
print(magma_siletti)

# ================================================================
# 6. Supplementary Figure 2
# ================================================================
ts("Building Supplementary Figure 2 (Siletti FLAMES + MAGMA)...")

library(ggplot2)
library(cowplot)

# ---- Colour/order helpers ----
family_colours <- c(
  "Brain size"                  = "#1F78B4",
  "Cortical thickness"          = "#33A02C",
  "Curvature"                   = "#FF7F00",
  "Microstructural coherence"   = "#C0006E",
  "Microstructural diffusivity" = "#6A3D9A",
  "Orientation dispersion"      = "#B15928")

# Siletti supercluster biological order
ct_hierarchy_siletti <- c(
  "Upper_layer_intratelencephalic",
  "Deep_layer_intratelencephalic",
  "Deep_layer_near_projecting",
  "Deep_layer_corticothalamic_and_6b",
  "Hippocampal_CA1_3",
  "Hippocampal_CA4",
  "Hippocampal_dentate_gyrus",
  "Amygdala_excitatory",
  "Thalamic_excitatory",
  "Mammillary_body",
  "MGE_interneuron",
  "CGE_interneuron",
  "LAMP5_LHX6_and_Chandelier",
  "Medium_spiny_neuron",
  "Eccentric_medium_spiny_neuron",
  "Midbrain_derived_inhibitory",
  "Lower_rhombic_lip",
  "Upper_rhombic_lip",
  "Cerebellar_inhibitory",
  "Astrocyte",
  "Microglia",
  "Oligodendrocyte",
  "Oligodendrocyte_precursor",
  "Committed_oligodendrocyte_precursor")

make_suppfig2_panel <- function(enrich_res, suffix) {
  sig <- enrich_res %>%
    dplyr::filter(significant) %>%
    dplyr::mutate(family = pheno_family_vec[phenotype_type],
                  family = factor(family, levels=family_order))
  
  if (nrow(sig) == 0) {
    cat(sprintf("No significant FLAMES enrichments: %s\n", suffix))
    return(NULL)
  }
  
  # Cell type order
  present_cts    <- unique(sig$celltype)
  celltype_order <- c(ct_hierarchy_siletti[
    ct_hierarchy_siletti %in% present_cts],
    sort(setdiff(present_cts,
                 ct_hierarchy_siletti)))
  
  # x-axis: phenotypes grouped by family
  x_order <- data.frame(
    phenotype_type = names(pheno_family_vec),
    family         = pheno_family_vec,
    stringsAsFactors=FALSE) %>%
    dplyr::filter(phenotype_type %in% unique(sig$phenotype_type)) %>%
    dplyr::mutate(family=factor(family, levels=family_order)) %>%
    dplyr::arrange(family, phenotype_type) %>%
    dplyr::pull(phenotype_type)
  
  fam_sep_x <- data.frame(
    phenotype_type=names(pheno_family_vec),
    family=pheno_family_vec, stringsAsFactors=FALSE) %>%
    dplyr::filter(phenotype_type %in% x_order) %>%
    dplyr::mutate(family=factor(family, levels=family_order),
                  xpos=match(phenotype_type, x_order)) %>%
    dplyr::group_by(family) %>%
    dplyr::summarise(x_end=max(xpos)+0.5, .groups="drop") %>%
    dplyr::pull(x_end)
  fam_sep_x <- fam_sep_x[-length(fam_sep_x)]
  
  plot_df <- sig %>%
    dplyr::mutate(
      celltype       = factor(celltype, levels=celltype_order),
      phenotype_type = factor(phenotype_type, levels=x_order))
  
  # MAGMA inner circle: cell type x phenotype significant
  magma_inner <- magma_siletti %>%
    dplyr::filter(celltype       %in% celltype_order,
                  phenotype_type %in% x_order) %>%
    dplyr::mutate(
      celltype       = factor(celltype,       levels=celltype_order),
      phenotype_type = factor(phenotype_type, levels=x_order))
  
  n_ct <- length(celltype_order)
  n_x  <- length(x_order)
  
  p <- ggplot(plot_df,
              aes(x=phenotype_type, y=celltype)) +
    geom_vline(xintercept=fam_sep_x,
               colour="grey60", linewidth=0.4, linetype="dashed") +
    # FLAMES: hollow circle (fixed size, no timepoint dimension)
    geom_point(size=4.5, shape=21,
               fill="white", colour="grey30",
               stroke=1.4, alpha=0.95) +
    # MAGMA inner filled circle (dark grey)
    { if (nrow(magma_inner) > 0)
      geom_point(data=magma_inner,
                 aes(x=phenotype_type, y=celltype),
                 size=2.0, shape=21,
                 fill="#444444", colour="#444444",
                 stroke=0, inherit.aes=FALSE)
      else list() } +
    scale_x_discrete(expand=expansion(add=0.4)) +
    scale_y_discrete(expand=expansion(add=0.4)) +
    labs(x="Phenotype (grouped by family)", y="Cell type") +
    theme_bw(base_size=12) +
    theme(
      axis.text.x      = element_text(angle=45, hjust=1,
                                      size=10, colour="black"),
      axis.text.y      = element_text(size=10, colour="black"),
      axis.title       = element_text(size=11),
      panel.grid.major = element_line(colour="grey94"),
      plot.margin      = margin(t=10, r=15, b=5, l=5),
      legend.position  = "none")
  
  # Legend
  leg_df <- data.frame(
    x     = c(1, 2),
    y     = c(1, 1),
    label = c("FLAMES enrichment\n(FDR<0.05, β>0)",
              "MAGMA enrichment\n(FDR<0.05, β>0)"),
    stringsAsFactors=FALSE)
  
  p_leg <- ggplot(leg_df, aes(x=x, y=y)) +
    geom_point(data=leg_df[1,],
               size=6, shape=21,
               fill="white", colour="grey30", stroke=1.4) +
    geom_point(data=leg_df[2,],
               size=2.5, shape=21,
               fill="#444444", colour="#444444", stroke=0) +
    geom_text(aes(label=label), vjust=-2.0, size=3.5,
              colour="grey20") +
    scale_x_continuous(limits=c(0.5, 2.5)) +
    scale_y_continuous(limits=c(0.5, 1.8)) +
    theme_void()
  
  p_final <- cowplot::plot_grid(
    p, p_leg, ncol=1, rel_heights=c(1, 0.08))
  
  h <- max(5, n_ct * 0.32 + 2.5)
  w <- max(6, n_x  * 0.60 + 3.0)
  # Main dataset ("full") is the manuscript's Supplementary Figure 2;
  # the outlier-excluded sensitivity check gets a suffixed filename.
  fig_name <- if (suffix == "full") "SupplementaryFigure2.pdf" else
    sprintf("SupplementaryFigure2_%s.pdf", suffix)
  out_path <- file.path(out_dir, fig_name)
  ggsave(out_path, p_final, width=w, height=h,
         device=pdf, limitsize=FALSE)
  ts(sprintf("  Saved: %s", out_path))
  p_final
}

suppfig2_full  <- make_suppfig2_panel(enrich_siletti_full,  "full")
suppfig2_noout <- make_suppfig2_panel(enrich_siletti_noout, "nooutlier")

# ================================================================
# 7. Supplementary tables
# ================================================================
ts("Writing supplementary tables...")

table_dir <- file.path(out_dir, "tables")
dir.create(table_dir, recursive=TRUE, showWarnings=FALSE)

# ---- FLAMES Siletti enrichment table ----
# All tested results (not just significant) for full dataset,
# after exclusions (SubCorVol/SubCor phenotypes removed,
# Vascular/Fibroblast/Bergmann glia/Choroid plexus/Ependymal/
# Miscellaneous/Splatter cell types removed).
flames_suppl <- enrich_siletti_full %>%
  dplyr::mutate(
    family = pheno_family_vec[phenotype_type]) %>%
  dplyr::select(
    Family         = family,
    Phenotype      = phenotype_type,
    Cell_type      = celltype,
    Beta           = beta,
    SE             = se,
    P_value        = p_value,
    FDR            = p_fdr,
    Significant    = significant) %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  dplyr::arrange(Family, Phenotype, FDR)

write.csv(flames_suppl,
          file.path(table_dir,
                    "SupplTable_Siletti_FLAMES_enrichment.csv"),
          row.names=FALSE)
cat(sprintf("  FLAMES table: %d rows (%d significant)\n",
            nrow(flames_suppl),
            sum(flames_suppl$Significant, na.rm=TRUE)))

# ---- MAGMA Siletti table ----
# All results after exclusions (beta > 0 filter NOT applied here
# so the full picture is in the table; significance column
# reflects beta > 0 AND FDR < 0.05).
magma_suppl <- magma_siletti_raw %>%
  dplyr::rename(phenotype_type = phenotype,
                celltype       = cell_type) %>%
  # Apply same cell type exclusions as FLAMES
  dplyr::filter(!celltype %in% exclude_types) %>%
  # Apply same phenotype exclusions
  dplyr::filter(!phenotype_type %in% exclude_phenos) %>%
  # Recompute FDR within phenotype on filtered set
  dplyr::group_by(phenotype_type) %>%
  dplyr::mutate(FDR = p.adjust(p, method="BH")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Significant = beta > 0 & FDR < 0.05,
    family      = pheno_family_vec[phenotype_type]) %>%
  dplyr::select(
    Family         = family,
    Phenotype      = phenotype_type,
    Cell_type      = celltype,
    Beta_raw       = beta_raw,
    Beta           = beta,
    SE,
    P_value        = p,
    FDR,
    Significant) %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  dplyr::arrange(Family, Phenotype, FDR)

write.csv(magma_suppl,
          file.path(table_dir,
                    "SupplTable_Siletti_MAGMA_enrichment.csv"),
          row.names=FALSE)
cat(sprintf("  MAGMA table: %d rows (%d significant)\n",
            nrow(magma_suppl),
            sum(magma_suppl$Significant, na.rm=TRUE)))

# ================================================================
# 8. Summary
# (BUGFIX: this section was previously mislabelled "7." again,
# duplicating the "7. Supplementary tables" heading above)
# ================================================================
ts("All done.")
cat(strrep("=", 55), "\n")
cat("Significant FLAMES enrichments (FDR<0.05, beta>0):\n\n")
for (suffix in c("full","nooutlier")) {
  res <- if (suffix=="full") enrich_siletti_full else
    enrich_siletti_noout
  cat(sprintf("  [%s]\n", suffix))
  sig <- res %>%
    dplyr::filter(significant) %>%
    dplyr::count(celltype, sort=TRUE)
  for (i in seq_len(nrow(sig)))
    cat(sprintf("    %-45s %d phenotypes\n",
                sig$celltype[i], sig$n[i]))
  cat("\n")
}
cat(strrep("=", 55), "\n")
cat("Outputs in", out_dir, ":\n")
cat("  SupplementaryFigure2.pdf\n")
cat("  SupplementaryFigure2_nooutlier.pdf\n")
cat("Tables in", table_dir, ":\n")
cat("  SupplTable_Siletti_FLAMES_enrichment.csv\n")
cat("  SupplTable_Siletti_MAGMA_enrichment.csv\n")
cat("Caches:\n")
# BUGFIX: filenames below previously said "enrich_siletti_noout.rds"
# but the actual cache written by run_siletti_enrichment() (see
# sprintf("enrich_siletti_%s.rds", suffix), suffix="nooutlier") is
# enrich_siletti_nooutlier.rds -- corrected to match.
cat("  enrich_siletti_full.rds\n")
cat("  enrich_siletti_nooutlier.rds\n")