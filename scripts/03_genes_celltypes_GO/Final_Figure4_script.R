# ===========================================================================
# Final_Figure4_script.R
# Manuscript item: Figure 3 and Supplementary Figure 5 (script is named Figure 4)
#
# What this does:
#   Main FLAMES figure: unique vs shared genes per family with chord
#   diagram, cell-type enrichment in the developing cortex (Wang et al.),
#   and GO:BP enrichment, plus the within/between-family similarity
#   supplements. Cell-type model per phenotype and cell type x time bin:
#   FLAMES hit ~ log relative expression + log gene length.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Figure 4: Combined FLAMES figure - FINAL SCRIPT
#
# MAIN FIGURE (Figure4_full.pdf)
#   Panel A: Bar chart (unique vs shared genes per family) +
#            chord diagram of between-family gene overlap.
#   Panel B: Cell type enrichment - concentric hollow circles.
#            One circle per significant timepoint per cell type
#            x phenotype. Size and colour = timepoint. No fill.
#            All axis labels in black. No family labels on top.
#   Panel C: GO BP enrichment dotplots - top 5 terms per family,
#            2 rows x 3 columns. All labels in black.
#
# SUPPLEMENTARY FIGURE 1 (Supp1_hclust_full.pdf)
#   Three dendrograms: FLAMES Jaccard | CT enrichment | GO
#   Phenotype labels coloured by family.
#
# SUPPLEMENTARY FIGURE 2 (Supp2_family_similarity_violins.pdf)
#   Violins: within vs between family similarity.
#   Permutation p annotated (family labels shuffled; see rationale
#   in compute_family_similarity for why Wilcoxon is not used).
#
# SUPPLEMENTARY FIGURE 3 (Supp3_family_consistency.pdf)
#   Grouped bar chart: mean within-family similarity per family.
#
# ---- STATISTICAL MODEL (cell type enrichment) ----
# For each phenotype p and cell type x timepoint bin b:
#
#   is_FLAMES_gene_p ~ log_rel_expr_b + log_length
#
#   is_FLAMES_gene_p = 1 if FLAMES hit for phenotype p, else 0
#                      (background = all protein-coding genes
#                       in Wang 2025 pseudobulk)
#
#   log_rel_expr_b   = log10(expr_b + 0.01)
#                    - log10(mean_expr_all_bins + 0.01)
#                    Bin-specific expression relative to each
#                    gene's own cross-bin mean. Prevents
#                    inflation by broadly expressed genes.
#
#   log_length       = log10(longest_transcript + 1)
#                    Controls for gene length bias in FLAMES.
#                    From Ensembl GRCh37 via biomaRt.
#
#   beta = log-odds of FLAMES hit per unit relative expression.
#   ONE-SIDED test (H1: beta > 0, positive enrichment only):
#   p = pnorm(beta/se, lower.tail=FALSE), FDR-corrected within
#   phenotype (BH). Significant = FDR < 0.05. Using the one-sided
#   p directly (rather than a two-sided Wald p filtered post-hoc
#   by sign) spends the FDR budget only on the positive-enrichment
#   hypothesis; negative-beta bins are never recorded as hits.
#
# ---- GO ENRICHMENT ----
# ORA via clusterProfiler::enrichGO (BP ontology). Hypergeometric
# over-representation test - inherently one-sided (tests only for
# enrichment, not depletion), so no separate sign filter is needed.
# Background = protein-coding genes in Wang 2025 pseudobulk.
# minGSSize=10, maxGSSize=500, pvalueCutoff=1 (all terms kept).
# simplify() run only on FDR<0.05 terms (cutoff=0.7).
# ALL terms stored in go_all cache; simplified significant terms
# stored separately for panel C display.
# GO similarity uses all tested terms (not just significant) with
# Spearman rank correlation between -log10(p.adjust) profiles.
#
# ---- JACCARD SIMILARITY ----
# J(p_i, p_j) = |genes_i  genes_j| / |genes_i  genes_j|
#
# ---- PAIRED SIMILARITY MATRICES ----
# CT: paired using intersection of bins that converged in both.
# GO: paired using union of all terms across both datasets,
#     with Spearman correlation on full ranked -log10(p) profiles.
#
# ---- FAMILY STRUCTURE PERMUTATION TEST ----
# Observed delta  = mean(within-family) - mean(between-family).
# Null: shuffle family labels n_perm=10,000 times.
# p = proportion of null delta  >= observed delta  (one-tailed).
#
# ---- OUTPUTS ----
# Figures:  Outputs/Figure_4/
# Tables:   Outputs/Figure_4/tables/
# Cache:    Outputs/cache/

source("config.R")

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(circlize)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(biomaRt)
library(ggdendro)
library(readxl)

select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange
rename  <- dplyr::rename

# ---- Paths ----
flames_path <- file.path(proj_dir, "Data", "FLAMES_gene_prioritisation_results.txt")
base_dir    <- file.path(proj_dir, "Outputs")
out_dir     <- file.path(base_dir, "Figure_4")
table_dir   <- file.path(out_dir,  "tables")
cache_dir   <- file.path(base_dir, "cache")
dir.create(out_dir,   recursive=TRUE, showWarnings=FALSE)
dir.create(table_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)

ts      <- function(...) cat(format(Sys.time(),"[%H:%M:%S]"), ..., "\n")
set.seed(42)
n_perm   <- 10000
top_n_go <- 5

# ================================================================
# WANG MAGMA DATA
# ================================================================
# Global GWAS MAGMA cell type enrichment (Wang 2025 neocortex CEPO).
# One row per phenotype x cell type. FDR computed within phenotype.
# Used in Panel B as a filled inner circle at FDR < 0.05.

wang_magma_path <- file.path(
  proj_dir,
  "Data/Wang_celltype_cleaned_MAGMA.xlsx")

ts("Loading Wang MAGMA data...")
wang_magma_raw <- readxl::read_excel(wang_magma_path)

# FDR correction within each phenotype
wang_magma <- wang_magma_raw %>%
  dplyr::group_by(phenotype) %>%
  dplyr::mutate(fdr_magma = p.adjust(p, method="BH")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(fdr_magma < 0.05)

# ---- Cell type name harmonisation ----
# Wang MAGMA uses underscores; FLAMES enrichment uses hyphens/dots.
# Map Wang MAGMA cell_type -> FLAMES celltype label.
wang_ct_map <- c(
  "RG_vRG"                  = "RG-vRG",
  "RG_oRG"                  = "RG-oRG",
  "RG_tRG"                  = "RG-tRG",
  "IPC_EN"                  = "IPC-EN",
  "Tri_IPC"                 = "Tri-IPC",
  "EN_IT_Immature"          = "EN-IT-Immature",
  "EN_L2_3_IT"              = "EN-L2_3-IT",
  "EN_L4_IT"                = "EN-L4-IT",
  "EN_L5_6_NP"              = "EN-L5_6-NP",
  "EN_L5_ET"                = "EN-L5-ET",
  "EN_L5_IT"                = "EN-L5-IT",
  "EN_L6_CT"                = "EN-L6-CT",
  "EN_L6_IT"                = "EN-L6-IT",
  "EN_L6b"                  = "EN-L6b",
  "EN_Newborn"              = "EN-Newborn",
  "EN_Non_IT_Immature"      = "EN-Non-IT-Immature",
  "IN_CGE_Immature"         = "IN-CGE-Immature",
  "IN_CGE_SNCG"             = "IN-CGE-SNCG",
  "IN_CGE_VIP"              = "IN-CGE-VIP",
  "IN_MGE_Immature"         = "IN-MGE-Immature",
  "IN_MGE_PV"               = "IN-MGE-PV",
  "IN_MGE_SST"              = "IN-MGE-SST",
  "IN_NCx_dGE_Immature"     = "IN-NCx_dGE-Immature",
  "Cajal_Retzius_cell"      = "Cajal-Retzius cell",
  "Astrocyte_Fibrous"       = "Astrocyte-Fibrous",
  "Astrocyte_Immature"      = "Astrocyte-Immature",
  "Astrocyte_Protoplasmic"  = "Astrocyte-Protoplasmic",
  "Oligodendrocyte"         = "Oligodendrocyte",
  "OPC"                     = "OPC"
)

# Build a lookup: FLAMES celltype x phenotype_type -> wang_sig (TRUE/FALSE)
wang_sig_lookup <- wang_magma %>%
  dplyr::mutate(
    celltype_flames = wang_ct_map[cell_type],
    phenotype_type  = phenotype) %>%
  dplyr::filter(!is.na(celltype_flames)) %>%
  dplyr::distinct(celltype_flames, phenotype_type) %>%
  dplyr::mutate(wang_sig = TRUE)
ts("Wang MAGMA loaded.")

# ================================================================
# COLOUR PALETTES (imaging_gen_v2)
# ================================================================
family_colours <- c(
  "Brain size"                  = "#00C28A",
  "Cortical thickness"          = "#0072B2",
  "Curvature"                   = "#E69F00",
  "Microstructural coherence"   = "#CC79A7",
  "Microstructural diffusivity" = "#56B4E9",
  "Orientation dispersion"      = "#D55E00"
)

measure_colors <- c(
  CT="#0072B2", FA="#C0006E", FI="#7DDCB8", GC="#6B4800",
  GMV="#003D28", IC="#2DC48E", ICVF="#E8A8CB", ISOVF="#A8D8F0",
  MC="#F5C842", MD="#0D3F5F", OD="#D55E00", SA="#007A52",
  SubCorVol="#9370DB"
)

# ================================================================
# GLOBAL TYPOGRAPHY (imaging_gen_v2 manuscript-wide convention)
# ================================================================
# Matches the universal typography used elsewhere in the paper
# (Figures 1-3). All theme_bw()/theme_classic() base sizes and
# axis/title text sizes below reference these constants so the
# whole figure set stays visually consistent.
BASE_SIZE  <- 10
TITLE_SIZE <- 10
LABEL_SIZE <- 9
TILE_SIZE  <- 2.8

# Timepoint colours and sizes for panel B concentric circles
group_order <- c("First_trimester","Second_trimester",
                 "Third_trimester","Infancy","Adolescence")
tp_sizes <- c(
  "First_trimester"=2.5, "Second_trimester"=4.5,
  "Third_trimester"=6.5, "Infancy"=8.5, "Adolescence"=10.5)
tp_colours <- c(
  "First_trimester"="#2166AC", "Second_trimester"="#4DAC26",
  "Third_trimester"="#D7191C", "Infancy"="#F4A582",
  "Adolescence"="#762A83")

# ================================================================
# FAMILY DEFINITIONS
# ================================================================
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
  "Microstructural coherence","Orientation dispersion"
)
pheno_to_family <- stack(family_map) %>%
  dplyr::rename(phenotype_type=values, family=ind) %>%
  dplyr::mutate(family=as.character(family), family_label=family)
pheno_family_vec <- setNames(
  rep(names(family_map), sapply(family_map, length)),
  unlist(family_map))


# ================================================================
# LOAD AND PREPARE FLAMES DATA
# ================================================================
ts("Loading FLAMES data...")

# Non-brain SubCorVol phenotypes excluded throughout the paper --
# identical list to Figure 2's exclude_patterns (choroid plexus,
# CSF, optic chiasm, the redundant SubCortGrayVol total, vessel,
# WM hypointensities). Figure 2 drops these before its clustering
# pipeline even runs; applying the same filter here, inside
# prep_flames(), keeps df_full consistent with that -- these were
# never legitimate brain structures, so they shouldn't contribute
# genes to Panel A (or any other panel) either.
nonbrain_subcor_patterns <- c(
  "SubCorVol_choroid.plexus", "SubCorVol_CSF", "SubCorVol_Optic.Chiasm",
  "SubCorVol_SubCortGrayVol", "SubCorVol_vessel", "SubCorVol_WM.hypointensities"
)

prep_flames <- function(df_raw) {
  df_raw %>%
    dplyr::mutate(
      phenotype_clean = phenotype %>%
        str_remove("^global_") %>% str_remove("^regional_"),
      phenotype_type = dplyr::case_when(
        str_detect(phenotype_clean,"^tract_") ~
          str_extract(phenotype_clean,"(?<=^tract_)[A-Za-z0-9]+"),
        TRUE ~ str_split(phenotype_clean,"_") %>% purrr::map_chr(1))) %>%
    dplyr::filter(!grepl(paste(nonbrain_subcor_patterns, collapse="|"),
                         phenotype_clean)) %>%
    dplyr::left_join(pheno_to_family, by="phenotype_type") %>%
    dplyr::filter(!is.na(family))
}

df_raw  <- read.table(flames_path, header=TRUE,
                      stringsAsFactors=FALSE)
df_full <- prep_flames(df_raw)
cat(sprintf("FLAMES rows: %d\n", nrow(df_full)))

# ---- Save FLAMES table ----
# Supplementary Table: one row per gene x family (collapsed).
# Includes all traits (SubCorVol + outliers).
save_flames_table <- function(df_in, suffix) {
  out <- df_in %>%
    dplyr::distinct(symbol, phenotype_type, family) %>%
    dplyr::group_by(symbol, family) %>%
    dplyr::summarise(
      n_phenotypes_in_family = n_distinct(phenotype_type),
      phenotypes_in_family   = paste(sort(unique(phenotype_type)),
                                     collapse="; "),
      .groups="drop") %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(
      n_families   = n_distinct(family),
      all_families = paste(sort(unique(family)), collapse="; ")) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(family, symbol) %>%
    dplyr::select(
      Gene            = symbol,
      Family          = family,
      N_families      = n_families,
      All_families    = all_families,
      N_phenotypes    = n_phenotypes_in_family,
      Phenotypes      = phenotypes_in_family)
  write.csv(out,
            file.path(table_dir,
                      sprintf("SupplTable_FLAMES_genes_%s.csv", suffix)),
            row.names=FALSE)
  cat(sprintf("  FLAMES suppl table saved (%s): %d rows\n",
              suffix, nrow(out)))
}
save_flames_table(df_full, "full")

# ================================================================
# PANEL A: Bar chart + chord diagram
# ================================================================
ts("Building Panel A...")

# Panel-A-only family scheme: SubCorVol stays part of "Brain size",
# but its corpus callosum and ventricle members are pulled out into
# their own entities, using the same colours and the same regex
# classification (^SubCorVol_CC_ ; SubCorVol_...Ventricle|Vent$) as
# Figure 2's get_subcor_group(). Scoped to make_panel_A() only --
# family_colours / family_order / pheno_family_vec stay untouched,
# so Panels B/C, the dendrograms, and the violin/bar consistency
# panels are unaffected.
panelA_family_colours <- c(
  family_colours,
  "Corpus callosum" = "#8B0000",
  "Ventricles"       = "#4169E1")

panelA_family_order <- c(
  "Brain size", "Corpus callosum", "Ventricles",
  family_order[family_order != "Brain size"])

refine_subcor_family <- function(phenotype_clean, phenotype_type, family) {
  dplyr::case_when(
    phenotype_type == "SubCorVol" &
      grepl("^SubCorVol_CC_", phenotype_clean) ~ "Corpus callosum",
    phenotype_type == "SubCorVol" &
      grepl("^SubCorVol_.*(Ventricle|Vent)$", phenotype_clean,
            ignore.case = TRUE) ~ "Ventricles",
    TRUE ~ family
  )
}

make_panel_A <- function(df_in, suffix) {
  # Split corpus callosum / ventricles out of Brain size for this
  # panel only (df_in is a local copy; df_full and pheno_family_vec
  # used elsewhere in the script are untouched).
  df_in <- df_in %>%
    dplyr::mutate(family = refine_subcor_family(
      phenotype_clean, phenotype_type, family))
  
  # Classify genes as shared (>=2 phenotypes in family) or unique
  gene_within <- df_in %>%
    dplyr::distinct(symbol, phenotype_type, family) %>%
    dplyr::group_by(symbol, family) %>%
    dplyr::summarise(
      n_pheno=n_distinct(phenotype_type),
      phenotypes=list(sort(unique(phenotype_type))),
      .groups="drop") %>%
    dplyr::mutate(
      sharing    = dplyr::if_else(n_pheno>=2,"Shared","Unique"),
      fill_label = dplyr::if_else(
        sharing=="Shared",
        paste0(family,"_Shared"),
        paste0(family,"_",purrr::map_chr(phenotypes,1))),
      bar_label  = dplyr::if_else(
        sharing=="Shared","Shared",purrr::map_chr(phenotypes,1)))
  
  phenos_per_fam <- df_in %>%
    dplyr::distinct(family, phenotype_type) %>%
    dplyr::group_by(family) %>%
    dplyr::summarise(phenotypes=list(sort(unique(phenotype_type))),
                     .groups="drop")
  
  all_fill_cols <- unlist(lapply(panelA_family_order, function(fam) {
    row <- which(phenos_per_fam$family==fam)
    if (length(row)==0) return(NULL)
    phenos <- phenos_per_fam$phenotypes[[row]]
    c(setNames(panelA_family_colours[fam], paste0(fam,"_Shared")),
      setNames(sapply(phenos, function(p)
        if (p %in% names(measure_colors)) measure_colors[p]
        else panelA_family_colours[fam]),
        paste0(fam,"_",phenos)))
  }))
  
  bar_df <- gene_within %>%
    dplyr::count(family, fill_label, sharing, bar_label) %>%
    dplyr::mutate(
      family     = factor(family, levels=rev(panelA_family_order)),
      fill_label = factor(fill_label, levels=names(all_fill_cols)))
  
  p_bar <- ggplot(bar_df, aes(x=n, y=family, fill=fill_label)) +
    geom_col(width=0.7, alpha=0.92, colour="white",
             linewidth=0.3) +
    # Only label segments large enough to hold readable text.
    # (Using ifelse() to blank the label -- rather than filtering
    # bar_df down to a subset -- keeps every row in this geom layer
    # so position_stack() still computes each segment's vjust
    # position from the full, correct stack; filtering rows here
    # would silently shift the remaining labels to the wrong
    # segment.)
    geom_text(aes(label=ifelse(n>=8, as.character(n), "")),
              position=position_stack(vjust=0.5),
              size=2.2, colour="black", fontface="bold") +
    geom_text(data=bar_df %>% dplyr::filter(n>=15),
              aes(label=bar_label),
              position=position_stack(vjust=0.5),
              size=1.8, colour="black", vjust=1.5,
              fontface="italic") +
    scale_fill_manual(values=all_fill_cols, guide="none") +
    scale_x_continuous(expand=expansion(mult=c(0,0.05))) +
    # Two-line family labels (str_wrap) instead of one long line --
    # frees up the horizontal room the label column was eating so
    # the bars themselves get more width.
    scale_y_discrete(labels=function(x) str_wrap(x, width=10)) +
    labs(x="Number of FLAMES genes", y=NULL, title="A") +
    theme_bw(base_size=BASE_SIZE) +
    theme(axis.text.y=element_text(size=LABEL_SIZE-1, lineheight=0.85),
          axis.text.x=element_text(size=LABEL_SIZE-1),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=TITLE_SIZE),
          plot.margin=margin(2,2,2,2))
  
  # Chord diagram
  genes_by_fam <- df_in %>%
    dplyr::distinct(symbol, family) %>%
    dplyr::group_by(family) %>%
    dplyr::summarise(genes=list(sort(unique(symbol))),
                     .groups="drop")
  fam_names <- panelA_family_order[panelA_family_order %in% genes_by_fam$family]
  n_fam     <- length(fam_names)
  shared_mat <- matrix(0, n_fam, n_fam,
                       dimnames=list(fam_names,fam_names))
  for (i in seq_len(n_fam)) for (j in seq_len(n_fam)) {
    gi <- genes_by_fam$genes[[match(fam_names[i],
                                    genes_by_fam$family)]]
    gj <- genes_by_fam$genes[[match(fam_names[j],
                                    genes_by_fam$family)]]
    shared_mat[i,j] <- length(intersect(gi,gj))
  }
  shared_plot <- shared_mat; diag(shared_plot) <- 0
  node_cols   <- panelA_family_colours[fam_names]
  
  chord_fn <- function() {
    circos.clear()
    circos.par(gap.after=rep(10,n_fam),
               track.margin=c(0.01,0.01),
               start.degree=90, clock.wise=TRUE)
    chordDiagram(
      shared_plot, order=fam_names, grid.col=node_cols,
      col=adjustcolor(node_cols[row(shared_plot)],alpha.f=0.4),
      transparency=0, annotationTrack="grid",
      preAllocateTracks=list(track.height=mm_h(5)),
      link.sort=TRUE, link.decreasing=TRUE)
    circos.trackPlotRegion(
      track.index=1,
      panel.fun=function(x,y) {
        xlim   <- get.cell.meta.data("xlim")
        ylim   <- get.cell.meta.data("ylim")
        sector <- get.cell.meta.data("sector.index")
        words  <- strsplit(sector," ")[[1]]
        if (length(words)>=2) {
          mid <- ceiling(length(words)/2)
          circos.text(mean(xlim),ylim[1]+(ylim[2]-ylim[1])*0.65,
                      paste(words[seq_len(mid)],collapse=" "),
                      facing="bending.outside",niceFacing=TRUE,
                      cex=0.65,font=2)
          circos.text(mean(xlim),ylim[1]+(ylim[2]-ylim[1])*0.2,
                      paste(words[seq(mid+1,length(words))],
                            collapse=" "),
                      facing="bending.outside",niceFacing=TRUE,
                      cex=0.65,font=2)
        } else {
          circos.text(mean(xlim),mean(ylim),sector,
                      facing="bending.outside",niceFacing=TRUE,
                      cex=0.65,font=2)
        }
      }, bg.border=NA)
    circos.clear()
  }
  
  chord_grob <- cowplot::as_grob(chord_fn)
  panel_A <- plot_grid(p_bar, chord_grob,
                       ncol=2, rel_widths=c(1.3, 1))
  ggsave(file.path(out_dir, sprintf("panelA_%s.pdf",suffix)),
         panel_A, width=11, height=6)
  list(p_bar=p_bar, chord_grob=chord_grob,
       panel_A=panel_A, shared_mat=shared_mat)
}

pa_full <- make_panel_A(df_full, "full")
ts("Panel A done.")

# ================================================================
# CELL TYPE ENRICHMENT
# ================================================================
# See model description at top of script.
# Pseudobulk = Wang 2025, protein-coding genes.
# Enrichment results cached to avoid re-running regressions.

ts("Loading pseudobulk and gene lengths...")
exclude_phenos <- c("SubCorVol","SubCor")

pseudobulk_pc <- readRDS(
  file.path(cache_dir,"pseudobulk_protein_coding.rds"))
pc_pool <- rownames(pseudobulk_pc)
cat(sprintf("Pseudobulk: %d genes x %d bins\n",
            nrow(pseudobulk_pc), ncol(pseudobulk_pc)))

gene_length_cache <- file.path(cache_dir,"gene_lengths_hg19.rds")
if (file.exists(gene_length_cache)) {
  gene_lengths_df <- readRDS(gene_length_cache)
  ts("Gene lengths loaded from cache.")
} else {
  ts("Querying Ensembl GRCh37...")
  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                           dataset="hsapiens_gene_ensembl",
                           host="https://grch37.ensembl.org")
  gene_lengths_raw <- biomaRt::getBM(
    attributes=c("hgnc_symbol","transcript_length"),
    filters="hgnc_symbol", values=pc_pool, mart=mart)
  gene_lengths_df <- gene_lengths_raw %>%
    dplyr::filter(hgnc_symbol!="") %>%
    dplyr::group_by(hgnc_symbol) %>%
    dplyr::summarise(gene_length=max(transcript_length,na.rm=TRUE),
                     .groups="drop") %>%
    dplyr::rename(symbol=hgnc_symbol)
  saveRDS(gene_lengths_df, gene_length_cache)
}

mean_expr_all     <- rowMeans(pseudobulk_pc, na.rm=TRUE)
log_mean_expr_all <- log10(mean_expr_all+0.01)
gene_covars <- data.frame(
  symbol=rownames(pseudobulk_pc),
  log_mean_expr_all=log_mean_expr_all,
  stringsAsFactors=FALSE) %>%
  dplyr::left_join(gene_lengths_df, by="symbol") %>%
  dplyr::mutate(log_length=log10(gene_length+1))
use_length <- !all(is.na(gene_covars$log_length))

run_enrichment <- function(df_in, suffix) {
  cache_file <- file.path(cache_dir,
                          sprintf("enrich_%s.rds",suffix))
  if (file.exists(cache_file)) {
    ts(sprintf("  Loading cached enrichment: %s",suffix))
    return(readRDS(cache_file))
  }
  ts(sprintf("  Running enrichment regressions: %s",suffix))
  df_filt <- df_in %>%
    dplyr::filter(!phenotype_type %in% exclude_phenos)
  genes_by_pheno <- df_filt %>%
    dplyr::distinct(symbol, phenotype_type) %>%
    dplyr::group_by(phenotype_type) %>%
    dplyr::summarise(genes=list(sort(unique(symbol))),
                     n_genes=n_distinct(symbol), .groups="drop")
  pheno_sets <- genes_by_pheno %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pc_genes=list(intersect(genes,pc_pool)),
                  n_pc=length(pc_genes)) %>%
    dplyr::filter(n_pc>=5) %>% dplyr::ungroup()
  
  results <- lapply(seq_len(nrow(pheno_sets)), function(i) {
    ptype  <- pheno_sets$phenotype_type[i]
    pgenes <- pheno_sets$pc_genes[[i]]
    ts(sprintf("    %s (n=%d)",ptype,length(pgenes)))
    lapply(colnames(pseudobulk_pc), function(bn) {
      expr_bin <- pseudobulk_pc[,bn,drop=TRUE]
      mdf <- gene_covars %>%
        dplyr::mutate(
          is_phenotype = as.integer(symbol %in% pgenes),
          log_rel_expr = log10(expr_bin[symbol]+0.01) -
            log_mean_expr_all)
      if (use_length) {
        mdf        <- mdf %>% dplyr::filter(complete.cases(log_length))
        formula_lm <- is_phenotype ~ log_rel_expr + log_length
      } else {
        formula_lm <- is_phenotype ~ log_rel_expr
      }
      fit <- tryCatch(glm(formula_lm,data=mdf,
                          family=binomial()), error=function(e) NULL)
      if (is.null(fit)) return(NULL)
      cf <- tryCatch(summary(fit)$coefficients,
                     error=function(e) NULL)
      if (is.null(cf)||!"log_rel_expr" %in% rownames(cf))
        return(NULL)
      beta_est <- cf["log_rel_expr","Estimate"]
      se_est   <- cf["log_rel_expr","Std. Error"]
      # One-sided test: H1: beta > 0 (positive enrichment only).
      # Using the one-sided p directly (rather than a two-sided
      # Wald p filtered post-hoc by sign) means FDR correction is
      # spent entirely on the positive-enrichment hypothesis, and
      # negative-beta bins are never recorded as candidate hits.
      z_stat <- beta_est / se_est
      data.frame(
        phenotype_type=ptype, bin=bn,
        celltype=gsub("__.*","",bn),
        timepoint=factor(gsub(".*__","",bn),levels=group_order),
        beta=beta_est,
        se=se_est,
        p_value=pnorm(z_stat, lower.tail=FALSE),
        stringsAsFactors=FALSE)
    }) %>% dplyr::bind_rows() %>%
      dplyr::mutate(p_fdr=p.adjust(p_value,method="BH"),
                    significant=!is.na(p_fdr)&p_fdr<0.05&beta>0)
  }) %>% dplyr::bind_rows()
  
  saveRDS(results, cache_file)
  ts(sprintf("  Done: %d significant",
             sum(results$significant,na.rm=TRUE)))
  results
}

enrich_full <- run_enrichment(df_full, "full")

write.csv(enrich_full %>%
            dplyr::mutate(across(where(is.numeric),~round(.x,4))),
          file.path(table_dir,"celltype_enrichment_full.csv"),
          row.names=FALSE)

# ---- Supplementary Table: Cell type enrichment (full dataset) ----
ct_suppl_table <- enrich_full %>%
  dplyr::filter(significant) %>%
  dplyr::select(
    Phenotype   = phenotype_type,
    Cell_type   = celltype,
    Timepoint   = timepoint,
    Beta        = beta,
    SE          = se,
    P_value     = p_value,
    FDR         = p_fdr) %>%
  dplyr::mutate(
    Family    = pheno_family_vec[Phenotype],
    across(where(is.numeric), ~round(.x, 4))) %>%
  dplyr::select(Family, everything()) %>%
  dplyr::arrange(Family, Phenotype, Cell_type, Timepoint)
write.csv(ct_suppl_table,
          file.path(table_dir,
                    "SupplTable_CellType_enrichment_full.csv"),
          row.names=FALSE)
cat(sprintf("  CT suppl table saved: %d rows\n", nrow(ct_suppl_table)))
ts("Enrichment tables saved.")

# ---- Panel B: concentric circles ----
# Each cell type x phenotype position shows one hollow circle
# per significant timepoint. Size and colour encode timepoint.
# Largest circles drawn first (background); smallest on top.
# No family labels; all axis text in black.
make_panel_B <- function(enrich_results, suffix) {
  sig <- enrich_results %>%
    dplyr::filter(significant) %>%
    dplyr::mutate(
      family    = pheno_family_vec[phenotype_type],
      family    = factor(family, levels=family_order),
      timepoint = factor(timepoint, levels=group_order))
  if (nrow(sig)==0) {
    cat(sprintf("No significant bins: %s\n",suffix))
    return(NULL)
  }
  
  # Biologically ordered cell types
  ct_hierarchy <- c(
    "RG-vRG","RG-oRG","RG-tRG","IPC-EN","Tri-IPC",
    "EN-IT-Immature","EN-L2_3-IT","EN-L4-IT","EN-IT",
    "IN-NCx_dGE-Immature","IN-MGE-Immature","IN-MGE-SST",
    "IN-MGE-PV","IN-CGE-SNCG","IN-CGE-VIP",
    "Cajal-Retzius cell",
    "Astrocyte-Immature","Astrocyte-Protoplasmic",
    "Astrocyte-Fibrous","Oligodendrocyte","OPC")
  present_cts    <- unique(sig$celltype)
  celltype_order <- c(ct_hierarchy[ct_hierarchy %in% present_cts],
                      sort(setdiff(present_cts, ct_hierarchy)))
  
  x_order <- data.frame(
    phenotype_type=names(pheno_family_vec),
    family=pheno_family_vec, stringsAsFactors=FALSE) %>%
    dplyr::filter(phenotype_type %in% unique(sig$phenotype_type)) %>%
    dplyr::mutate(family=factor(family,levels=family_order)) %>%
    dplyr::arrange(family, phenotype_type) %>%
    dplyr::pull(phenotype_type)
  
  fam_bounds <- data.frame(
    phenotype_type=names(pheno_family_vec),
    family=pheno_family_vec, stringsAsFactors=FALSE) %>%
    dplyr::filter(phenotype_type %in% x_order) %>%
    dplyr::mutate(family=factor(family,levels=family_order),
                  x_pos=match(phenotype_type,x_order)) %>%
    dplyr::group_by(family) %>%
    dplyr::summarise(x_end=max(x_pos)+0.5,
                     x_mid=mean(x_pos),.groups="drop")
  fam_sep_x <- fam_bounds$x_end[seq_len(nrow(fam_bounds)-1)]
  n_ct <- length(celltype_order)
  n_x  <- length(x_order)
  
  plot_df <- sig %>%
    dplyr::mutate(
      celltype       = factor(celltype, levels=celltype_order),
      phenotype_type = factor(phenotype_type, levels=x_order),
      tp_size        = tp_sizes[as.character(timepoint)],
      tp_colour      = tp_colours[as.character(timepoint)]) %>%
    dplyr::arrange(desc(tp_size))
  
  # Wang MAGMA inner circle: one dot per cell type x phenotype
  # (no timepoint dimension - sits at centre of concentric stack)
  # Size = smallest timepoint circle / 2, fixed dark grey fill.
  wang_inner <- wang_sig_lookup %>%
    dplyr::filter(
      celltype_flames %in% celltype_order,
      phenotype_type  %in% x_order) %>%
    dplyr::mutate(
      celltype       = factor(celltype_flames, levels=celltype_order),
      phenotype_type = factor(phenotype_type,  levels=x_order))
  
  p <- ggplot(plot_df, aes(x=phenotype_type, y=celltype)) +
    geom_vline(xintercept=fam_sep_x,
               colour="grey60", linewidth=0.4, linetype="dashed") +
    # Concentric hollow circles - largest (oldest) drawn first
    geom_point(aes(size=I(tp_size), colour=I(tp_colour)),
               shape=21, fill="white", stroke=1.2, alpha=0.95) +
    # Wang MAGMA inner filled circle (dark grey, fixed small size)
    { if (nrow(wang_inner) > 0)
      geom_point(data=wang_inner,
                 aes(x=phenotype_type, y=celltype),
                 size=2.5, shape=21,
                 fill="#444444", colour="#444444",
                 stroke=0, inherit.aes=FALSE)
      else list() } +
    scale_x_discrete(expand=expansion(add=0.3)) +
    scale_y_discrete(expand=expansion(add=0.3)) +
    labs(x="Phenotype (grouped by family)", y="Cell type", title="B") +
    theme_bw(base_size=BASE_SIZE) +
    theme(
      axis.text.x      = element_text(angle=45, hjust=1,
                                      size=LABEL_SIZE, colour="black"),
      axis.text.y      = element_text(size=LABEL_SIZE, colour="black"),
      axis.title       = element_text(size=TITLE_SIZE),
      plot.title       = element_text(face="bold", size=TITLE_SIZE),
      panel.grid.major = element_line(colour="grey94"),
      plot.margin      = margin(t=10, r=15, b=5, l=5),
      legend.position  = "none")
  
  # Timepoint legend
  leg_df <- data.frame(
    size=tp_sizes[group_order], col=tp_colours[group_order],
    label=c("First trimester","Second trimester",
            "Third trimester","Infancy","Adolescence"),
    x=seq_along(group_order), y=1, stringsAsFactors=FALSE)
  p_leg <- ggplot(leg_df, aes(x=x,y=y)) +
    geom_point(aes(size=I(size*1.5), colour=I(col)),
               shape=21, fill="white", stroke=1.5) +
    geom_text(aes(label=label), vjust=-1.8, size=3.5,
              colour="grey20") +
    # Wang MAGMA legend entry
    annotate("point", x=length(group_order)+1.2, y=1,
             size=4.0, shape=21, fill="#444444", colour="#444444") +
    annotate("text",  x=length(group_order)+1.2, y=1,
             label="Wang MAGMA\n(FDR<0.05)", vjust=-1.8,
             size=3.5, colour="grey20") +
    scale_x_continuous(limits=c(0.3, length(group_order)+2.2)) +
    scale_y_continuous(limits=c(0.5,1.9)) +
    theme_void()
  
  p_with_leg <- cowplot::plot_grid(
    p, p_leg, ncol=1, rel_heights=c(1,0.06))
  
  h <- max(5, n_ct*0.28+2.5)
  w <- max(6, n_x*0.55+3)
  ggsave(file.path(out_dir, sprintf("panelB_%s.pdf",suffix)),
         p_with_leg, width=w, height=h, device=pdf,
         limitsize=FALSE)
  list(plot=p, legend=p_leg, n_y=n_ct, n_x=n_x)
}

pb_full  <- make_panel_B(enrich_full,  "full")
ts("Panel B done.")

# ================================================================
# GO ENRICHMENT
# ================================================================
# ORA via clusterProfiler::enrichGO.
# pvalueCutoff=1 keeps ALL terms for similarity matrices.
# simplify() run only on significant terms (fast).
# Two caches per suffix:
#   go_simplified_SUFFIX.rds - significant simplified terms
#                              (used for panel C)
#   go_all_SUFFIX.rds        - all tested terms
#                              (used for similarity matrices)

ts("Running GO enrichment...")
background_pc   <- rownames(pseudobulk_pc)
map_entrez <- function(symbols) {
  if (length(symbols)==0) return(character(0))
  suppressMessages(bitr(symbols, fromType="SYMBOL",
                        toType="ENTREZID",
                        OrgDb=org.Hs.eg.db))$ENTREZID
}
entrez_universe <- map_entrez(background_pc)

run_go <- function(df_in, suffix) {
  cache_sig <- file.path(cache_dir,
                         sprintf("go_simplified_%s.rds",suffix))
  cache_all <- file.path(cache_dir,
                         sprintf("go_all_%s.rds",suffix))
  
  if (file.exists(cache_sig) && file.exists(cache_all)) {
    ts(sprintf("  Loading cached GO: %s",suffix))
    return(list(simplified=readRDS(cache_sig),
                all=readRDS(cache_all)))
  }
  ts(sprintf("  Running GO ORA: %s",suffix))
  
  genes_by_pheno <- df_in %>%
    dplyr::distinct(symbol, phenotype_type, family) %>%
    dplyr::group_by(phenotype_type, family) %>%
    dplyr::summarise(genes=list(sort(unique(symbol))),
                     n_genes=n_distinct(symbol),
                     .groups="drop") %>%
    dplyr::mutate(family=factor(family,levels=family_order)) %>%
    dplyr::arrange(family, phenotype_type) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(entrez_genes=list(map_entrez(genes)),
                  n_entrez=length(entrez_genes)) %>%
    dplyr::ungroup()
  
  ora_results <- lapply(seq_len(nrow(genes_by_pheno)), function(i) {
    ptype <- genes_by_pheno$phenotype_type[i]
    eg    <- genes_by_pheno$entrez_genes[[i]]
    fam   <- as.character(genes_by_pheno$family[i])
    if (length(eg)<5) return(NULL)
    ts(sprintf("    ORA: %s (n=%d)",ptype,length(eg)))
    res <- tryCatch(
      enrichGO(gene=eg, universe=entrez_universe,
               OrgDb=org.Hs.eg.db, ont="BP",
               pAdjustMethod="BH",
               pvalueCutoff=1, qvalueCutoff=1,
               minGSSize=10, maxGSSize=500,
               readable=TRUE),
      error=function(e) NULL)
    if (is.null(res)||nrow(res@result)==0) return(NULL)
    
    # simplify() only on significant terms - much faster
    res_sig        <- res
    res_sig@result <- res_sig@result %>%
      dplyr::filter(p.adjust<0.05)
    simplified_result <- if (nrow(res_sig@result)>0) {
      tryCatch(
        clusterProfiler::simplify(res_sig, cutoff=0.7,
                                  by="p.adjust",
                                  select_fun=min)@result,
        error=function(e) res_sig@result)
    } else data.frame()
    
    annotate_df <- function(df) df %>%
      dplyr::mutate(
        phenotype_type = ptype, family=fam,
        significant    = p.adjust<0.05,
        GeneRatio_num  = sapply(strsplit(GeneRatio,"/"),
                                function(v) as.numeric(v[1])/
                                  as.numeric(v[2])),
        neg_log10_fdr  = -log10(pmax(p.adjust,1e-300)))
    
    list(
      simplified = if (nrow(simplified_result)>0)
        annotate_df(simplified_result) else NULL,
      all        = annotate_df(res@result))
  })
  ora_results <- Filter(Negate(is.null), ora_results)
  
  go_simplified <- dplyr::bind_rows(
    lapply(ora_results,function(x) x$simplified)) %>%
    dplyr::filter(!is.na(ID))
  go_all <- dplyr::bind_rows(
    lapply(ora_results,function(x) x$all)) %>%
    dplyr::filter(!is.na(ID))
  
  saveRDS(go_simplified, cache_sig)
  saveRDS(go_all,        cache_all)
  ts(sprintf("  Done: %d simplified, %d total, %d significant",
             nrow(go_simplified), nrow(go_all),
             sum(go_all$significant,na.rm=TRUE)))
  list(simplified=go_simplified, all=go_all)
}

go_res_full <- run_go(df_full, "full")

# Unpack: go_full for panel C; go_full_all for similarity matrix
go_full     <- go_res_full$simplified
go_full_all <- go_res_full$all

write.csv(go_full %>%
            dplyr::filter(significant) %>%
            dplyr::mutate(across(where(is.numeric),~round(.x,4))),
          file.path(table_dir,"go_enrichment_significant_full.csv"),
          row.names=FALSE)

# ---- Supplementary Table: GO enrichment (full dataset) ----
go_suppl_table <- go_full %>%
  dplyr::filter(significant) %>%
  dplyr::select(
    Family          = family,
    Phenotype       = phenotype_type,
    GO_ID           = ID,
    Term            = Description,
    GeneRatio,
    BgRatio,
    FDR             = p.adjust,
    neg_log10_FDR   = neg_log10_fdr,
    Genes           = geneID) %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  dplyr::arrange(Family, Phenotype, FDR)
write.csv(go_suppl_table,
          file.path(table_dir, "SupplTable_GO_enrichment_full.csv"),
          row.names=FALSE)
cat(sprintf("  GO suppl table saved: %d rows\n", nrow(go_suppl_table)))
ts("GO tables saved.")

# ---- Panel C: Unified GO dotplot with semantic deduplication ----
# Term selection pipeline:
#   1. Build GOSemSim godata object (cached).
#   2. Discover shared terms: significant in >=4 families,
#      cluster by Wang semantic similarity (threshold=0.5),
#      pick top-FDR representative per cluster,
#      select 3 most broadly significant cluster reps.
#   3. Unique terms per family: from remaining significant terms,
#      cluster by similarity (threshold=0.5, also vs shared terms),
#      pick top-FDR rep per cluster, select 2 most dissimilar.
#   4. Single panel: phenotypes on x, terms on y, dot size =
#      -log10(FDR), fill = family colour (unique) / grey (shared).
#      Dashed vertical lines separate families; dashed horizontal
#      lines separate unique-term family blocks and shared block.

ts("Building Panel C with semantic deduplication...")

# ---- GOSemSim setup (cached) ----
godata_cache <- file.path(cache_dir, "godata_BP_human.rds")
if (file.exists(godata_cache)) {
  ts("  Loading cached GOSemSim godata...")
  go_data_obj <- readRDS(godata_cache)
} else {
  ts("  Building GOSemSim godata (BP, human)...")
  go_data_obj <- GOSemSim::godata("org.Hs.eg.db",
                                  ont="BP", computeIC=FALSE)
  saveRDS(go_data_obj, godata_cache)
}

# ---- Semantic similarity helpers ----
# Returns pairwise Wang similarity matrix for a vector of GO IDs.
# Terms not in the ontology get NA rows/cols.
go_sim_mat <- function(go_ids) {
  n <- length(go_ids)
  if (n < 2) return(matrix(1, 1, 1, dimnames=list(go_ids, go_ids)))
  mat <- matrix(NA_real_, n, n, dimnames=list(go_ids, go_ids))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i == j) { mat[i,j] <- 1; next }
    if (!is.na(mat[j,i])) { mat[i,j] <- mat[j,i]; next }
    mat[i,j] <- tryCatch(
      GOSemSim::goSim(go_ids[i], go_ids[j],
                      semData=go_data_obj, measure="Wang"),
      error=function(e) NA_real_)
  }
  mat
}

# ---- Greedy dissimilar term picker ----
# Given a data frame of candidates (ID, Description, mean_fdr)
# and a set of already-selected IDs to avoid similarity with,
# greedily pick n_pick terms such that each new term has
# Wang similarity < sim_threshold to all already-selected terms
# AND to all previously picked terms.
# Falls back to top mean_fdr if not enough dissimilar terms exist.
pick_dissimilar <- function(candidates, already_selected_ids,
                            n_pick=2, sim_threshold=0.5) {
  if (nrow(candidates) == 0) return(candidates[0,])
  
  # Pre-compute pairwise sim matrix for all candidates + seeds
  all_ids <- unique(c(already_selected_ids, candidates$ID))
  sim_mat <- go_sim_mat(all_ids)
  
  selected_ids <- already_selected_ids
  selected_rows <- candidates[0,]
  
  # Sort candidates by mean_fdr descending - take strongest first
  cands_sorted <- candidates %>% dplyr::arrange(desc(mean_fdr))
  
  for (i in seq_len(nrow(cands_sorted))) {
    cid <- cands_sorted$ID[i]
    # Similarity of this candidate to everything already selected
    sims_to_selected <- sapply(selected_ids, function(sid) {
      v <- sim_mat[cid, sid]
      if (is.na(v)) 0 else v
    })
    if (length(sims_to_selected) == 0 ||
        max(sims_to_selected) < sim_threshold) {
      selected_rows <- dplyr::bind_rows(selected_rows,
                                        cands_sorted[i,])
      selected_ids  <- c(selected_ids, cid)
    }
    if (nrow(selected_rows) >= n_pick) break
  }
  
  # Fallback: fill remaining slots with top-FDR terms not yet picked
  if (nrow(selected_rows) < n_pick) {
    remaining <- cands_sorted %>%
      dplyr::filter(!ID %in% selected_rows$ID)
    n_need <- n_pick - nrow(selected_rows)
    selected_rows <- dplyr::bind_rows(
      selected_rows,
      remaining[seq_len(min(n_need, nrow(remaining))),])
  }
  selected_rows
}

# ---- Build term summary per family ----
all_fams <- family_order

term_summary <- go_full %>%
  dplyr::filter(significant) %>%
  dplyr::group_by(family, ID, Description) %>%
  dplyr::summarise(mean_fdr = mean(neg_log10_fdr, na.rm=TRUE),
                   n_phenos = n_distinct(phenotype_type),
                   .groups="drop")

# ---- Step 1: Shared terms ----
# From SupplTable_GO_semantic_clusters.csv: terms significant in
# >=4 families, deduplicated by pairwise Wang similarity.
# Full >=4-family terms sorted by mean_fdr:
#   GO:0016055  Wnt signaling pathway          4 fams  FDR 3.46
#   GO:0198738  cell-cell signaling by wnt     5 fams  FDR 3.09  <- sim 0.623 to GO:0016055
#   GO:0048589  developmental growth           5 fams  FDR 2.03
#   GO:0001667  ameboidal-type cell migration  5 fams  FDR 1.93
#   GO:0048588  developmental cell growth      5 fams  FDR 1.57
#
# GO:0016055 and GO:0198738 are both Wnt terms with Wang sim 0.623
# - above threshold, so we keep only GO:0016055 (higher FDR, canonical).
# That leaves 4 candidates; pick top 3 by mean_fdr that are mutually
# dissimilar (all pairs < 0.5):

shared_candidates <- term_summary %>%
  dplyr::group_by(ID, Description) %>%
  dplyr::summarise(n_fams   = n_distinct(family),
                   mean_fdr = mean(mean_fdr, na.rm=TRUE),
                   .groups="drop") %>%
  dplyr::filter(n_fams >= 4) %>%
  dplyr::arrange(desc(mean_fdr))

shared_terms_df <- pick_dissimilar(shared_candidates,
                                   already_selected_ids=character(0),
                                   n_pick=3, sim_threshold=0.5)

shared_ids   <- shared_terms_df$ID
shared_descs <- shared_terms_df$Description

cat(sprintf("  Shared terms (%d):\n", length(shared_descs)))
for (i in seq_along(shared_descs))
  cat(sprintf("    %d. %s (n_fams=%d, mean_fdr=%.2f)\n",
              i, shared_descs[i],
              shared_terms_df$n_fams[i],
              shared_terms_df$mean_fdr[i]))

# ---- Step 2: Unique terms per family ----
# For each family: take significant terms not already shared,
# pick 2 that are mutually dissimilar AND dissimilar to shared terms.

unique_terms_by_fam <- lapply(all_fams, function(fam) {
  candidates <- term_summary %>%
    dplyr::filter(family == fam, !ID %in% shared_ids) %>%
    dplyr::arrange(desc(mean_fdr))
  pick_dissimilar(candidates,
                  already_selected_ids=shared_ids,
                  n_pick=2, sim_threshold=0.5)
})
names(unique_terms_by_fam) <- all_fams

cat("  Unique terms per family:\n")
for (f in all_fams) {
  descs <- unique_terms_by_fam[[f]]$Description
  cat(sprintf("    %s: %s\n", f, paste(descs, collapse=" | ")))
}

# ---- Step 3: Build plot data ----
shared_plot_data <- go_full %>%
  dplyr::filter(significant, ID %in% shared_ids)

unique_plot_data <- dplyr::bind_rows(lapply(all_fams, function(fam) {
  ids <- unique_terms_by_fam[[fam]]$ID
  go_full %>% dplyr::filter(significant, family==fam, ID %in% ids)
}))

plot_data <- dplyr::bind_rows(shared_plot_data, unique_plot_data)
if (nrow(plot_data) == 0) {
  cat("No data for Panel C\n")
  pc_full <- NULL
} else {
  
  # ---- Step 4: Y-axis order ----
  # Bottom -> top: unique terms (family order, low FDR first within
  # family) -> dashed separator -> 3 shared terms (weakest at bottom)
  unique_y <- unlist(lapply(all_fams, function(fam) {
    df <- unique_terms_by_fam[[fam]]
    if (nrow(df) == 0) return(character(0))
    df %>% dplyr::arrange(mean_fdr) %>% dplyr::pull(Description)
  }))
  shared_y <- shared_terms_df %>%
    dplyr::arrange(mean_fdr) %>% dplyr::pull(Description)
  
  y_order_raw <- unique(c(unique_y, shared_y))
  
  # ---- Step 5: Separator positions ----
  # Separators sit between family unique blocks and between the
  # unique block and the shared block.
  # geom_hline on a discrete y-axis uses the numeric index of the
  # factor level, so we look up the position of each family's
  # last unique term in y_order_raw, then add 0.5.
  
  # Position of each term in y_order_raw (1 = bottom row)
  term_pos <- setNames(seq_along(y_order_raw), y_order_raw)
  
  # For each family with unique terms, find the position of its
  # last (highest) unique term in y_order_raw
  hsep_y <- c()
  for (f in all_fams) {
    fam_unique <- unique_terms_by_fam[[f]]$Description
    if (length(fam_unique) == 0) next
    last_pos <- max(term_pos[fam_unique], na.rm=TRUE)
    hsep_y <- c(hsep_y, last_pos + 0.5)
  }
  # Add separator between unique block and shared block
  # (position of last unique term overall + 0.5)
  n_unique_total <- length(unique_y)
  if (n_unique_total > 0)
    hsep_y <- c(hsep_y, term_pos[unique_y[length(unique_y)]] + 0.5)
  # Deduplicate (in case last family separator == shared separator)
  hsep_y <- sort(unique(hsep_y))
  
  # ---- Step 6: X-axis ----
  x_order <- data.frame(
    phenotype_type = names(pheno_family_vec),
    family         = pheno_family_vec,
    stringsAsFactors=FALSE) %>%
    dplyr::filter(phenotype_type %in% unique(plot_data$phenotype_type)) %>%
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
  
  # ---- Step 7: Dot colours by phenotype measure ----
  # Each dot coloured by its phenotype_type using measure_colors.
  # No colour legend shown - phenotype labels on x-axis suffice.
  
  # ---- Step 8: Assemble plot data with factor levels ----
  wrap_w <- 30
  y_wrapped <- str_wrap(y_order_raw, width=wrap_w)
  
  plot_data <- plot_data %>%
    dplyr::mutate(
      Description    = str_wrap(Description, width=wrap_w),
      Description    = factor(Description, levels=y_wrapped),
      # IMPORTANT: compute dot_colour from the character phenotype_type
      # BEFORE converting it to a factor below. measure_colors[] is a
      # named lookup -- once phenotype_type becomes a factor, indexing
      # switches to positional (factor integer code), silently mapping
      # each phenotype to the wrong colour instead of its own.
      dot_colour     = measure_colors[as.character(phenotype_type)],
      phenotype_type = factor(phenotype_type, levels=x_order))
  
  fdr_range <- range(plot_data$neg_log10_fdr, na.rm=TRUE)
  
  # ---- Step 9: Plot ----
  make_panel_C_plot <- function(go_in_sub, suffix) {
    pd <- plot_data %>%
      dplyr::filter(phenotype_type %in%
                      unique(go_in_sub$phenotype_type[
                        go_in_sub$significant]))
    if (nrow(pd) == 0) return(NULL)
    
    p <- ggplot(pd,
                aes(x=phenotype_type, y=Description,
                    size=neg_log10_fdr, fill=I(dot_colour))) +
      geom_vline(xintercept=fam_sep_x,
                 colour="grey55", linewidth=0.4, linetype="dashed") +
      geom_hline(yintercept=hsep_y,
                 colour="grey55", linewidth=0.4, linetype="dashed") +
      geom_point(shape=21, colour="white", stroke=0.3, alpha=0.92) +
      scale_size_continuous(
        name   = expression(-log[10](FDR)),
        range  = c(1.5, 6),
        breaks = pretty(fdr_range, n=4),
        limits = fdr_range) +
      scale_x_discrete(expand=expansion(add=0.4)) +
      scale_y_discrete(expand=expansion(add=0.4)) +
      labs(x="Phenotype (grouped by family)", y=NULL, title="C") +
      theme_bw(base_size=BASE_SIZE) +
      theme(
        axis.text.x      = element_text(angle=45, hjust=1,
                                        size=LABEL_SIZE, colour="black"),
        axis.text.y      = element_text(size=LABEL_SIZE-0.5, lineheight=0.9,
                                        colour="black"),
        axis.title.x     = element_text(size=TITLE_SIZE),
        plot.title       = element_text(face="bold", size=TITLE_SIZE),
        panel.grid.major = element_line(colour="grey94",
                                        linewidth=0.25),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(linewidth=0.5),
        legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.title     = element_text(size=LABEL_SIZE),
        legend.text      = element_text(size=LABEL_SIZE-1),
        plot.margin      = margin(5,5,5,5)) +
      guides(
        size = guide_legend(title.position="left",
                            override.aes=list(fill="grey40")),
        fill = "none")
    
    n_y <- length(y_order_raw)
    n_x <- length(x_order)
    go_h <- max(5, n_y * 0.38 + 2.5)
    go_w <- max(8, n_x * 0.52 + 4)
    ggsave(file.path(out_dir, sprintf("panelC_%s.pdf", suffix)),
           p, width=go_w, height=go_h, limitsize=FALSE)
    p
  }
  
  pc_full <- make_panel_C_plot(go_full, "full")
  
} # end nrow(plot_data) > 0
ts("Panel C done.")

# ================================================================
# MAIN FIGURE ASSEMBLY
# ================================================================
ts("Assembling Figure 4...")

assemble_fig4 <- function(pa, pb_obj, pc, suffix) {
  if (is.null(pa)||is.null(pb_obj)||is.null(pc)) {
    cat(sprintf("Skipping Figure4_%s — missing panel\n",suffix))
    return(invisible(NULL))
  }
  pb_with_leg <- cowplot::plot_grid(
    pb_obj$plot, pb_obj$legend,
    ncol=1, rel_heights=c(1,0.08))
  # Panel A is back to its original short/wide (bar | chord)
  # shape. Trimmed its height share slightly (vs. Panel C) now
  # that its labels and margins are tighter and it needs less
  # room to breathe than before.
  left_col <- plot_grid(pa$panel_A, pc, ncol=1,
                        rel_heights=c(0.85, 1.6))
  # Right column (Panel B) width: still scales mildly with the
  # number of phenotype columns it has to display, but capped
  # tighter than before since the figure's overall footprint is
  # smaller now (12x10in rather than 26x20in) -- an uncapped B
  # panel was starving the left column of the width it needs for
  # readable bar/chord/GO-term labels.
  b_width <- max(0.34, min(0.46, pb_obj$n_x*0.045))
  fig4 <- plot_grid(left_col, pb_with_leg, nrow=1,
                    rel_widths=c(1-b_width, b_width))
  ggsave(file.path(out_dir,sprintf("Figure4_%s.pdf",suffix)),
         fig4, width=12, height=10, limitsize=FALSE)
  ggsave(file.path(out_dir,sprintf("Figure4_%s.png",suffix)),
         fig4, width=12, height=10, dpi=300, limitsize=FALSE)
  cat(sprintf("Figure4_%s saved.\n",suffix))
}

assemble_fig4(pa_full, pb_full, pc_full, "full")

# ================================================================
# SIMILARITY MATRIX HELPER FUNCTIONS
# ================================================================

# FLAMES Jaccard similarity matrix
make_jac_mat <- function(df_in) {
  phenos <- unique(df_in$phenotype_type)
  mat <- outer(phenos,phenos,FUN=Vectorize(function(a,b) {
    ga <- unique(df_in$symbol[df_in$phenotype_type==a])
    gb <- unique(df_in$symbol[df_in$phenotype_type==b])
    length(intersect(ga,gb))/length(union(ga,gb))
  }))
  dimnames(mat) <- list(phenos,phenos); mat
}

# CT: paired using intersection of converged bins
make_ct_mat <- function(enrich_in) {
  beta_mat <- enrich_in %>%
    dplyr::select(phenotype_type,bin,beta) %>%
    tidyr::pivot_wider(names_from=bin,values_from=beta,
                       values_fn=mean) %>%
    tibble::column_to_rownames("phenotype_type") %>%
    as.matrix()
  beta_mat <- beta_mat[,colSums(is.na(beta_mat))==0,drop=FALSE]
  if (ncol(beta_mat)<2||nrow(beta_mat)<2) return(NULL)
  cor(t(beta_mat),use="complete.obs")
}

# GO: Spearman rank correlation on ALL tested terms, using ranks
# within each phenotype's profile across all tested GO terms.
# Zero-variance rows (all-zero profiles) removed before correlating.
make_go_mat <- function(go_in_all) {
  go_in_all %>%
    dplyr::filter(!is.na(ID),!is.na(p.adjust)) %>%
    dplyr::mutate(nlp=-log10(pmax(p.adjust,1e-300))) %>%
    dplyr::group_by(ID) %>%
    dplyr::filter(n_distinct(phenotype_type)>=2) %>%
    dplyr::ungroup() %>%
    dplyr::select(phenotype_type,ID,nlp) %>%
    tidyr::pivot_wider(names_from=ID,values_from=nlp,
                       values_fill=0,values_fn=mean) %>%
    tibble::column_to_rownames("phenotype_type") %>%
    as.matrix() %>%
    {if (nrow(.)>=2) cor(t(.),use="complete.obs",method="spearman") else NULL}
}

# Compute similarity matrices now - used by Supp1, Supp2, Supp3
ct_mat_full <- make_ct_mat(enrich_full)
go_mat_full <- make_go_mat(go_full_all)

cat("CT similarity matrix — rows:",
    length(rownames(ct_mat_full)), "\n")
cat("GO similarity matrix — rows:",
    length(rownames(go_mat_full)), "\n")

# ================================================================
# SUPPLEMENTARY FIGURE 1: Phenotype clustering dendrograms
# ================================================================
ts("Building Supp1 (hclust)...")

make_dendro_plot <- function(hc_obj, title_str) {
  dd   <- ggdendro::dendro_data(hc_obj)
  labs <- dd$labels
  labs$label <- as.character(labs$label)
  lab_cols <- family_colours[pheno_family_vec[labs$label]]
  lab_cols[is.na(lab_cols)] <- "grey30"
  ggplot() +
    geom_segment(data=dd$segments,
                 aes(x=x,y=y,xend=xend,yend=yend),
                 colour="#666666",linewidth=0.4) +
    geom_text(data=labs,
              aes(x=x,y=y-0.01,label=label),
              angle=90,hjust=1,size=3.5,colour=lab_cols) +
    scale_y_continuous(expand=expansion(mult=c(0.25,0.05))) +
    theme_classic(base_size=BASE_SIZE) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          plot.title=element_text(size=TITLE_SIZE,face="bold"),
          plot.margin=margin(5,5,5,5)) +
    labs(title=title_str,x=NULL,y="Distance")
}

make_supp1 <- function(df_in, enrich_in, go_in_all, suffix) {
  # i. FLAMES Jaccard
  phenos  <- unique(df_in$phenotype_type)
  jac_mat <- outer(phenos,phenos,FUN=Vectorize(function(a,b) {
    ga <- unique(df_in$symbol[df_in$phenotype_type==a])
    gb <- unique(df_in$symbol[df_in$phenotype_type==b])
    length(intersect(ga,gb))/length(union(ga,gb))
  }))
  dimnames(jac_mat) <- list(phenos,phenos)
  hc1 <- hclust(as.dist(1-jac_mat),method="average")
  p1  <- make_dendro_plot(hc1,
                          "i. FLAMES gene overlap\n(Jaccard distance)")
  
  # ii. Cell type enrichment
  ct_cor <- make_ct_mat(enrich_in)
  if (!is.null(ct_cor)&&nrow(ct_cor)>1) {
    hc2 <- hclust(as.dist(1-ct_cor),method="average")
    p2  <- make_dendro_plot(hc2,
                            "ii. Cell type enrichment\n(1 - Pearson r)")
  } else {
    p2 <- ggplot()+annotate("text",x=0.5,y=0.5,
                            label="Insufficient data")+theme_void()
  }
  
  # iii. GO enrichment (all terms, Spearman)
  go_cor <- make_go_mat(go_in_all)
  if (!is.null(go_cor)&&nrow(go_cor)>1) {
    hc3 <- hclust(as.dist(1-go_cor),method="average")
    p3  <- make_dendro_plot(hc3,
                            "iii. GO enrichment\n(1 - Pearson r)")
  } else {
    p3 <- ggplot()+annotate("text",x=0.5,y=0.5,
                            label="Insufficient data")+theme_void()
  }
  
  leg_df <- data.frame(
    family=names(family_colours),fill=unname(family_colours),
    x=seq_along(family_colours),y=1)
  p_leg <- ggplot(leg_df,aes(x=x,y=y,fill=fill,label=family)) +
    geom_tile(width=0.8,height=0.6) +
    geom_text(aes(y=0.55),size=2.8,hjust=0.5) +
    scale_fill_identity() +
    theme_void()+theme(plot.margin=margin(2,2,2,2))
  
  p_row   <- plot_grid(p1,p2,p3,nrow=1,
                       labels=c("i","ii","iii"),label_size=10)
  p_supp1 <- plot_grid(p_row,p_leg,ncol=1,rel_heights=c(1,0.1))
  ggsave(file.path(out_dir,
                   sprintf("Supp1_hclust_%s.pdf",suffix)),
         p_supp1, width=18, height=7)
  cat(sprintf("Supp1_hclust_%s.pdf saved.\n",suffix))
}

make_supp1(df_full, enrich_full, go_full_all, "full")
ts("Supp1 done.")

# ================================================================
# FAMILY SIMILARITY STATISTICS
# ================================================================
# compute_family_similarity(): for a similarity matrix, classifies
# all unique phenotype pairs as within or between family, and tests
# whether within-family similarity exceeds between-family via
# permutation (family labels shuffled n_perm times). A Wilcoxon
# rank-sum test is NOT used here: it assumes independent samples,
# but pairwise similarity values are not independent (each phenotype
# appears in multiple pairs), so its p-value would be anti-conservative.
# The permutation test respects this structure by shuffling the
# independent units (phenotype family labels) rather than the
# derived pairwise values themselves.

compute_family_similarity <- function(sim_mat, label) {
  if (is.null(sim_mat)) return(NULL)
  phenos <- rownames(sim_mat)
  pair_df <- expand.grid(a=phenos,b=phenos,
                         stringsAsFactors=FALSE) %>%
    dplyr::filter(a<b) %>%
    dplyr::mutate(
      fam_a      = pheno_family_vec[a],
      fam_b      = pheno_family_vec[b],
      within     = fam_a==fam_b,
      pair_type  = dplyr::if_else(within,"Within family",
                                  "Between family"),
      similarity = mapply(function(x,y) sim_mat[x,y],a,b),
      version    = label) %>%
    dplyr::filter(!is.na(similarity))
  
  within_s  <- pair_df$similarity[pair_df$within]
  between_s <- pair_df$similarity[!pair_df$within]
  if (length(within_s)<2||length(between_s)<2) return(NULL)
  obs_diff <- mean(within_s,na.rm=TRUE)-mean(between_s,na.rm=TRUE)
  fam_vec    <- pheno_family_vec[phenos]
  null_diffs <- replicate(n_perm, {
    perm <- sample(fam_vec)
    pp   <- expand.grid(i=seq_along(perm),j=seq_along(perm)) %>%
      dplyr::filter(i<j) %>%
      dplyr::mutate(
        within=perm[i]==perm[j],
        sim=mapply(function(a,b) sim_mat[phenos[a],phenos[b]],i,j))
    mean(pp$sim[pp$within],na.rm=TRUE)-
      mean(pp$sim[!pp$within],na.rm=TRUE)
  })
  perm_p <- mean(null_diffs>=obs_diff)
  list(pair_df=pair_df,
       perm_p=perm_p,obs_diff=obs_diff)
}


# ================================================================
# EXTENDED FIGURE: Family similarity violins + consistency bars
# ================================================================
# Layout:
#   A: Violins - with outliers     (FLAMES | CT | GO)
#   B: Bar charts - with outliers  (FLAMES | CT | GO)
#   C: Violins - outliers excluded
#   D: Bar charts - outliers excluded
# ================================================================
ts("Building Extended Figure (A-D)...")

# ---- Violin helper ----
# Returns a single combined violin (within vs between family)
# for one analysis type and one version. x-axis = pair type.
make_violin_single <- function(stats_obj, analysis_label, y_label) {
  if (is.null(stats_obj)) {
    return(ggplot() +
             annotate("text", x=0.5, y=0.5,
                      label=paste(analysis_label,"— insufficient data"),
                      size=3) +
             theme_void())
  }
  pd <- stats_obj$pair_df %>%
    dplyr::mutate(pair_type = factor(pair_type,
                                     levels=c("Within family",
                                              "Between family")))
  annot_label <- sprintf(
    "\u0394=%.3f\nPerm p=%s",
    stats_obj$obs_diff,
    formatC(stats_obj$perm_p,   format="e", digits=2))
  
  ggplot(pd, aes(x=pair_type, y=similarity, fill=pair_type)) +
    geom_violin(alpha=0.55, colour=NA) +
    geom_boxplot(width=0.12, outlier.size=0.4, colour="grey30") +
    geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
    annotate("text", x=1.5, y=Inf, label=annot_label,
             vjust=1.3, size=2.6) +
    scale_fill_manual(
      values=c("Within family"="#B2182B","Between family"="#2166AC"),
      guide="none") +
    theme_bw(base_size=BASE_SIZE) +
    theme(axis.text.x  = element_text(size=8),
          plot.title   = element_text(size=9, face="bold"),
          plot.margin  = margin(4,6,4,4)) +
    labs(title=analysis_label, x=NULL, y=y_label)
}

# ---- Bar helper ----
make_consistency_single <- function(sim_mat, family_label,
                                    version_label, x_label) {
  if (is.null(sim_mat)) return(NULL)
  phenos <- rownames(sim_mat)
  df <- expand.grid(a=phenos, b=phenos, stringsAsFactors=FALSE) %>%
    dplyr::filter(a < b) %>%
    dplyr::mutate(
      fam_a      = pheno_family_vec[a],
      fam_b      = pheno_family_vec[b],
      within     = fam_a == fam_b,
      similarity = mapply(function(x,y) sim_mat[x,y], a, b)) %>%
    dplyr::filter(within, !is.na(similarity)) %>%
    dplyr::group_by(fam_a) %>%
    dplyr::summarise(mean_sim = mean(similarity, na.rm=TRUE),
                     .groups="drop") %>%
    dplyr::rename(family=fam_a) %>%
    dplyr::mutate(
      family  = factor(family, levels=rev(family_order)),
      version = version_label)
  
  ggplot(df, aes(x=mean_sim, y=family)) +
    geom_col(width=0.6, alpha=0.88,
             fill=if (grepl("excl",version_label,ignore.case=TRUE))
               "#B2182B" else "#AAAAAA") +
    scale_x_continuous(expand=expansion(mult=c(0, 0.05))) +
    theme_bw(base_size=BASE_SIZE) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y        = element_text(size=9),
          plot.title         = element_text(size=9, face="bold"),
          plot.margin        = margin(4,6,4,4)) +
    labs(title=family_label, x=x_label, y=NULL)
}

# ---- Build row of three violin panels ----
make_violin_row <- function(s_flames, s_ct, s_go, panel_label) {
  pv1 <- make_violin_single(s_flames, "FLAMES gene overlap",
                            "Jaccard similarity")
  pv2 <- make_violin_single(s_ct,     "Cell type enrichment",
                            "Pearson r")
  pv3 <- make_violin_single(s_go,     "GO enrichment",
                            "Spearman r")
  plot_grid(pv1, pv2, pv3, nrow=1,
            labels=c(panel_label,"",""), label_size=12)
}

# ---- Build row of three bar panels ----
make_bar_row <- function(jac_mat, ct_mat, go_mat,
                         version_label, panel_label) {
  pb1 <- make_consistency_single(
    jac_mat, "FLAMES gene overlap", version_label,
    "Mean within-family Jaccard")
  pb2 <- make_consistency_single(
    ct_mat,  "Cell type enrichment", version_label,
    "Mean within-family Pearson r")
  pb3 <- make_consistency_single(
    go_mat,  "GO enrichment", version_label,
    "Mean within-family Spearman r")
  # Handle NULLs
  panels <- list(pb1, pb2, pb3)
  panels[sapply(panels, is.null)] <- list(
    ggplot() + theme_void())
  plot_grid(panels[[1]], panels[[2]], panels[[3]], nrow=1,
            labels=c(panel_label,"",""), label_size=12)
}

# ---- Compute similarity matrices ----
jac_mat_full <- make_jac_mat(df_full)
# ct_mat_full / go_mat_full already computed earlier (used by Supp1)

# ---- Similarity stats for violins ----
compute_sim_stats <- function(sim_mat, version_label) {
  # Was previously a separate implementation that shuffled raw
  # similarity VALUES into pools of the observed within/between
  # sizes, rather than shuffling family LABELS as described in the
  # methods text (and as compute_family_similarity, used for the
  # Supp2 violin check, correctly does). Delegating here so the
  # headline Extended Figure numbers use the same, text-matching
  # permutation method as the supplementary check.
  compute_family_similarity(sim_mat, version_label)
}

sim_stats <- list(
  flames_full = compute_sim_stats(jac_mat_full, "All parcels"),
  ct_full     = compute_sim_stats(ct_mat_full,   "All parcels"),
  go_full     = compute_sim_stats(go_mat_full,   "All parcels"))

# ---- SENSITIVITY: leave-one-family-out ----
# Brain size has 4-5 members vs 1-2 for every other family, so it
# alone contributes the large majority of within-family pairs
# (C(4,2)=6 or C(5,2)=10, out of ~9-13 total). A significant overall
# result is therefore consistent with either (a) family structure
# holding broadly, or (b) Brain size's internal similarity alone
# driving the effect, with the other families contributing little.
# Re-running the same permutation test with Brain size excluded
# distinguishes these: if within > between still holds among the
# remaining families, the effect is not solely a Brain size artifact.
run_loo_check <- function(sim_mat, exclude_family = "Brain size") {
  if (is.null(sim_mat)) return(NULL)
  phenos <- rownames(sim_mat)
  keep <- phenos[pheno_family_vec[phenos] != exclude_family]
  if (length(keep) < 4) return(NULL)  # too few remaining to test meaningfully
  sub_mat <- sim_mat[keep, keep, drop = FALSE]
  compute_family_similarity(sub_mat, paste0("excl. ", exclude_family))
}

loo_stats <- list(
  flames_loo = run_loo_check(jac_mat_full),
  ct_loo     = run_loo_check(ct_mat_full),
  go_loo     = run_loo_check(go_mat_full))

cat("\n-- Leave-one-family-out sensitivity (Brain size excluded) --\n")
for (nm in names(loo_stats)) {
  s <- loo_stats[[nm]]
  if (is.null(s)) {
    cat(sprintf("  %s: insufficient remaining data for this check\n", nm))
  } else {
    cat(sprintf("  %s: obs_diff=%.3f, perm_p=%.4f\n", nm, s$obs_diff, s$perm_p))
  }
}

loo_summary <- purrr::map_dfr(names(loo_stats), function(nm) {
  s <- loo_stats[[nm]]
  if (is.null(s)) return(tibble(analysis = nm, obs_diff = NA_real_, perm_p = NA_real_))
  tibble(analysis = nm, obs_diff = s$obs_diff, perm_p = s$perm_p)
})
write.csv(loo_summary, file.path(table_dir, "SupplTable_leave_one_family_out.csv"), row.names = FALSE)

# ---- Assemble Extended Figure ----
row_A <- make_violin_row(sim_stats$flames_full,
                         sim_stats$ct_full,
                         sim_stats$go_full,
                         panel_label="A")

row_B <- make_bar_row(jac_mat_full,
                      ct_mat_full,
                      go_mat_full,
                      version_label="All parcels",
                      panel_label="B")

extended_fig <- plot_grid(
  row_A, row_B,
  ncol=1, rel_heights=c(1,1))

ggsave(file.path(out_dir, "ExtendedFigure_similarity.pdf"),
       extended_fig, width=10, height=7)
ggsave(file.path(out_dir, "ExtendedFigure_similarity.png"),
       extended_fig, width=10, height=7, dpi=200)
cat("Extended Figure saved.\n")
ts("Extended Figure done.")

# ================================================================
# SUMMARY
# ================================================================
ts("All complete.")
cat("\n", strrep("=",55), "\n")
cat("OUTPUTS — Figure_4/\n")
cat(strrep("=",55), "\n\n")
cat("Main figure:\n")
cat("  Figure4_full.pdf\n\n")
cat("Individual panels:\n")
cat("  panelA_full.pdf\n")
cat("  panelB_full.pdf\n")
cat("  panelC_full.pdf\n\n")
cat("Extended figure:\n")
cat("  ExtendedFigure_similarity.pdf / .png\n")
cat("  (A: violins | B: bars — within vs between family similarity)\n\n")
cat("Supplementary figure:\n")
cat("  Supp1_hclust_full.pdf\n\n")
cat("Tables (Figure_4/tables/):\n")
cat("  SupplTable_FLAMES_genes_full.csv\n")
cat("  SupplTable_GO_enrichment_full.csv\n")
cat("  SupplTable_CellType_enrichment_full.csv\n")
cat("  SupplTable_GO_semantic_clusters.csv\n")
cat("  SupplTable_GO_pairwise_similarity.csv\n")
cat("  celltype_enrichment_full.csv\n")
cat("  go_enrichment_significant_full.csv\n")
cat("  (go_all_*.rds in cache/ for similarity matrices)\n")