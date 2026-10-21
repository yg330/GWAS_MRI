# ===========================================================================
# Final_Figure6_full_and_GRNsupp.R
# Manuscript item: Figure 5 (final layout) and Supplementary Figure 10
#
# What this does:
#   Rebuilds the final arealisation figure with the mirrored panel A and the
#   eRegulon panel D, and the TCF7L2 / SPRY2 gene regulatory network
#   supplement.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Figure 6 (FULL, v2) + Supplementary GRN Figure
#
# Rebuilds the complete main Figure 6 with the revised Panel A and
# a new Panel D, plus a separate supplementary figure of gene
# regulatory network (GRN) diagrams for TCF7L2 and SPRY2.
#
# MAIN FIGURE 6 PANELS
#   A -- exemplar genes (FOXP1, TCF7L2, SPRY2), EACH shown as a
#        mirrored pair of views: (i) binary FDR-significant spatial
#        cluster (as in the original figure) and (ii) the lead
#        SNP's continuous per-parcel beta effect, diverging colour
#        scale. Both views as lateral+medial pairs, left hemisphere
#        (bilaterally-averaged parcellation).
#   B -- within- vs between-family Jaccard overlap of arealisation
#        gene parcels (unchanged from the original figure; rebuilt
#        from the already-computed pairwise table + a fresh
#        permutation test using the same procedure as before).
#   C -- external validation: enrichment of the 79 arealisation
#        genes against the Qian and Bhaduri candidate gene sets
#        (unchanged).
#   D -- NEW. eRegulon TF-target enrichment among the "All
#        arealisation" gene set (FDR<0.05 hits), i.e. which TFs'
#        predicted target genes are enriched for our arealisation
#        genes. TFs that are THEMSELVES one of the 79 arealisation
#        genes (tf_is_arealisation == TRUE) are flagged with a
#        filled point/asterisk -- this is the direct visual evidence
#        for the tier-one/tier-two argument in the text (a
#        regulator that is itself regulated).
#   E -- anterior-posterior morphogen pathway enrichment (Wnt/FGF/
#        Notch) by spatial cluster (this was "Panel D" in the prior
#        version of the script; relabelled E now that D is taken).
#
# SUPPLEMENTARY FIGURE: GRN diagrams
#   Two small networks -- upstream TFs (nodes) -> TCF7L2, and
#   upstream TFs -> SPRY2 -- built directly from the raw eRegulon
#   target-gene lists (Wang_eregulonlist.xlsx), NOT from the
#   pre-filtered "7 arealization TFs" text claim, so this is an
#   independent check of that claim as well as a figure. Each
#   upstream TF is coloured by whether ITS target-gene set is
#   significantly enriched for the 79 arealisation genes overall
#   (red = FDR<0.05 in the "All arealisation" eRegulon test, grey =
#   not significant), and labelled with an asterisk if that TF is
#   itself one of the 79 arealisation genes.
#
# Inputs (all pre-computed by Final_Figure6_arealizationscript.R,
# not re-derived here):
#   Outputs/Figure6_arealisation/tables/arealisation_lookup.csv
#   Outputs/Figure6_arealisation/tables/jaccard_within_between_pairs.csv
#   Outputs/Figure6_arealisation/tables/enrichment_results.csv
#   Outputs/Figure6_arealisation/tables/arealisation_binary_gsea.csv
#   Outputs/Figure6_arealisation/tables/arealisation_eregulon_enrichment.csv
#   Outputs/Figure6_arealisation/tables/arealisation_filtered.csv
#   Outputs/ExtDataFig4_PCA_master/arealisation_beta_profiles.csv
#   Data/Wang_eregulonlist.xlsx   (raw eRegulon target-gene lists,
#                                  for the GRN supplement only)
#
# Outputs:
#   Outputs/Figure6_arealisation/Figure6_arealisation_v2.pdf / .png
#   Outputs/Figure6_arealisation/SuppFig_GRN_TCF7L2_SPRY2.pdf

source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readxl)
  library(ggplot2)
  library(cowplot)
  library(ggseg)
  library(ggsegGlasser)
  library(sf)
  library(igraph)
  library(scales)
  library(RColorBrewer)
})

ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

base_dir    <- proj_dir
tables_dir  <- file.path(base_dir, "Outputs/Figure6_arealisation/tables")
out_dir     <- file.path(base_dir, "Outputs/Figure6_arealisation")
figures_dir <- file.path(out_dir, "figures")
beta_path   <- file.path(base_dir, "Outputs/ExtDataFig4_PCA_master/arealisation_beta_profiles.csv")
eregulon_xlsx <- file.path(base_dir, "Data/Wang_eregulonlist.xlsx")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

BASE_SIZE  <- 10
TITLE_SIZE <- 10
LABEL_SIZE <- 9

theme_fig6 <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          plot.title = element_text(size = TITLE_SIZE, face = "bold"),
          strip.background = element_rect(fill = "grey95", colour = NA))
}

phenotype_colours <- c(
  CT = "#0072B2", FA = "#C0006E", FI = "#7DDCB8",
  GC = "#6B4800", GMV = "#003D28", IC = "#2DC48E",
  ICVF = "#E8A8CB", ISOVF = "#A8D8F0", MC = "#F5C842",
  MD = "#0D3F5F", OD = "#D55E00", SA = "#007A52"
)

family_map <- c(SA = "BrainSize", GMV = "BrainSize",
                 GC = "Curvature", MC = "Curvature", FI = "Curvature", IC = "Curvature",
                 CT = "Thickness",
                 FA = "Coherence", ICVF = "Coherence",
                 MD = "Diffusivity", ISOVF = "Diffusivity",
                 OD = "OrientationDispersion")

# ================================================================
# PART 0: load all pre-computed tables
# ================================================================
ts("Loading pre-computed tables...")
lookup           <- read.csv(file.path(tables_dir, "arealisation_lookup.csv"), stringsAsFactors = FALSE)
pair_df_jac       <- read.csv(file.path(tables_dir, "jaccard_within_between_pairs.csv"), stringsAsFactors = FALSE)
enrichment_results <- read.csv(file.path(tables_dir, "enrichment_results.csv"), stringsAsFactors = FALSE)
pathway_combined  <- read.csv(file.path(tables_dir, "arealisation_binary_gsea.csv"), stringsAsFactors = FALSE)
eregulon_enrich   <- read.csv(file.path(tables_dir, "arealisation_eregulon_enrichment.csv"), stringsAsFactors = FALSE)
beta_profiles     <- read.csv(beta_path, stringsAsFactors = FALSE)

# ================================================================
# PART 1: Glasser geometry (shared by Panel A and nothing else)
# ================================================================
ts("Loading Glasser geometry...")
glasser_atlas <- ggsegGlasser::glasser()

lh_sf_lat <- glasser_atlas$data$sf %>%
  dplyr::filter(grepl("^lh_L_", label), view == "lateral") %>%
  mutate(parcel = gsub("^lh_L_(.*)_ROI$", "\\1", label))
lh_sf_med <- glasser_atlas$data$sf %>%
  dplyr::filter(grepl("^lh_L_", label), view == "medial") %>%
  mutate(parcel = gsub("^lh_L_(.*)_ROI$", "\\1", label))

# ================================================================
# PART 2: PANEL A -- mirrored exemplar brain plots
# ================================================================
ts("Building Panel A (combined beta + FLAMES-outline exemplars)...")
# Combined into ONE brain map per gene: continuous fill = lead-SNP
# beta (as before), with a bold black outline traced around the
# FDR-significant FLAMES cluster parcels overlaid on top. This
# collapses the previous 4-image-per-gene layout (binary pair +
# gradient pair) down to a single lateral+medial pair per gene,
# freeing enough vertical space to show a 4th exemplar gene.
exemplar_genes <- c("FOXP1", "TCF7L2", "SPRY2", "BCL11B")

exemplar_rows <- lookup %>%
  dplyr::filter(symbol %in% exemplar_genes) %>%
  dplyr::filter(!(symbol == "SPRY2" & phenotype_type != "SA")) %>%
  group_by(symbol) %>%
  slice_min(fdr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(symbol = factor(symbol, levels = exemplar_genes)) %>%
  arrange(symbol)
stopifnot(nrow(exemplar_rows) == length(exemplar_genes))

# Dissolves the hit-parcel set into its outer boundary (via union +
# st_boundary) so the outline traces the CLUSTER perimeter rather
# than every individual parcel edge -- much cleaner at this scale.
make_combined_brain <- function(hit_p, sf_df, beta_df, title_str, limits) {
  plot_df <- sf_df %>% left_join(beta_df, by = "parcel")
  outline_sf <- plot_df %>% dplyr::filter(parcel %in% hit_p)

  p <- ggplot(plot_df) +
    geom_sf(aes(geometry = geometry, fill = beta), colour = "grey85", linewidth = 0.08) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                          midpoint = 0, limits = limits, oob = scales::squish, name = expression(beta[lead~SNP]))

  if (nrow(outline_sf) > 0) {
    outline_union <- sf::st_union(outline_sf$geometry) %>% sf::st_boundary()
    p <- p + geom_sf(data = outline_union, fill = NA, colour = "black", linewidth = 0.9)
  }

  p + coord_sf(expand = FALSE) +
    theme_void(base_size = BASE_SIZE) +
    labs(title = title_str) +
    theme(plot.title = element_text(size = LABEL_SIZE, hjust = 0.5, face = "bold", margin = margin(b = 0)),
          plot.margin = margin(0, 0, 0, 0), legend.position = "none")
}

make_gradient_legend <- function(beta_df, limits) {
  ggplot(beta_df, aes(x = 1, y = 1, fill = beta)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                          midpoint = 0, limits = limits, oob = scales::squish, name = expression(beta[lead~SNP])) +
    guides(fill = guide_colourbar(barwidth = 4, barheight = 0.4, title.position = "top", title.hjust = 0.5)) +
    theme_void(base_size = 7) +
    theme(legend.position = "bottom", legend.title = element_text(size = 7), legend.text = element_text(size = 6))
}

rowA_plots <- pmap(
  list(exemplar_rows$symbol, exemplar_rows$phenotype_type, exemplar_rows$hit_parcels),
  function(sym, phen, hp_str) {
    hit_p <- str_split(hp_str, ";")[[1]]

    beta_df <- beta_profiles %>%
      dplyr::filter(symbol == as.character(sym), phenotype_type == phen) %>%
      dplyr::select(parcel, beta)
    lim <- c(-1, 1) * max(abs(beta_df$beta), na.rm = TRUE)

    # Title pared down to just the gene symbol -- the fill/outline
    # legend (beta scale, "black outline = FDR-sig FLAMES cluster")
    # belongs in the figure caption, not repeated on every panel.
    p_lat <- make_combined_brain(hit_p, lh_sf_lat, beta_df, as.character(sym), lim)
    p_med <- make_combined_brain(hit_p, lh_sf_med, beta_df, "", lim)
    p_legend <- make_gradient_legend(beta_df, lim)

    brain_pair <- plot_grid(p_lat, p_med, nrow = 1, align = "hv", rel_widths = c(1, 1))
    plot_grid(brain_pair, p_legend, ncol = 1, rel_heights = c(1, 0.16))
  })

# 2 rows x 2 cols (plot_grid fills row-major by default, so order
# is FOXP1|TCF7L2 on row 1, SPRY2|BCL11B on row 2, matching
# exemplar_genes order above). scale slightly >1 lets each brain
# image fill more of its cell, eating into the dead space that a
# fixed-aspect coord_sf() map otherwise leaves inside an
# equal-width grid cell.
panelA_ncol <- 2
panelA_nrow <- ceiling(length(exemplar_genes) / panelA_ncol)
panel_A <- plot_grid(plotlist = rowA_plots, nrow = panelA_nrow, ncol = panelA_ncol,
                      align = "hv", scale = 1.08)
panelA_w <- 9
panelA_h <- 1.7 * panelA_nrow
# add the "A" tag as a title-row overlay
panel_A <- plot_grid(ggdraw() + draw_label("A", fontface = "bold", size = 11, x = 0, hjust = 0),
                      panel_A, ncol = 1, rel_heights = c(0.05, 1))

# ================================================================
# PART 3: PANEL B -- within/between family Jaccard (re-permuted)
# ================================================================
ts("Building Panel B (Jaccard)...")
set.seed(42)
N_PERM_JAC <- 10000

within_vals  <- pair_df_jac$jaccard[pair_df_jac$within]
between_vals <- pair_df_jac$jaccard[!pair_df_jac$within]
obs_diff_jac <- mean(within_vals, na.rm = TRUE) - mean(between_vals, na.rm = TRUE)

phenos_present_jac <- names(family_map)
null_diffs_jac <- replicate(N_PERM_JAC, {
  perm_fam <- setNames(sample(family_map), phenos_present_jac)
  w <- perm_fam[pair_df_jac$phen_a] == perm_fam[pair_df_jac$phen_b]
  mean(pair_df_jac$jaccard[w], na.rm = TRUE) - mean(pair_df_jac$jaccard[!w], na.rm = TRUE)
})
perm_p_jac <- mean(null_diffs_jac >= obs_diff_jac)

jaccard_df <- bind_rows(
  tibble(group = "Within",  jaccard = within_vals),
  tibble(group = "Between", jaccard = between_vals)
) %>% mutate(group = factor(group, levels = c("Within", "Between")))

p_B <- ggplot(jaccard_df, aes(x = group, y = jaccard, fill = group)) +
  geom_violin(alpha = 0.7, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("Within" = "#B2182B", "Between" = "#AAAAAA"), guide = "none") +
  labs(x = NULL, y = "Jaccard", tag = "B",
       subtitle = sprintf("Δ = %.3f, perm p = %s", obs_diff_jac, formatC(perm_p_jac, format = "e", digits = 2))) +
  theme_fig6(BASE_SIZE - 2) +
  theme(axis.text = element_text(size = LABEL_SIZE - 2),
        plot.tag = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = LABEL_SIZE - 2))

# ================================================================
# PART 4: PANEL C -- Qian / Bhaduri overall enrichment
# ================================================================
ts("Building Panel C (external enrichment)...")
panelC_df <- enrichment_results %>%
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
  theme(axis.text = element_text(size = LABEL_SIZE - 2), plot.tag = element_text(face = "bold", size = 11))

# ================================================================
# PART 5: PANEL D -- NEW. eRegulon TF-target enrichment ("All
# arealisation" gene set), FDR<0.05, TFs that are themselves
# arealisation genes flagged.
# ================================================================
ts("Building Panel D (eRegulon TF enrichment, new)...")

panelD_df <- eregulon_enrich %>%
  dplyr::filter(gene_set == "All arealisation", fdr < 0.05) %>%
  arrange(fdr) %>%
  slice_head(n = 20) %>%   # top 20 for readability; full list in Supp Table
  mutate(log2OR = log2(pmax(OR, 0.01)),
         tf_label = ifelse(tf_is_arealisation, paste0(tf, " *"), tf),
         tf_label = factor(tf_label, levels = rev(tf_label)))

p_D <- ggplot(panelD_df, aes(x = log2OR, y = tf_label, colour = tf_is_arealisation, size = n_overlap)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_point() +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#4393C3"),
                       labels = c(`TRUE` = "TF is itself\nan arealisation gene", `FALSE` = "TF only"),
                       name = NULL) +
  scale_size_continuous(range = c(1.5, 5), name = "n overlap") +
  labs(x = "log2 OR", y = NULL, tag = "D",
       title = "eRegulon TF-target enrichment\n(All arealisation genes, FDR<0.05)") +
  theme_fig6(BASE_SIZE - 2) +
  theme(axis.text = element_text(size = LABEL_SIZE - 2.5),
        plot.title = element_text(size = TITLE_SIZE - 1),
        legend.position = "bottom", legend.text = element_text(size = 6.5), legend.title = element_text(size = 7),
        plot.tag = element_text(face = "bold", size = 11))

# ================================================================
# PART 6: PANEL E -- AP morphogen pathway enrichment (was "D")
# ================================================================
ts("Building Panel E (AP pathway enrichment)...")
pathway_df <- pathway_combined %>% dplyr::filter(threshold == 0.50)

panelE_df <- pathway_df %>%
  mutate(label = gsub("^(GOBP_|KEGG_|REACTOME_)", "", pathway),
         label = gsub("_", " ", label), label = tools::toTitleCase(tolower(label)),
         label = str_wrap(label, width = 22),
         log2OR = log2(pmax(OR, 0.01)),
         sig = case_when(fdr < 0.05 ~ "FDR < 0.05", p_value < 0.05 ~ "p < 0.05", TRUE ~ "ns"),
         group = factor(group, levels = c("frontal", "posterior", "intermediate")))

sig_pathways_e <- panelE_df %>% dplyr::filter(p_value < 0.05) %>% pull(pathway) %>% unique()
panelE_df <- panelE_df %>% dplyr::filter(pathway %in% sig_pathways_e)

p_E <- ggplot(panelE_df, aes(x = group, y = reorder(label, log2OR), size = OR, colour = sig)) +
  geom_point() +
  scale_colour_manual(values = c("FDR < 0.05" = "#B2182B", "p < 0.05" = "#F46D43", "ns" = "#AAAAAA"), name = NULL) +
  scale_size_continuous(range = c(2, 7), name = "Odds Ratio") +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, tag = "E", title = "Pathway enrichment\n(>50% territory)") +
  theme_fig6(BASE_SIZE - 1) +
  theme(axis.text = element_text(size = LABEL_SIZE - 1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        plot.title = element_text(size = TITLE_SIZE - 1),
        legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 11))

# ================================================================
# PART 7: assembly -- Figure 6 (A on top, B/C/D/E in a bottom row)
# ================================================================
ts("Assembling Figure 6...")
ggsave(file.path(figures_dir, "PanelA_exemplars_mirrored.pdf"), panel_A, width = panelA_w, height = panelA_h)
ggsave(file.path(figures_dir, "PanelB_jaccard.pdf"), p_B, width = 3, height = 3)
ggsave(file.path(figures_dir, "PanelC_DEG_overall.pdf"), p_C, width = 4, height = 2.5)
ggsave(file.path(figures_dir, "PanelD_eregulon_TF.pdf"), p_D, width = 4.5,
       height = max(3, 0.3 * nrow(panelD_df) + 1.8), limitsize = FALSE)
ggsave(file.path(figures_dir, "PanelE_pathway.pdf"), p_E, width = 4,
       height = max(3, 0.5 * n_distinct(panelE_df$label) + 1.5), limitsize = FALSE)

bottom_row_h <- 4.5
# Panel A and the bottom row share one total width when stacked
# below, so this is the main lever on how "wide" the assembled
# figure looks -- reduced from 16 (which stretched Panel A's now-
# compact brain grid back out with dead space) to 13, the
# narrowest width that still keeps Panel D's TF labels and Panel
# E's pathway labels legible.
fig6_w <- max(panelA_w, 13)
fig6_h <- panelA_h + bottom_row_h

# B and C stacked in a single first column (B gets 60% of the
# column's height, C gets 40%) rather than each taking a full
# column of its own -- neither needs that much width, and stacking
# frees room for D/E to breathe.
col_BC <- plot_grid(p_B, p_C, ncol = 1, rel_heights = c(0.6, 0.4), align = "v")
bottom_row <- plot_grid(col_BC, p_D, p_E, nrow = 1, rel_widths = c(1, 1.3, 1.3))
fig6 <- plot_grid(panel_A, bottom_row, ncol = 1, rel_heights = c(panelA_h, bottom_row_h))

ggsave(file.path(out_dir, "Figure6_arealisation_v2.pdf"), fig6, width = 0.8*fig6_w, height = 0.8*fig6_h, limitsize = FALSE)
ggsave(file.path(out_dir, "Figure6_arealisation_v2.png"), fig6, width = fig6_w, height = fig6_h, dpi = 300, limitsize = FALSE)
ts("Saved main figure:", file.path(out_dir, "Figure6_arealisation_v2.pdf"))

# ================================================================
# PART 8: SUPPLEMENTARY FIGURE -- eGRN, SCENIC+-style radial network
# ================================================================
# Reformatted to match the SCENIC+ eGRN convention (e.g. Gonzalez-
# Blas et al. 2023 Fig 5e / your reference image): a small number of
# HUB TFs sit near the centre, with straight coloured edges (one
# colour per hub) radiating out to their TARGET GENES, which are
# arranged and labelled around the outer circumference. Hub TFs
# that are themselves a target of another hub (regulatory cross-
# talk) sit both centrally AND appear as a labelled point on the
# circumference.
#
# Hub TFs used here: the union of TCF7L2's and SPRY2's significant
# (FDR<0.05 in the "All arealisation" eRegulon test) upstream
# regulators -- i.e. exactly the TFs already named in the manuscript
# text for these two genes, not an arbitrary top-N cut. Target genes
# shown are each hub's predicted eRegulon targets RESTRICTED to the
# 79 arealisation genes (not the full genome-wide target list), so
# this is a focused "arealisation eGRN", not the full eRegulon.
ts("Building eGRN supplementary figure (SCENIC+ radial format, two-tier)...")

areal_genes_79 <- unique(lookup$symbol)

sheets_er <- readxl::excel_sheets(eregulon_xlsx)
eregulon_raw <- map_dfr(sheets_er, function(sh) {
  df <- readxl::read_excel(eregulon_xlsx, sheet = sh,
                            col_names = c("eregulon", "eregulon_n", "target_genes"), skip = 1)
  df %>% dplyr::filter(!is.na(eregulon), !is.na(target_genes)) %>%
    mutate(tf = str_remove(eregulon, "_extended.*|_\\+.*|_-.*"),
           targets = str_split(target_genes, ",\\s*")) %>%
    dplyr::select(eregulon, tf, targets)
})

tf_targets_flat <- eregulon_raw %>%
  group_by(tf) %>%
  summarise(targets = list(unique(unlist(targets))), .groups = "drop")
all_tfs_with_eregulon <- unique(tf_targets_flat$tf)

targets_of <- function(tf) {
  idx <- which(tf_targets_flat$tf == tf)
  if (length(idx) == 0) return(character(0))
  tf_targets_flat$targets[[idx]]
}
regulators_of <- function(gene) {
  tf_targets_flat %>% dplyr::filter(map_lgl(targets, ~ gene %in% .x)) %>% pull(tf) %>% unique()
}

all_areal_enrich <- eregulon_enrich %>% dplyr::filter(gene_set == "All arealisation")
sig_tfs   <- all_areal_enrich %>% dplyr::filter(fdr < 0.05) %>% pull(tf)

# ---- TIER 1: hub TFs found via TCF7L2/SPRY2's significant regulators ----
tcf7l2_tfs_sig <- intersect(regulators_of("TCF7L2"), sig_tfs)
spry2_tfs_sig  <- intersect(regulators_of("SPRY2"),  sig_tfs)
tier1_hubs <- union(tcf7l2_tfs_sig, spry2_tfs_sig)

tier1_targets <- setNames(lapply(tier1_hubs, function(tf) setdiff(intersect(targets_of(tf), areal_genes_79), tf)), tier1_hubs)
tier1_edges <- imap_dfr(tier1_targets, function(targs, tf) if (length(targs)) tibble(hub = tf, target = targs, tier = 1L) else NULL)
tier1_leaves <- unique(tier1_edges$target)

cat(sprintf("  Tier 1: %d hub TFs (TCF7L2's %d + SPRY2's %d significant regulators): %s\n",
            length(tier1_hubs), length(tcf7l2_tfs_sig), length(spry2_tfs_sig), paste(tier1_hubs, collapse = ", ")))

# ---- TIER 2: any tier-1 TARGET that is itself a TF with its own eRegulon ----
# ("are any of the target genes TFs themselves? add their downstream
# targets too") -- added regardless of the FDR<0.05 pooled-test filter,
# since the point here is structural completeness of the cascade, not
# a significance claim; each tier-2 hub's significance status is
# reported to the console and can be re-applied if you'd rather gate it.
tier2_hubs <- setdiff(intersect(tier1_leaves, all_tfs_with_eregulon), tier1_hubs)
tier2_targets <- setNames(lapply(tier2_hubs, function(tf) setdiff(intersect(targets_of(tf), areal_genes_79), tf)), tier2_hubs)
tier2_edges <- imap_dfr(tier2_targets, function(targs, tf) if (length(targs)) tibble(hub = tf, target = targs, tier = 2L) else NULL)

cat(sprintf("  Tier 2: %d additional hub TFs, found among tier-1 targets (%s of which are FDR<0.05 in the pooled test): %s\n",
            length(tier2_hubs), sum(tier2_hubs %in% sig_tfs), paste(tier2_hubs, collapse = ", ")))

edge_df  <- bind_rows(tier1_edges, tier2_edges) %>% distinct(hub, target, .keep_all = TRUE)
all_hubs <- union(tier1_hubs, tier2_hubs)
hub_tier <- setNames(ifelse(all_hubs %in% tier1_hubs, 1L, 2L), all_hubs)
leaf_genes <- setdiff(unique(edge_df$target), all_hubs)   # pure leaves: targeted, never themselves a hub here

self_reg_hubs <- all_hubs[map_lgl(all_hubs, ~ .x %in% intersect(targets_of(.x), areal_genes_79))]

cat(sprintf("  Combined: %d hub TFs, %d pure target genes, %d edges. Self-loops: %s\n",
            length(all_hubs), length(leaf_genes), nrow(edge_df),
            if (length(self_reg_hubs)) paste(self_reg_hubs, collapse = ", ") else "none"))

# Every node in this network is either (a) a hub TF found via the
# TCF7L2/SPRY2 regulator search, which may or may not itself be one
# of the 79 arealisation genes, or (b) a target gene, which -- by
# construction (targets_of(hub) INTERSECT areal_genes_79) -- always
# IS one of the 79. So node fill only varies among the hub set.
is_areal_gene <- function(gene) gene %in% areal_genes_79
areal_fill  <- "#B2182B"
other_fill  <- "white"

# ---- Layout: leaves on outer ring; tier-1 hubs on inner ring;
# tier-2 hubs on a middle ring -- visually separates "regulators
# found via TCF7L2/SPRY2" from "regulators found because a target
# turned out to be a TF itself", i.e. shows the cascade depth. ----
hub_sizes <- map_int(setNames(lapply(all_hubs, function(h) intersect(targets_of(h), areal_genes_79)), all_hubs), length)
primary_hub_of <- function(gene) {
  hubs_reg <- edge_df$hub[edge_df$target == gene]
  hubs_reg[which.min(hub_sizes[hubs_reg])]
}
leaf_df <- tibble(gene = sort(leaf_genes),
                   primary_hub = map_chr(sort(leaf_genes), primary_hub_of)) %>%
  arrange(primary_hub, gene) %>%
  mutate(angle = 2 * pi * (row_number() - 1) / n(),
         x = cos(angle), y = sin(angle))

# Iterative angular relaxation so labels/points at a shared radius
# don't sit on top of each other: nudges any pair closer than
# min_gap apart, circularly, over a few hundred iterations.
relax_angles <- function(angles, min_gap, iters = 400) {
  a <- angles %% (2 * pi)
  n <- length(a)
  if (n < 2) return(a)
  for (it in seq_len(iters)) {
    o <- order(a)
    a_sorted <- a[o]
    gaps <- c(diff(a_sorted), a_sorted[1] + 2 * pi - a_sorted[n])
    moved <- FALSE
    for (i in seq_len(n)) {
      if (gaps[i] < min_gap) {
        push <- (min_gap - gaps[i]) / 2
        j <- if (i == n) 1 else i + 1
        a_sorted[i] <- (a_sorted[i] - push) %% (2 * pi)
        a_sorted[j] <- (a_sorted[j] + push) %% (2 * pi)
        moved <- TRUE
      }
    }
    a[o] <- a_sorted
    if (!moved) break
  }
  a
}

build_hub_ring <- function(hubs, radius, min_gap) {
  if (length(hubs) == 0) return(tibble(hub = character(0)))
  raw_angle <- map_dbl(hubs, function(h) {
    targs <- edge_df$target[edge_df$hub == h]
    leaf_ang <- leaf_df$angle[leaf_df$gene %in% targs]
    if (length(leaf_ang) > 0) atan2(mean(sin(leaf_ang)), mean(cos(leaf_ang))) else 0
  })
  final_angle <- relax_angles(raw_angle, min_gap = min_gap)
  tibble(hub = hubs, angle = final_angle, radius = radius, x = radius * cos(final_angle), y = radius * sin(final_angle))
}

tier1_df <- build_hub_ring(tier1_hubs, radius = 0.32, min_gap = 2 * pi / (length(tier1_hubs) * 1.4))
tier2_df <- build_hub_ring(tier2_hubs, radius = 0.62, min_gap = 2 * pi / max(length(tier2_hubs) * 1.4, 6))
hub_df <- bind_rows(tier1_df, tier2_df)

node_pos <- bind_rows(hub_df %>% dplyr::select(gene = hub, x, y),
                       leaf_df %>% dplyr::select(gene, x, y))

hub_palette <- setNames(
  scales::hue_pal(l = 55, c = 100)(length(all_hubs)),
  all_hubs
)

focal_genes <- c("TCF7L2", "SPRY2")

pdf(file.path(out_dir, "SuppFig_GRN_TCF7L2_SPRY2.pdf"), width = 16, height = 16)
par(mar = c(1, 1, 4, 1))
plot(NULL, xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35), asp = 1, axes = FALSE, xlab = "", ylab = "",
     main = paste0("Arealisation eGRN: hub TFs and their target arealisation genes\n",
                    "(two-tier SCENIC+-style radial layout; inner ring = TFs found via TCF7L2/SPRY2; ",
                    "middle ring = TFs found among their targets; edge colour = source hub)"))

# edges, coloured by source hub, drawn first
for (i in seq_len(nrow(edge_df))) {
  h <- edge_df$hub[i]; g <- edge_df$target[i]
  p1 <- node_pos[node_pos$gene == h, ]; p2 <- node_pos[node_pos$gene == g, ]
  if (nrow(p1) == 0 || nrow(p2) == 0) next
  segments(p1$x, p1$y, p2$x, p2$y, col = scales::alpha(hub_palette[h], 0.55), lwd = 1.1)
}

# leaf points (always filled -- every leaf IS an arealisation gene by construction) + radial labels
leaf_is_focal <- leaf_df$gene %in% focal_genes
points(leaf_df$x, leaf_df$y, pch = 21, bg = areal_fill,
       col = "grey40", cex = ifelse(leaf_is_focal, 1.9, 1.2), lwd = ifelse(leaf_is_focal, 1.8, 1))

label_r <- 1.12
for (i in seq_len(nrow(leaf_df))) {
  ang_deg <- leaf_df$angle[i] * 180 / pi
  flip <- ang_deg > 90 && ang_deg < 270
  srt_val <- if (flip) ang_deg + 180 else ang_deg
  adj_val <- if (flip) c(1, 0.5) else c(0, 0.5)
  is_foc <- leaf_df$gene[i] %in% focal_genes
  text(label_r * cos(leaf_df$angle[i]), label_r * sin(leaf_df$angle[i]), labels = leaf_df$gene[i],
       srt = srt_val, adj = adj_val, cex = if (is_foc) 0.85 else 0.62,
       font = if (is_foc) 2 else 1, col = "black")
}

# hub points -- fill = red if this TF is itself one of the 79
# arealisation genes, white otherwise; border = hub identity colour
# (matches edge colour, for tracing); tier-2 hubs drawn slightly
# smaller/thinner to distinguish from tier-1 at a glance.
hub_plot_df <- hub_df %>%
  mutate(tier = hub_tier[hub], is_areal = map_lgl(hub, is_areal_gene),
         fill = ifelse(is_areal, areal_fill, other_fill),
         is_focal = hub %in% focal_genes,
         cex_pt = ifelse(is_focal, 3.0, ifelse(tier == 1, 2.3, 1.9)),
         lwd_pt = ifelse(tier == 1, 2.2, 1.4))
points(hub_plot_df$x, hub_plot_df$y, pch = 21, bg = hub_plot_df$fill,
       col = hub_palette[hub_plot_df$hub], cex = hub_plot_df$cex_pt, lwd = hub_plot_df$lwd_pt)
text(hub_plot_df$x, hub_plot_df$y,
     labels = paste0(hub_plot_df$hub, ifelse(hub_plot_df$hub %in% self_reg_hubs, " ⟲", "")),
     cex = ifelse(hub_plot_df$is_focal, 0.85, 0.68), font = 2, col = "black")

legend("bottomleft",
       legend = c("Hub TF -- is itself an arealisation gene", "Hub TF -- regulator only (not in the 79-gene set)",
                  "Target gene (always an arealisation gene)", "TCF7L2 / SPRY2 (focal genes)",
                  "Inner ring = tier-1 hub (found via TCF7L2/SPRY2)", "Middle ring = tier-2 hub (found among tier-1 targets)",
                  "⟲ = predicted self-loop", "Edge colour = source hub identity"),
       pch = c(21, 21, 21, 21, NA, NA, NA, NA),
       pt.bg = c(areal_fill, other_fill, areal_fill, areal_fill, NA, NA, NA, NA),
       pt.cex = c(1.6, 1.6, 1.2, 2, NA, NA, NA, NA),
       col = c("grey30", "grey30", "grey40", "grey30", NA, NA, NA, NA),
       bty = "n", cex = 0.72)

dev.off()
ts("Saved eGRN supplement:", file.path(out_dir, "SuppFig_GRN_TCF7L2_SPRY2.pdf"))

cat("\n=== ALL DONE ===\n")
cat("Main figure: ", file.path(out_dir, "Figure6_arealisation_v2.pdf"), "\n")
cat("eGRN supp:   ", file.path(out_dir, "SuppFig_GRN_TCF7L2_SPRY2.pdf"), "\n")
