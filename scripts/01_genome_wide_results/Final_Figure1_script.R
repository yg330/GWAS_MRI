# ===========================================================================
# Final_Figure1_script.R
# Manuscript item: Figure 1
#
# What this does:
#   Builds Figure 1: (A) Miami plot of genome-wide significant loci for UKB
#   (top) and the UKB-ABCD meta-analysis (bottom), (B) SNP heritability by
#   measure, (C) UKB/ABCD effect-direction concordance by measure. Panel C
#   reads Data/UKB_sig_info.txt without a threshold filter, so it shows
#   genome-wide (P <= 5e-8) concordance; the experiment-wide numbers are
#   quoted in the text.
#
# Run from the repository root (setwd() there, or open the folder as a project).
# Paths come from config.R via proj_dir; inputs are in Data/, results go to Outputs/.
# ===========================================================================
#
# ---- Notes from the original script header ----
# Figure 1 - multipanel
#   Panel A: Miami plot (UKB top / Meta-analysis bottom)
#   Panel B: Heritability violin plot, by MRI phenotype
#   Panel C: UKB-ABCD concordance/discordance bar chart, by phenotype
#
# Phenotype colours: imaging_gen_v2 measure_colors palette (used in all 3 panels)
# Typography:        Helvetica throughout, Nature-scale sizing

source("config.R")

library(data.table)
library(ggplot2)
library(stringr)
library(patchwork)
library(colorspace)   # for lighten() used in Panel C
library(ggrastr)      # rasterises Panel A's point cloud so the combined
                       # Figure 1 PDF stays a manageable file size -- Panels
                       # B and C, all axes/text/legends, remain vector.
                       # Needs the `ragg` package installed as its backend
                       # (install.packages("ragg")) -- ragg is self-contained
                       # and does NOT require system Cairo/X11, unlike the
                       # default grDevices::png(), which is why a plain
                       # ggsave(..., dpi=600) PNG came back blank too.

out_dir <- file.path(proj_dir, "Outputs/Figure_1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Shared phenotype colour palette (imaging_gen_v2) ----
measure_colors <- c(
  CT             = "#0072B2",
  FA             = "#C0006E",
  FI             = "#7DDCB8",
  GC             = "#6B4800",
  GMV            = "#003D28",
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

# "Vol" (MRI_scale == "SubCor") lumps together corpus callosum segments,
# ventricles, and other subcortical structures. Split it into those three
# constituents by region name rather than plotting it as one category.
classify_subcor_region <- function(region) {
  fifelse(
    grepl("^CC\\.", region), "CorpusCallosum",
    fifelse(grepl("Ventricle|\\.Vent$", region), "Ventricle", "SubCorVol")
  )
}

# ---- Shared theme: Helvetica, Nature-scale text sizes ----
# Nature figures are typically set at 89 mm (single column) or
# 183 mm (double column) final width - base_size 8 with 7pt axis
# text reads cleanly once scaled down to that width at 300 dpi.
nature_theme <- theme_minimal(base_family = "Helvetica", base_size = 8) +
  theme(
    text          = element_text(colour = "black"),
    axis.text     = element_text(colour = "black", size = 7),
    axis.title    = element_text(size = 8),
    strip.text    = element_text(face = "bold", size = 8),
    legend.text   = element_text(size = 7),
    legend.title  = element_text(size = 8),
    plot.tag      = element_text(face = "bold", size = 10)
  )

# ============================================================
# PANEL A - Miami plot
# ============================================================

data1 <- fread(file.path(proj_dir, "Data/UKB_sig_info.txt"))
data2 <- data1[, c("GWAS_name", "SNP", "CHR", "POS", "UKB_P")]

data1 <- fread(file.path(proj_dir, "Data/Meta_sig_info.txt"))
data2_meta <- data1[, c("GWAS_name", "SNP", "CHR", "POS", "meta_P")]

# Count number of phenotypes each SNP is significant for in UKB
snp_pheno_count_ukb <- data2[, .(n_pheno = uniqueN(GWAS_name)), by = SNP]
data2 <- merge(data2, snp_pheno_count_ukb, by = "SNP")

# Count number of phenotypes each SNP is significant for in meta
snp_pheno_count_meta <- data2_meta[, .(n_pheno = uniqueN(GWAS_name)), by = SNP]
data2_meta <- merge(data2_meta, snp_pheno_count_meta, by = "SNP")

# Prepare UKB data - positive log10P
data2[, log10P := -log10(UKB_P)]
data2[, cohort := "UKB"]

# Prepare meta data - negative log10P for bottom panel
data2_meta[, log10P := log10(meta_P)]
data2_meta[, cohort := "Meta"]

plot_data <- rbind(
  data2[, .(SNP, CHR, POS, log10P, n_pheno, cohort)],
  data2_meta[, .(SNP, CHR, POS, log10P, n_pheno, cohort)]
)

setorder(plot_data, CHR, POS)

# --- Add a gap between consecutive chromosomes on the x-axis ---
gap <- 3e7  # bp of blank space inserted between chromosomes; adjust to taste

chr_lengths <- plot_data[, .(chr_len = max(POS)), by = CHR]
setorder(chr_lengths, CHR)
chr_lengths[, cum_start := cumsum(c(0, chr_len[-.N] + gap))]

plot_data <- merge(plot_data, chr_lengths[, .(CHR, cum_start)], by = "CHR")
plot_data[, cum_pos := POS + cum_start]

# Chromosome midpoints, for the number labels
chr_mids <- plot_data[, .(mid = mean(cum_pos)), by = CHR][order(CHR)]

# Truncate log10P for plotting
plot_data[, log10P_trunc := pmax(pmin(log10P, 200), -200)]

# --- Significance thresholds ---
genome_wide_p <- 5e-8
exp_wide_p    <- 4.65e-11

genome_wide_y <- -log10(genome_wide_p)
exp_wide_y    <- -log10(exp_wide_p)

x_max <- max(plot_data$cum_pos)
x_min <- min(plot_data$cum_pos)

panelA <- ggplot(plot_data, aes(x = cum_pos, y = log10P_trunc, colour = n_pheno)) +
  # dev = "ragg" pins the rasteriser to the ragg backend explicitly, so it
  # doesn't fall back to trying Cairo (unavailable) or base grDevices::png()
  # (which returned blank without a working Cairo/X11 stack on this
  # machine). Text, axes, thresholds and legend below are untouched and
  # stay fully vector/editable -- only this point layer is rasterised.
  geom_point_rast(size = 0.8, alpha = 0.8, raster.dpi = 300, dev = "ragg") +

  # significance thresholds only - mirrored top (UKB) and bottom (Meta)
  geom_hline(yintercept = c(genome_wide_y, -genome_wide_y),
             linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_hline(yintercept = c(exp_wide_y, -exp_wide_y),
             linetype = "dashed", colour = "firebrick3", linewidth = 0.4) +
  
  annotate("text", x = x_max, y = genome_wide_y,
           label = "Genome-wide significant (P = 5e-08)",
           hjust = 1, vjust = -0.5, size = 2.4, colour = "grey45",
           family = "Helvetica") +
  annotate("text", x = x_max, y = exp_wide_y,
           label = "Experiment-wide significant (P = 4.65e-11)",
           hjust = 1, vjust = -0.5, size = 2.4, colour = "firebrick3",
           family = "Helvetica") +
  
  scale_colour_viridis_c(option = "plasma", name = "N phenotypes") +
  
  # Chromosome number labels, placed at y = 0
  annotate("text", x = chr_mids$mid, y = 0,
           label = chr_mids$CHR, size = 2.5, vjust = 0.5,
           fontface = "bold", colour = "black", family = "Helvetica") +
  
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(
    limits = c(-200, 200),
    labels = function(x) abs(x),
    name = expression(-log[10](P))
  ) +
  
  annotate("text", x = x_min, y = 190, label = "UK Biobank",
           hjust = 0, fontface = "bold", size = 3.2, family = "Helvetica") +
  annotate("text", x = x_min, y = -190, label = "Meta-analysis",
           hjust = 0, fontface = "bold", size = 3.2, family = "Helvetica") +
  
  nature_theme +
  theme(
    panel.grid       = element_blank(),   # removes ALL horizontal + vertical gridlines
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    legend.position  = "bottom"
  )

# ============================================================
# PANEL B - Heritability violin plot, by MRI phenotype
# ============================================================

h2_data <- fread(file.path(proj_dir, "Data/heritability_table.txt"))
# expected columns: MRI_scale  MRI_pheno  region  UKB_Lambda_GC  UKB_Mean_chi_square  UKB_h2

# Split "Vol" into its three constituents by region name - everywhere below
# uses this column instead of the original MRI_pheno.
h2_data[, MRI_pheno := fifelse(MRI_pheno == "Vol",
                               classify_subcor_region(region),
                               MRI_pheno)]

# Scale labels that carry a per-region/per-structure/per-tract distribution
# (i.e. something a violin can be built from), as opposed to a single
# whole-brain "global" summary value. Rather than hardcoding which labels
# count ("regional", "SubCor", "tract", ...), treat everything that isn't
# "global" as regional-like - keeps this robust if more scale labels show
# up later.
regional_like <- setdiff(unique(h2_data$MRI_scale), "global")

# Order phenotypes by median regional h2 (descending) for a cleaner plot.
# Built from ALL phenotypes present (not just ones filtered to
# MRI_scale == "regional") - otherwise a phenotype with only global-level
# rows (e.g. Vol, if it has no per-region breakdown) would be missing from
# this vector entirely, and factor(..., levels = pheno_order) would then
# silently coerce every one of its values to NA.
pheno_order_dt <- h2_data[, .(
  med_h2    = median(UKB_h2[MRI_scale %in% regional_like], na.rm = TRUE),
  global_h2 = UKB_h2[MRI_scale == "global"][1]
), by = MRI_pheno]

# fall back to the global h2 value when there's no regional data to average
pheno_order_dt[is.na(med_h2), med_h2 := global_h2]
setorder(pheno_order_dt, -med_h2)
pheno_order <- pheno_order_dt$MRI_pheno

h2_data[, MRI_pheno := factor(MRI_pheno, levels = pheno_order)]

# Diagnostic: flag any phenotype with literally no plottable data (no
# regional-like rows AND an NA/missing global value) - ggplot's discrete
# scales drop such categories entirely by default, which looks identical
# to a colour/mapping bug but isn't one.
check_dt <- h2_data[, .(
  n_regional = sum(MRI_scale %in% regional_like),
  global_val = UKB_h2[MRI_scale == "global"][1]
), by = MRI_pheno]
flagged <- check_dt[n_regional == 0 & is.na(global_val), MRI_pheno]
if (length(flagged) > 0) {
  message("NOTE: these phenotypes have NO regional-like rows and an NA global h2 - ",
          "there is nothing to plot for them and they will not appear on the ",
          "axis regardless of drop = FALSE: ", paste(flagged, collapse = ", "))
}

pheno_cols_present <- measure_colors[levels(h2_data$MRI_pheno)]

panelB <- ggplot(h2_data[MRI_scale %in% regional_like],
                 aes(x = MRI_pheno, y = UKB_h2, fill = MRI_pheno, colour = MRI_pheno)) +
  # colour == fill -> outline is present (violins need a boundary to
  # render cleanly) but is invisible against its own fill
  geom_violin(linewidth = 0.3, trim = TRUE, alpha = 0.95) +
  geom_jitter(width = 0.08, size = 0.5, alpha = 0.35, colour = "grey20") +
  # overlay the single global h2 value per phenotype - solid diamond,
  # no separate outline aesthetic to worry about
  geom_point(data = h2_data[MRI_scale == "global"],
             aes(x = MRI_pheno, y = UKB_h2),
             shape = 18, size = 2.6, colour = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = pheno_cols_present, guide = "none") +
  scale_colour_manual(values = pheno_cols_present, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = expression(SNP~heritability~(h^2))) +
  nature_theme +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(face = "bold")
  )

# ============================================================
# PANEL C - UKB / ABCD concordance-discordance bar chart, by phenotype
# ============================================================

rep_data <- fread(file.path(proj_dir, "Data/UKB_sig_info.txt"))
# expected columns include: GWAS_name, concordance ("Concordant"/"Discordant")
#


# Phenotype code sits as the 2nd underscore-delimited token of GWAS_name,
# region is everything after that, e.g. "regional_SA_TE2p" -> phenotype
# "SA", region "TE2p"; "SubCor_Vol_CC.Anterior" -> phenotype "Vol",
# region "CC.Anterior"
gwas_parts <- str_split_fixed(rep_data$GWAS_name, "_", 3)
rep_data[, phenotype := gwas_parts[, 2]]
rep_data[, region    := gwas_parts[, 3]]

# Split "Vol" into its three constituents, same as Panel B
rep_data[, phenotype := fifelse(phenotype == "Vol",
                                classify_subcor_region(region),
                                phenotype)]


rep_data_valid <- rep_data[!is.na(concordance)]

conc_summary <- rep_data_valid[, .(
  n_total      = .N,
  n_concordant = sum(concordance == "Concordant"),
  n_discordant = sum(concordance == "Discordant")
), by = .(phenotype)]

conc_summary[, pct_concordant := 100 * n_concordant / n_total]
conc_summary[, pct_discordant := 100 * n_discordant / n_total]

# order phenotypes by % concordant, ascending (so most-concordant plots at top)
pheno_order_C <- conc_summary[order(pct_concordant), phenotype]
conc_summary[, phenotype := factor(phenotype, levels = pheno_order_C)]

# long format: one row per phenotype x side (Concordant / Discordant)
bar_data <- rbind(
  conc_summary[, .(phenotype, side = "Concordant",
                   value = pct_concordant, n = n_concordant)],
  conc_summary[, .(phenotype, side = "Discordant",
                   value = -pct_discordant, n = n_discordant)]
)

# Concordant bars = full phenotype colour; discordant bars = lightened
# version of that SAME phenotype colour, so hue still maps to phenotype
# and saturation maps to concordance
pheno_cols_C <- measure_colors[levels(bar_data$phenotype)]
fill_lookup <- c(
  setNames(pheno_cols_C, paste0(names(pheno_cols_C), "_Concordant")),
  setNames(colorspace::lighten(pheno_cols_C, 0.55),
           paste0(names(pheno_cols_C), "_Discordant"))
)
bar_data[, fill_key := paste0(phenotype, "_", side)]
bar_data[, label_hjust := ifelse(side == "Concordant", -0.15, 1.15)]

panelC <- ggplot(bar_data, aes(y = phenotype, x = value, fill = fill_key)) +
  geom_col(width = 0.7, colour = NA) +   # no outline
  geom_text(aes(label = sprintf("%.1f%% (n=%d)", abs(value), n),
                hjust = label_hjust),
            size = 2.3, family = "Helvetica", colour = "black") +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = fill_lookup, guide = "none") +
  # right-hand (concordant) side truncated to 75 - the previous 115 left a
  # lot of unused white space since concordant % doesn't reach that high
  scale_x_continuous(limits = c(-75, 100),
                     breaks = seq(-100, 75, 25),
                     labels = function(x) abs(x)) +
  labs(x = "% of loci        Discordant  |  Concordant", y = NULL) +
  nature_theme +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(face = "bold")
  )

# ============================================================
# COMBINE + SAVE
# ============================================================

figure1 <- panelA / (panelB | panelC) +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Helvetica", face = "bold", size = 10))

# NOTE: Helvetica is one of the 14 standard PDF fonts (built into every PDF
# viewer, same as Times/Courier), so the default pdf() device can write it
# directly - no Cairo dependency needed. (cairo_pdf() was tried here first
# but throws "failed to load cairo DLL" on some Windows R installations
# where the Cairo library isn't properly linked; this sidesteps that
# entirely rather than requiring a Cairo/Ghostscript install.)
ggsave(file.path(out_dir, "Figure1.pdf"), figure1,
       width = 14, height = 12, dpi = 300)