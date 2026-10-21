# Shared settings, sourced at the top of every script.
# proj_dir is the repository root. By default it is the current working directory,
# so setwd() to the repo root first, or point GWAS_MRI_DIR at it.
proj_dir <- normalizePath(Sys.getenv("GWAS_MRI_DIR", unset = getwd()), mustWork = FALSE)
if (!file.exists(file.path(proj_dir, "config.R"))) {
  stop("Run from the repository root (config.R not found in ", proj_dir, ")")
}
dir.create(file.path(proj_dir, "Outputs"), showWarnings = FALSE)
