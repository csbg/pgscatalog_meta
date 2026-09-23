#!/usr/bin/env Rscript

# ============================================================
# Script: 01_Installation_API_and_perfomance_pull.R
# Project: pgscatalog_meta
# Purpose: Download PGS Catalog performance metrics and save normalized cache tables.
# Inputs: PGS Catalog API via quincunx
# Outputs: data/pgs_cache/pm_<DATE>.rds; normalized performance, sample-set, sample, and metric tables
# Run after: none
# Run before: 02_loading_PGS_matrix.R
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

# ---- Helper: install missing packages (silent) ----
.install_if_missing <- function(pkgs) {
  to_get <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_get)) install.packages(to_get, quiet = TRUE)
}
.install_if_missing(c("quincunx","dplyr","tidyr","readr","purrr","stringr"))

suppressPackageStartupMessages({
  library(quincunx)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
})

have_arrow   <- requireNamespace("arrow", quietly = TRUE)
have_jsonlite<- requireNamespace("jsonlite", quietly = TRUE)

if (have_arrow)   library(arrow)

if (!file.exists("R/load.R")) {
  stop("R/load.R not found. Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

STAMP     <- format(Sys.Date(), "%Y%m%d")
CACHE_DIR <- file.path("data", "pgs_cache")
OUT_DIR   <- file.path("results", "pgs_top15_perf")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)

log_step <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# ---- 1) Fetch from API and cache (S4 + normalized tables) ----
log_step("Pulling Performance Metrics bundle from PGS Catalog REST API…")
pm <- quincunx::get_performance_metrics(progress_bar = TRUE)

rds_path <- file.path(CACHE_DIR, paste0("pm_", STAMP, ".rds"))
saveRDS(pm, rds_path, compress = "xz")
log_step("Saved S4 cache ⇒ ", rds_path)

# Extract normalized tables (slots)
ppm       <- pm@performance_metrics                       # PPM-level (links to PGS)
pss_links <- pm@sample_sets %>% distinct(ppm_id, pss_id)  # mapping PPM→PSS
samples   <- pm@samples                                   # ancestry + sizes
classm    <- pm@pgs_classification_metrics                # AUROC / C-index
effectm   <- pm@pgs_effect_sizes                          # OR / HR / beta
otherm    <- pm@pgs_other_metrics                         # R^2 etc.

# 2) Map ancestry with the shared display mapper

# 3) Map ONCE
samples <- samples %>%
  mutate(
    ancestry_display  = to_display_cat(ancestry_category),
    ancestry_category = ancestry_display
  )

# 4) Sanity check (African should now be present)
samples %>% filter(stage == "eval") %>%
  count(ancestry_category, sort = TRUE) %>% print(n = 50)

# Write CSVs (always)
write_csv(ppm,       file.path(CACHE_DIR, paste0("ppm_", STAMP, ".csv")))
write_csv(pss_links, file.path(CACHE_DIR, paste0("pss_links_", STAMP, ".csv")))
write_csv(samples,   file.path(CACHE_DIR, paste0("samples_", STAMP, ".csv")))
write_csv(classm,    file.path(CACHE_DIR, paste0("classm_", STAMP, ".csv")))
write_csv(effectm,   file.path(CACHE_DIR, paste0("effectm_", STAMP, ".csv")))
write_csv(otherm,    file.path(CACHE_DIR, paste0("otherm_", STAMP, ".csv")))
log_step("Saved normalized CSV tables to: ", CACHE_DIR)

# Parquet (if available)
if (have_arrow) {
  write_parquet(ppm,       file.path(CACHE_DIR, paste0("ppm_", STAMP, ".parquet")))
  write_parquet(pss_links, file.path(CACHE_DIR, paste0("pss_links_", STAMP, ".parquet")))
  write_parquet(samples,   file.path(CACHE_DIR, paste0("samples_", STAMP, ".parquet")))
  write_parquet(classm,    file.path(CACHE_DIR, paste0("classm_", STAMP, ".parquet")))
  write_parquet(effectm,   file.path(CACHE_DIR, paste0("effectm_", STAMP, ".parquet")))
  write_parquet(otherm,    file.path(CACHE_DIR, paste0("otherm_", STAMP, ".parquet")))
  log_step("Saved Parquet tables (arrow) to: ", CACHE_DIR)
}

# ---- 2) Provenance (for methods) ----
prov <- list(
  pulled_at_utc = format(Sys.time(), tz = "UTC"),
  quincunx_version = as.character(utils::packageVersion("quincunx")),
  R_version = R.version.string,
  counts = list(
    n_ppm       = nrow(ppm),
    n_pss_links = nrow(pss_links),
    n_samples   = nrow(samples),
    n_classm    = nrow(classm),
    n_effectm   = nrow(effectm),
    n_otherm    = nrow(otherm)
  )
)
if (have_jsonlite) {
  jsonlite::write_json(prov, file.path(CACHE_DIR, paste0("pm_provenance_", STAMP, ".json")),
                       pretty = TRUE, auto_unbox = TRUE)
  log_step("Wrote provenance JSON.")
} else {
  log_step("jsonlite not installed; skipping provenance JSON.")
}

