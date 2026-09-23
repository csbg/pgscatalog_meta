#!/usr/bin/env Rscript

# ============================================================
# Script: 03_bulk_metadata_pull.R
# Project: pgscatalog_meta
# Purpose: Derive training ancestry categories from the frozen Catalog download.
# Inputs: data/catalog_bulk/<pipeline_stamp>/bulk_score_development_samples_<pipeline_stamp>.csv
# Outputs: data/pgs_cache/train_bucket_<pipeline_stamp>.csv
# Run after: 02_loading_PGS_matrix.R
# Run before: 04_PGS_systemic_portability_unique_pss.R
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

.install_if_missing <- function(pkgs) {
  miss <- setdiff(pkgs, rownames(installed.packages()))
  if (length(miss)) install.packages(miss, quiet = TRUE)
}
.install_if_missing(c("readr","dplyr","stringr","tibble"))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tibble)
})

if (!file.exists("R/load.R")) {
  stop("R/load.R not found. Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

stamp <- pipeline_stamp
cache_dir <- file.path("data", "pgs_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

log_step <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

dev_fp <- catalog_bulk_file("score_development_samples", stamp)
log_step("Reading development samples: ", dev_fp)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr); library(janitor)
})

dev <- readr::read_csv(dev_fp, show_col_types = FALSE) |> janitor::clean_names()
# Expected cleaned columns include:
# polygenic_score_pgs_id, stage_of_pgs_development, broad_ancestry_category,
# ancestry_e_g_french_chinese, additional_sample_cohort_information

bucket_from_set <- function(anc_set) {
  anc_set <- unique(stats::na.omit(anc_set))
  if (length(anc_set) == 0) return(NA_character_)
  if (length(anc_set) == 1 && anc_set == "European")    return("Train: European-only")
  if (length(anc_set) == 1 && anc_set == "African")     return("Train: African-only")
  if (length(anc_set) == 1 && anc_set == "East Asian")  return("Train: East Asian-only")
  if (length(anc_set) == 1 && anc_set == "South Asian") return("Train: South Asian-only")
  if ("European" %in% anc_set) return("Train: Multi incl. EUR")
  "Train: Multi excl. EUR"
}

# Keep GWAS / training (dev) rows.
# The dev sheet may encode “training” either in stage_of_pgs_development or in additional_sample_cohort_information.
dev_train <- dev |>
  mutate(
    stage_low = tolower(stage_of_pgs_development),
    info_low  = tolower(coalesce(additional_sample_cohort_information, "")),
    is_train  = str_detect(stage_low, "gwas|source of variant|deriv|develop|train") |
      str_detect(info_low,  "train"),
    anc_raw   = coalesce(broad_ancestry_category, ancestry_e_g_french_chinese)
  ) |>
  filter(is_train, !is.na(anc_raw), !is.na(polygenic_score_pgs_id)) |>
  transmute(pgs_id = polygenic_score_pgs_id,
            ancestry_train = anc_raw)

train_bucket <- dev_train |>
  mutate(ancestry_train_disp = to_display_cat(ancestry_train)) |>
  distinct(pgs_id, ancestry_train_disp) |>
  group_by(pgs_id) |>
  summarise(train_bucket = bucket_from_set(ancestry_train_disp), .groups = "drop")

out_tb <- file.path(cache_dir, paste0("train_bucket_", stamp, ".csv"))
readr::write_csv(train_bucket, out_tb)
message("Wrote: ", out_tb, " (", nrow(train_bucket), " rows)")

# sanity checks:
train_bucket |> count(train_bucket, sort=TRUE) |> print(n=50)
dev |> count(stage_of_pgs_development, sort=TRUE) |> print(n=50)
