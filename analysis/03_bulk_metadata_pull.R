#!/usr/bin/env Rscript

# ============================================================
# Script: 03_bulk_metada_pull.R
# Project: pgscatalog_meta
# Purpose: Download PGS Catalog bulk metadata and derive training ancestry categories.
# Inputs: PGS Catalog bulk metadata files from official FTP
# Outputs: data/pgs_cache/bulk_*_<DATE>.csv; data/pgs_cache/train_bucket_<DATE>.csv
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

STAMP     <- format(Sys.Date(), "%Y%m%d")
CACHE_DIR <- file.path("data","pgs_cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

log_step <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# --- Bulk metadata file URLs (official FTP) ---
base <- "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata"
files <- c(
  "pgs_all_metadata_scores.csv",
  "pgs_all_metadata_score_development_samples.csv",
  "pgs_all_metadata_performance_metrics.csv",
  "pgs_all_metadata_evaluation_sample_sets.csv",
  "pgs_all_metadata_cohorts.csv",
  "pgs_all_metadata_publications.csv",
  "pgs_all_metadata_efo_traits.csv"
)

dest <- file.path(CACHE_DIR, sub("^pgs_all_metadata_", "bulk_", files))
dest <- sub("\\.csv$", paste0("_", STAMP, ".csv"), dest)

# --- Download if missing (idempotent) ---
for (i in seq_along(files)) {
  src <- paste0(base, "/", files[i])
  if (!file.exists(dest[i])) {
    log_step("Downloading: ", files[i])
    utils::download.file(src, dest[i], mode = "wb", quiet = TRUE)
  } else {
    log_step("Exists, skipping: ", basename(dest[i]))
  }
}

# --- Load the two we need right now for training ancestry ---
scores_fp <- dest[basename(dest) == paste0("bulk_scores_", STAMP, ".csv")]
dev_fp    <- dest[basename(dest) == paste0("bulk_score_development_samples_", STAMP, ".csv")]

scores <- readr::read_csv(scores_fp, show_col_types = FALSE)
dev    <- readr::read_csv(dev_fp,    show_col_types = FALSE)
# ---- Build train_bucket from the bulk score-development sample file ----

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr); library(janitor)
})

cache_dir <- "data/pgs_cache"
stamp <- format(Sys.Date(), "%Y%m%d")

cache_dir <- "data/pgs_cache"
dev_fp <- list.files(cache_dir, pattern="^bulk_score_development_samples_\\d{8}\\.csv$", full.names=TRUE) |>
  sort(decreasing = TRUE) |>
  head(1)
stopifnot(length(dev_fp) == 1, file.exists(dev_fp))

dev <- readr::read_csv(dev_fp, show_col_types = FALSE) |> janitor::clean_names()
# Expected cleaned columns include:
# polygenic_score_pgs_id, stage_of_pgs_development, broad_ancestry_category,
# ancestry_e_g_french_chinese, additional_sample_cohort_information

# Map ancestry to your display buckets
to_display_cat <- function(x) {
  x <- gsub("\\s*,\\s*", ",", x)
  sapply(strsplit(x, ","), function(v) {
    v <- unique(trimws(v))
    if (length(v) > 1) {
      if ("European" %in% v) "Multi-ancestry including European" else "Multi-ancestry excluding European"
    } else {
      vv <- v[1]
      dplyr::case_when(
        vv %in% c("European","African","East Asian","South Asian",
                  "Hispanic or Latin American","Middle Eastern or North African",
                  "Other/Mixed","Not reported",
                  "Multi-ancestry including European","Multi-ancestry excluding European") ~ vv,
        vv %in% c("African American or Afro-Caribbean","African unspecified","Sub-Saharan African") ~ "African",
        vv == "Hispanic or Latin American" ~ "Hispanic or Latin American",
        vv == "Greater Middle Eastern (Middle Eastern, North African or Persian)" ~ "Middle Eastern or North African",
        vv %in% c("Central Asian","South East Asian","Asian unspecified","Oceanian","Native American",
                  "Aboriginal Australian","Other","Other admixed ancestry") ~ "Other/Mixed",
        vv == "Not reported" ~ "Not reported",
        TRUE ~ "Other/Mixed"
      )
    }
  }, USE.NAMES = FALSE)
}

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
