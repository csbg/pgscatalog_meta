#!/usr/bin/env Rscript

# ============================================================
# Script: 02_loading_PGS_matrix.R
# Project: pgscatalog_meta
# Purpose: Load cached PGS Catalog tables and validate relationships between performance metrics, sample sets, samples, and ancestry.
# Inputs: data/pgs_cache/ppm_<DATE>.csv, pss_links_<DATE>.csv, samples_<DATE>.csv, classm_<DATE>.csv
# Outputs: Console checks and validation summaries
# Run after: 01_Installation_API_and_perfomance_pull.R
# Run before: 03_bulk_metada_pull.R
# ============================================================

library(stringr)
library(purrr)
library(dplyr)
library(readr)

STAMP <- format(Sys.Date(), "%Y%m%d")
setwd("~/pgscatalog")

cache_dir <- "data/pgs_cache"

# ---- Helper: find latest COMPLETE cache stamp ----
latest_complete_stamp <- function(
    cache_dir = "data/pgs_cache",
    required_prefixes = c("ppm_", "pss_links_", "samples_", "classm_")
) {
  files <- list.files(
    cache_dir,
    pattern = "_\\d{8}\\.(csv|rds|parquet)$",
    full.names = FALSE
  )
  
  stamps <- stringr::str_extract(files, "(?<=_)\\d{8}(?=\\.)")
  stamps <- unique(stamps[!is.na(stamps)])
  
  if (!length(stamps)) {
    stop("No stamped cache files found in ", cache_dir)
  }
  
  stamp_ok <- function(stamp) {
    all(vapply(
      required_prefixes,
      function(pref) {
        any(grepl(
          paste0("^", pref, stamp, "\\.(csv|rds|parquet)$"),
          files
        ))
      },
      logical(1)
    ))
  }
  
  candidate_stamps <- sort(stamps, decreasing = TRUE)
  
  for (st in candidate_stamps) {
    if (stamp_ok(st)) return(st)
  }
  
  stop(
    "No complete cache stamp found in ", cache_dir,
    ". Need all of: ", paste(required_prefixes, collapse = ", ")
  )
}

# ---- Helper: read cache table in whatever format exists ----
read_cache_table <- function(prefix, stamp, cache_dir = "data/pgs_cache") {
  exts <- c("parquet", "rds", "csv")
  
  for (ext in exts) {
    fp <- file.path(cache_dir, paste0(prefix, stamp, ".", ext))
    if (file.exists(fp)) {
      message("Reading ", fp)
      
      if (ext == "csv") {
        return(readr::read_csv(fp, show_col_types = FALSE))
      }
      
      if (ext == "rds") {
        return(readRDS(fp))
      }
      
      if (ext == "parquet") {
        if (!requireNamespace("arrow", quietly = TRUE)) {
          stop("Package 'arrow' is required to read parquet files: ", fp)
        }
        return(arrow::read_parquet(fp))
      }
    }
  }
  
  stop(
    "No file found for prefix '", prefix,
    "' and stamp '", stamp,
    "' in ", cache_dir
  )
}

# ---- Choose the latest COMPLETE cache snapshot ----
latest_stamp <- latest_complete_stamp(
  cache_dir = cache_dir,
  required_prefixes = c("ppm_", "pss_links_", "samples_", "classm_")
)

message("Using cache stamp: ", latest_stamp)

# ---- Load near-raw tables ----
ppm       <- read_cache_table("ppm_",       latest_stamp, cache_dir)
pss_links <- read_cache_table("pss_links_", latest_stamp, cache_dir)
samples   <- read_cache_table("samples_",   latest_stamp, cache_dir)
classm    <- read_cache_table("classm_",    latest_stamp, cache_dir)

cat("\n-- ROW COUNTS --\n")
print(list(
  ppm = nrow(ppm),
  pss_links = nrow(pss_links),
  samples = nrow(samples),
  classm = nrow(classm)
))


# 3) Key sanity: ppm_id uniqueness in ppm, mapping multiplicities
stopifnot(!any(duplicated(ppm$ppm_id)))
mult_ppm_to_pss <- pss_links |> count(ppm_id) |> arrange(desc(n)) |> head()
cat("\nTop ppm_id with many pss links (expected many-to-many):\n"); print(mult_ppm_to_pss)

# 4) Stage + ancestry quick view (eval only)
samples_eval <- samples |>
  mutate(sample_size = suppressWarnings(as.numeric(sample_size)),
         ancestry_category = ancestry_display) |>
  filter(stage == "eval", !is.na(sample_size))

cat("\nEval samples by ancestry:\n")
print(samples_eval |> count(ancestry_category, sort=TRUE), n=50)

# 5) Tidy (NO top-15 filter), minimal columns, and duplication check
class_keep <- classm |>
  filter(estimate_type %in% c("AUROC","AUC","C-index")) |>
  select(ppm_id, estimate_type, estimate, interval_type, interval_lower, interval_upper)

tidy_all <- ppm |>
  select(ppm_id, pgs_id, reported_trait, covariates) |>
  left_join(pss_links, by="ppm_id") |>
  left_join(samples_eval |> 
              transmute(pss_id, sample_id,
                        ancestry_category = ancestry_category,
                        n = sample_size,
                        cases = suppressWarnings(as.numeric(sample_cases)),
                        ctrls = suppressWarnings(as.numeric(sample_controls))),
            by="pss_id") |>
  left_join(class_keep, by="ppm_id") |>
  filter(!is.na(estimate))

cat("\nTidy rows (all PGS, eval only):\n"); print(nrow(tidy_all))

# 6) Duplicate guard: there should be at most 1 row per ppm_id-pss_id-sample_id-ancestry in tidy_all
dups <- tidy_all |>
  count(ppm_id, pss_id, sample_id, ancestry_category) |>
  filter(n > 1)
if (nrow(dups)) {
  cat("\nWARNING: duplicated rows after joins (showing a few):\n")
  print(head(dups, 10))
} else {
  cat("\nNo duplicated ppm-pss-sample-ancestry rows detected. ✅\n")
}

