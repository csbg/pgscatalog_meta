library(stringr)

cache_dir <- "data/pgs_cache"

# Look for ANY of the cache tables and extract their YYYYMMDD stamps
files  <- list.files(cache_dir, pattern = "_\\d{8}\\.csv$|_\\d{8}\\.rds$|_\\d{8}\\.parquet$", full.names = FALSE)
stamps <- str_extract(files, "(?<=_)\\d{8}(?=\\.)")
stamps <- stamps[!is.na(stamps)]
stopifnot(length(stamps) > 0)

latest_stamp <- sort(unique(stamps), decreasing = TRUE)[1]
latest_stamp

# 2) Load near-raw tables
library(readr)

ppm       <- read_csv(file.path(cache_dir, paste0("ppm_",       latest_stamp, ".csv")), show_col_types = FALSE)
pss_links <- read_csv(file.path(cache_dir, paste0("pss_links_", latest_stamp, ".csv")), show_col_types = FALSE)
samples   <- read_csv(file.path(cache_dir, paste0("samples_",   latest_stamp, ".csv")), show_col_types = FALSE)
classm    <- read_csv(file.path(cache_dir, paste0("classm_",    latest_stamp, ".csv")), show_col_types = FALSE)

cat("\n-- ROW COUNTS --\n")
print(list(ppm=nrow(ppm), pss_links=nrow(pss_links), samples=nrow(samples), classm=nrow(classm)))

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

