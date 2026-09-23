#!/usr/bin/env Rscript

# ============================================================
# Script: 04_PGS_systemic_portability_unique_pss.R
# Project: pgscatalog_meta
# Purpose: Build the evaluation-level portability dataset using unique evaluation sample sets.
# Inputs: data/pgs_cache/bulk_performance_metrics_<DATE>.csv, bulk_evaluation_sample_sets_<DATE>.csv, train_bucket_<DATE>.csv
# Outputs: results/pgs_systemic/ evaluation-level summary and plotting source tables
# Run after: 03_bulk_metadata_pull.R
# Run before: 05_pgs_auc_ci_audit.R
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(janitor); library(forcats); library(ggplot2); library(scales)
  library(purrr)
})

if (!file.exists("R/load.R")) {
  stop("R/load.R not found. Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

## -------------------------
## Parameters (defaults)
## -------------------------
if (!exists("min_eval_ancestries"))  min_eval_ancestries  <- 3
if (!exists("top_n_efo"))            top_n_efo            <- 10
if (!exists("top_n_pgs_barplot"))    top_n_pgs_barplot    <- 10
if (!exists("perf_metric"))          perf_metric          <- "auc"

# Gate plots via env var RUN_PLOTS=1 (default FALSE)
run_plots <- identical(Sys.getenv("RUN_PLOTS", unset = "0"), "1")

message(sprintf(
  "Params → min_eval_ancestries=%s | top_n_efo=%s | top_n_pgs_barplot=%s | metric=%s | run_plots=%s",
  min_eval_ancestries, top_n_efo, top_n_pgs_barplot, perf_metric, run_plots
))

## -------------------------
## Paths, output, stamp
## -------------------------
cache_dir <- "data/pgs_cache"
out_dir   <- file.path("results","pgs_systemic_unique_pss")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp     <- format(Sys.Date(), "%Y%m%d")

## -------------------------
## Helpers
## -------------------------
latest_of <- function(pattern, dir = cache_dir, fail_if_missing = TRUE) {
  fs <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(fs)) {
    msg <- paste0("No files matching '", pattern, "' in ", dir)
    if (fail_if_missing) stop(msg) else return(NA_character_)
  }
  head(sort(fs, decreasing = TRUE), 1)
}
calc_Neff_vec <- function(cases, ctrls) {
  cases <- as.numeric(cases); ctrls <- as.numeric(ctrls)
  out <- rep(NA_real_, length(cases))
  ok  <- !is.na(cases) & !is.na(ctrls) & cases > 0 & ctrls > 0
  out[ok] <- 4 / (1/cases[ok] + 1/ctrls[ok])
  out
}
num <- function(x) suppressWarnings(readr::parse_number(as.character(x)))

display_levels <- c("European","African","East Asian","South Asian",
                    "Hispanic or Latin American","Middle Eastern or North African",
                    "Other/Mixed","Not reported",
                    "Multi-ancestry including European","Multi-ancestry excluding European")

# Label helper (replacement for defunct label_number_si)
lab_si <- scales::label_number(accuracy = 1, scale_cut = scales::cut_short_scale())

## -------------------------
## Load newest cache files
## -------------------------
pm_fp     <- latest_of("^bulk_performance_metrics_\\d{8}\\.csv$")
pss_fp    <- latest_of("^bulk_evaluation_sample_sets_\\d{8}\\.csv$")
tb_fp     <- latest_of("^train_bucket_\\d{8}\\.csv$")
scores_fp <- latest_of("^bulk_scores_\\d{8}\\.csv$")

pm     <- read_csv(pm_fp,     show_col_types = FALSE) |> clean_names()
pss    <- read_csv(pss_fp,    show_col_types = FALSE) |> clean_names()
tb     <- read_csv(tb_fp,     show_col_types = FALSE)
scores <- read_csv(scores_fp, show_col_types = FALSE) |> clean_names()

message(
  "Using:\n  ", basename(pm_fp),
  "\n  ", basename(pss_fp),
  "\n  ", basename(tb_fp),
  "\n  ", basename(scores_fp)
)

## -------------------------
## Build eval_df (enriched) — R-tight (PSS×PPM aware)
## -------------------------
# Metrics (prefer AUC, fallback C-index)
auc_p  <- parse_est_ci(pm$area_under_the_receiver_operating_characteristic_curve_auroc)
cidx_p <- parse_est_ci(pm$concordance_statistic_c_index)

pm_class <- pm %>%
  transmute(
    pgs_id          = evaluated_score,
    ppm_id          = pgs_performance_metric_ppm_id,
    pss_id          = pgs_sample_set_pss,
    reported_trait  = reported_trait,
    metric_type_std = case_when(!is.na(auc_p$est) ~ "AUC",
                                !is.na(cidx_p$est) ~ "C-index",
                                TRUE ~ NA_character_),
    estimate        = coalesce(auc_p$est,  cidx_p$est),
    interval_lower  = coalesce(auc_p$lo,   cidx_p$lo),
    interval_upper  = coalesce(auc_p$hi,   cidx_p$hi)
  )

pss_min <- pss %>%
  transmute(
    pss_id        = pgs_sample_set_pss,
    ancestry_eval = broad_ancestry_category,
    n             = num(number_of_individuals),
    cases         = num(number_of_cases),
    ctrls         = num(number_of_controls)
  )

# PGS → EFO mapping (may be multi-EFO)
pgs_traits <- scores %>%
  transmute(
    pgs_id    = polygenic_score_pgs_id,
    efo_id    = mapped_trait_s_efo_id,
    efo_label = mapped_trait_s_efo_label
  ) %>%
  separate_rows(efo_id, efo_label, sep = ";") %>%
  mutate(
    efo_id    = str_trim(efo_id),
    efo_label = str_trim(efo_label)
  ) %>%
  distinct(pgs_id, efo_id, efo_label)

stopifnot(all(c("pgs_id","train_bucket") %in% names(tb)))
stopifnot(all(c("pgs_id","efo_id","efo_label") %in% names(pgs_traits)))

# ------- R-tight collapse logic --------
# Join & compute evaluation weights/metrics
eval_pre <- pm_class %>%
  inner_join(pss_min, by = "pss_id") %>%
  mutate(
    Neff   = calc_Neff_vec(cases, ctrls),     # evaluation-set Neff
    weight = ifelse(!is.na(Neff), Neff, n),
    auc    = estimate
  )

# Rank PPM variants within the same (PGS × PSS × ancestry_eval) and keep ONE
eval_ranked <- eval_pre %>%
  group_by(pgs_id, pss_id, ancestry_eval) %>%
  arrange(
    desc(!is.na(Neff)),   # prefer having Neff
    desc(Neff),           # higher Neff
    desc(n),              # then higher n
    desc(auc),            # then higher AUC (NA last)
    .by_group = TRUE
  ) %>%
  mutate(
    rtight_rank = dplyr::row_number(),
    kept_rtight = (rtight_rank == 1L)
  ) %>%
  ungroup()

# Log: all PPM rows dropped within the same PSS (BESS) after R-tight
bess_filtered <- eval_ranked %>%
  filter(!kept_rtight) %>%
  select(
    pgs_id, pss_id, ancestry_eval,
    ppm_id, metric_type_std, estimate, interval_lower, interval_upper,
    Neff, n, cases, ctrls, rtight_rank
  ) %>%
  arrange(pgs_id, pss_id, ancestry_eval, rtight_rank)

bess_log_fp <- file.path(out_dir, paste0("BESS_filtered_RTight_", stamp, ".csv"))
readr::write_csv(bess_filtered, bess_log_fp)

# Keep the R-tight representative rows and continue exactly as before
eval_df <- eval_ranked %>%
  filter(kept_rtight) %>%
  select(-rtight_rank, -kept_rtight) %>%
  left_join(tb,         by = "pgs_id") %>%
  left_join(pgs_traits, by = "pgs_id") %>%
  mutate(
    ancestry_raw     = as.character(ancestry_eval),
    ancestry_display = to_display_cat(ancestry_raw),
    ancestry_display = forcats::fct_relevel(ancestry_display, display_levels, after = 0)
  )

message(
  "[R-tight] Groups with >1 PPM on the same PSS: ",
  nrow(
    eval_pre %>%
      count(pgs_id, pss_id, ancestry_eval, name = "n_ppm") %>% filter(n_ppm > 1)
  ),
  " | Dropped rows: ", nrow(bess_filtered),
  " | Log: ", bess_log_fp
)

# Preview & full dump (unchanged)
write_csv(head(eval_df, 50), file.path(out_dir, paste0("eval_df_preview_", stamp, ".csv")))
write_csv(eval_df,          file.path(out_dir, paste0("eval_df_full_",    stamp, ".csv")))

## -------------------------
## TABLES (analysis sources)
## -------------------------

# Keep AUC-only for analysis tables
eval_auc <- eval_df %>% filter(metric_type_std == "AUC")

# 0) Mapping audit
tbl_ancestry_mapping <- eval_auc %>%
  count(ancestry_raw, ancestry_display, sort = TRUE)
write_csv(tbl_ancestry_mapping, file.path(out_dir, paste0("tbl_ancestry_mapping_", stamp, ".csv")))

# 1) Coverage per PGS (distinct display ancestries)
tbl_pgs_coverage <- eval_auc %>%
  group_by(pgs_id) %>%
  summarise(
    n_eval_rows          = n(),
    n_ancestries_display = n_distinct(ancestry_display),
    n_ancestries_raw     = n_distinct(ancestry_raw),
    efo_any              = paste(unique(na.omit(efo_id)), collapse = ";"),
    efo_label_any        = paste(unique(na.omit(efo_label)), collapse = ";"),
    train_bucket         = first(train_bucket),
    .groups = "drop"
  ) %>%
  arrange(desc(n_ancestries_display), desc(n_eval_rows))
write_csv(tbl_pgs_coverage, file.path(out_dir, paste0("tbl_pgs_coverage_", stamp, ".csv")))

# 2) Filter by coverage (≥ min_eval_ancestries in DISPLAY categories)
pgs_keep  <- tbl_pgs_coverage %>% filter(n_ancestries_display >= min_eval_ancestries) %>% pull(pgs_id)
eval_work <- eval_auc %>% filter(pgs_id %in% pgs_keep)

# 3A) Bar-plot source (RAW ancestries; collapsing multiple PSS rows)
tbl_bar_source_RAW <- eval_work %>%
  group_by(pgs_id, ancestry_raw, train_bucket) %>%
  summarise(
    n_eval       = sum(n, na.rm = TRUE),
    cases_total  = sum(cases, na.rm = TRUE),
    ctrls_total  = sum(ctrls, na.rm = TRUE),
    auc_mean_w   = as.numeric(weighted.mean(auc, w = ifelse(is.na(Neff), n, Neff), na.rm = TRUE)),
    auc_mean_unw = mean(auc, na.rm = TRUE),
    n_records    = n(),
    .groups = "drop"
  ) %>% arrange(pgs_id, ancestry_raw)
write_csv(tbl_bar_source_RAW, file.path(out_dir, paste0("tbl_barplot_source_RAW_", stamp, ".csv")))

# 3B) Bar-plot source (DISPLAY buckets)
tbl_bar_source_DISPLAY <- eval_work %>%
  group_by(pgs_id, ancestry_display, train_bucket) %>%
  summarise(
    n_eval       = sum(n, na.rm = TRUE),
    cases_total  = sum(cases, na.rm = TRUE),
    ctrls_total  = sum(ctrls, na.rm = TRUE),
    auc_mean_w   = as.numeric(weighted.mean(auc, w = ifelse(is.na(Neff), n, Neff), na.rm = TRUE)),
    auc_mean_unw = mean(auc, na.rm = TRUE),
    n_records    = n(),
    .groups = "drop"
  ) %>% arrange(pgs_id, ancestry_display) %>%
  mutate(n_label = lab_si(n_eval))
write_csv(tbl_bar_source_DISPLAY, file.path(out_dir, paste0("tbl_barplot_source_DISPLAY_", stamp, ".csv")))

# 4) Choose illustrative PGS (top_n_pgs_barplot by display coverage)
pgs_for_barplot <- tbl_pgs_coverage %>%
  filter(pgs_id %in% pgs_keep) %>%
  arrange(desc(n_ancestries_display), desc(n_eval_rows)) %>%
  slice_head(n = top_n_pgs_barplot) %>%
  pull(pgs_id)
write_csv(
  tibble(pgs_id = pgs_for_barplot),
  file.path(out_dir, paste0("tbl_barplot_pgs_selected_", stamp, ".csv"))
)

tbl_bar_illustrative_RAW     <- tbl_bar_source_RAW     %>% filter(pgs_id %in% pgs_for_barplot)
tbl_bar_illustrative_DISPLAY <- tbl_bar_source_DISPLAY %>% filter(pgs_id %in% pgs_for_barplot)
write_csv(tbl_bar_illustrative_RAW,     file.path(out_dir, paste0("tbl_barplot_illustrative_RAW_", stamp, ".csv")))
write_csv(tbl_bar_illustrative_DISPLAY, file.path(out_dir, paste0("tbl_barplot_illustrative_DISPLAY_", stamp, ".csv")))

# 5) Pairwise ΔAUC on DISPLAY buckets — define pairs you care about
target_pairs <- list(
  c("European","African"),
  c("East Asian","African"),
  c("South Asian","African"),
  c("East Asian","Hispanic or Latin American"),
  c("South Asian","Hispanic or Latin American"),
  c("European","Hispanic or Latin American")
)

wide_display <- tbl_bar_source_DISPLAY %>%
  select(pgs_id, ancestry_display, auc_mean_w, train_bucket) %>%
  pivot_wider(names_from = ancestry_display, values_from = auc_mean_w)

tbl_delta_pairs_DISPLAY <- map_dfr(target_pairs, function(pair) {
  a <- pair[1]; b <- pair[2]
  if (!a %in% names(wide_display) || !b %in% names(wide_display)) return(tibble())
  present <- wide_display %>% filter(!is.na(.data[[a]]) & !is.na(.data[[b]]))
  if (nrow(present) == 0) return(tibble())
  present %>%
    transmute(
      pgs_id, train_bucket,
      pair      = paste(b, "vs", a),
      delta_auc = .data[[b]] - .data[[a]],
      auc_A     = .data[[a]],
      auc_B     = .data[[b]]
    )
})
write_csv(tbl_delta_pairs_DISPLAY, file.path(out_dir, paste0("tbl_delta_pairs_DISPLAY_", stamp, ".csv")))

# 6) Top traits and ΔAUC restricted to them
tbl_trait_coverage <- eval_work %>%
  distinct(efo_id, efo_label, pgs_id) %>%
  count(efo_id, efo_label, name = "n_pgs") %>%
  arrange(desc(n_pgs))
write_csv(tbl_trait_coverage, file.path(out_dir, paste0("tbl_trait_coverage_", stamp, ".csv")))

top_efo_ids <- tbl_trait_coverage %>% slice_head(n = top_n_efo) %>% pull(efo_id)

tbl_delta_top_traits_DISPLAY <- tbl_delta_pairs_DISPLAY %>%
  inner_join(eval_work %>% distinct(pgs_id, efo_id, efo_label), by = "pgs_id") %>%
  filter(efo_id %in% top_efo_ids)
write_csv(tbl_delta_top_traits_DISPLAY, file.path(out_dir, paste0("tbl_delta_top_traits_DISPLAY_", stamp, ".csv")))

## -------------------------
## DIAGNOSTICS: summarize tables (printed + saved)
## -------------------------
diag_lines <- c(
  sprintf("PGS systemic portability — diagnostics (%s)", stamp),
  sprintf("eval_df: %s rows | %s PGS | %s EFO",
          nrow(eval_df), dplyr::n_distinct(eval_df$pgs_id), dplyr::n_distinct(eval_df$efo_id)),
  "",
  sprintf("AUC-only rows (eval_auc): %s", nrow(eval_auc)),
  sprintf("PGS with ≥%d display ancestries (pgs_keep): %s", min_eval_ancestries, length(pgs_keep)),
  "",
  "Top 10 PGS by display-ancestry coverage:",
  paste0(
    capture.output(
      tbl_pgs_coverage %>%
        select(pgs_id, n_ancestries_display, n_eval_rows, train_bucket) %>%
        slice_head(n = 10) %>% print(n = 10)
    ),
    collapse = "\n"
  ),
  "",
  "Evaluation counts per display ancestry (eval_work):",
  paste0(
    capture.output(
      eval_work %>% count(ancestry_display, sort = TRUE) %>% print(n = Inf)
    ),
    collapse = "\n"
  ),
  "",
  "Target pairs coverage (non-empty pairs only):",
  paste0(
    capture.output(
      tbl_delta_pairs_DISPLAY %>% count(pair, sort = TRUE) %>% print(n = Inf)
    ),
    collapse = "\n"
  ),
  "",
  sprintf("Selected for barplots (n=%d): %s",
          length(pgs_for_barplot), paste(pgs_for_barplot, collapse = ", "))
)
readme_fp <- file.path(out_dir, paste0("README_tables_", stamp, ".txt"))
writeLines(diag_lines, readme_fp)

message(
  "\nTables created in ", out_dir, ":\n",
  " - eval_df_full_", stamp, ".csv\n",
  " - tbl_ancestry_mapping_", stamp, ".csv\n",
  " - tbl_pgs_coverage_", stamp, ".csv\n",
  " - tbl_barplot_source_RAW_", stamp, ".csv\n",
  " - tbl_barplot_source_DISPLAY_", stamp, ".csv\n",
  " - tbl_barplot_illustrative_RAW_", stamp, ".csv\n",
  " - tbl_barplot_illustrative_DISPLAY_", stamp, ".csv\n",
  " - tbl_delta_pairs_DISPLAY_", stamp, ".csv\n",
  " - tbl_trait_coverage_", stamp, ".csv\n",
  " - tbl_delta_top_traits_DISPLAY_", stamp, ".csv\n",
  " - README_tables_", stamp, ".txt\n"
)


## =========================
## PLOTS — optional diagnostic portability plots
## expects the TABLES section above has already run
## =========================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(forcats); library(ggplot2); library(scales); library(purrr); library(readr)
})

# Helpers (same as in tables)
if (!exists("display_levels")) display_levels <- c(
  "European","African","East Asian","South Asian",
  "Hispanic or Latin American","Middle Eastern or North African",
  "Other/Mixed","Not reported",
  "Multi-ancestry including European","Multi-ancestry excluding European"
)
lab_si <- scales::label_number(accuracy = 1, scale_cut = scales::cut_short_scale())

# Load stamped tables generated above

load_csv <- function(fp) readr::read_csv(fp, show_col_types = FALSE)

tbl_pgs_cov_fp       <- file.path(out_dir, paste0("tbl_pgs_coverage_",               stamp, ".csv"))
tbl_src_disp_fp      <- file.path(out_dir, paste0("tbl_barplot_source_DISPLAY_",     stamp, ".csv"))
eval_df_full_fp      <- file.path(out_dir, paste0("eval_df_full_",                   stamp, ".csv"))
tbl_trait_cov_fp     <- file.path(out_dir, paste0("tbl_trait_coverage_",             stamp, ".csv"))

if (!exists("tbl_pgs_coverage"))       tbl_pgs_coverage       <- load_csv(tbl_pgs_cov_fp)
if (!exists("tbl_bar_source_DISPLAY")) tbl_bar_source_DISPLAY <- load_csv(tbl_src_disp_fp)
if (!exists("eval_df"))                eval_df                <- load_csv(eval_df_full_fp)
if (!exists("tbl_trait_coverage"))     tbl_trait_coverage     <- load_csv(tbl_trait_cov_fp)

## ---- Visual defaults (source once) ----
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

# Colorblind-safe Okabe–Ito
pal_okabe <- c(
  "Train: European-only"          = "#0072B2",
  "Train: Multi incl. EUR"        = "#009E73",
  "Train: Multi excl. EUR"        = "#E69F00",
  "Train: African-only"           = "#D55E00",
  "NA"                            = "#666666"
)

# Shapes mirror the same mapping (helps when printed in grayscale)
shape_map <- c(
  "Train: European-only"   = 16,  # filled circle
  "Train: Multi incl. EUR" = 17,  # filled triangle
  "Train: Multi excl. EUR" = 15,  # filled square
  "Train: African-only"    = 18,  # filled diamond
  "NA"                     = 1
)

# Unified minimal theme tuned for dense labels
theme_portability <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 9, angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title.position = "plot"
  )

# Helper for compact N labels
lab_si <- scales::label_number(accuracy = 1, scale_cut = scales::cut_short_scale())

## ---------------------------------------------------------------
## Step 2 — Single-PGS illustration (bar plots), top 10 by breadth + spread
## ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(stringr)
})

# ---- Visual defaults ----
pal_okabe_vec <- c(
  "#0072B2", "#009E73", "#E69F00", "#D55E00",
           "#CC79A7", "#F0E442", "#56B4E9", "#000000", "#999999"
)

theme_portability <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90"),
    strip.text         = element_text(face = "bold"),
    axis.text.x        = element_text(size = 9, angle = 25, hjust = 1),
    legend.position    = "bottom",
    legend.box         = "vertical",
    plot.title.position= "plot"
  )

lab_si <- scales::label_number(accuracy = 1, scale_cut = scales::cut_short_scale())

# ---- Load stamped tables (from your TABLES section) ----
tbl_pgs_cov_fp  <- file.path(out_dir, paste0("tbl_pgs_coverage_",           stamp, ".csv"))
tbl_src_disp_fp <- file.path(out_dir, paste0("tbl_barplot_source_DISPLAY_", stamp, ".csv"))

tbl_pgs_coverage       <- readr::read_csv(tbl_pgs_cov_fp,  show_col_types = FALSE)
tbl_bar_source_DISPLAY <- readr::read_csv(tbl_src_disp_fp, show_col_types = FALSE)

# ---- Display order fallback ----
if (!exists("display_levels")) {
  display_levels <- c(
    "European", "African", "East Asian", "South Asian",
    "Hispanic or Latin American", "Middle Eastern or North African",
    "Other/Mixed", "Not reported",
    "Multi-ancestry including European", "Multi-ancestry excluding European"
  )
}

# ---- Parameters fallback ----
if (!exists("min_eval_ancestries")) min_eval_ancestries <- 3
if (!exists("top_n_pgs_barplot"))   top_n_pgs_barplot   <- 10

# ---- Select top PGS (breadth + spread) ----
pgs_eligible <- tbl_pgs_coverage %>%
  filter(n_ancestries_display >= min_eval_ancestries) %>%
  select(pgs_id, n_eval_rows, train_bucket)

pgs_spread <- tbl_bar_source_DISPLAY %>%
  group_by(pgs_id) %>%
  summarise(
    n_ancestries_display = n_distinct(ancestry_display),
    auc_spread           = diff(range(auc_mean_w, na.rm = TRUE)),
    .groups              = "drop"
  )

pgs_for_barplot <- pgs_eligible %>%
  inner_join(pgs_spread, by = "pgs_id") %>%
  arrange(desc(n_ancestries_display), desc(auc_spread), desc(n_eval_rows)) %>%
  slice_head(n = top_n_pgs_barplot) %>%
  pull(pgs_id)

# ---- Build phenotype (EFO) captions per PGS ----
# Use the concatenated efo_label_any column; show up to 2 labels, then "+N more"
pgs_meta <- tbl_pgs_coverage %>%
  filter(pgs_id %in% pgs_for_barplot) %>%
  transmute(
    pgs_id,
    efo_label_any = coalesce(efo_label_any, "Unknown phenotype")
  ) %>%
  rowwise() %>%
  mutate(
    efo_vec   = {
      v <- str_split(efo_label_any, ";", simplify = FALSE)[[1]]
      v <- unique(trimws(v))
      v[nzchar(v)]
    },
    efo_caption = {
      shown <- paste(head(efo_vec, 2), collapse = ", ")
      if (length(efo_vec) > 2) paste0(shown, " +", length(efo_vec) - 2, " more") else shown
    }
  ) %>%
  ungroup() %>%
  mutate(
    # Final facet label: "PGS000123 — Breast carcinoma, XYZ (+more)"
    facet_label = paste0(pgs_id, " — ", efo_caption)
  ) %>%
  select(pgs_id, facet_label)

# Keep facet order aligned with pgs_for_barplot
facet_levels <- pgs_meta %>%
  right_join(tibble(pgs_id = pgs_for_barplot), by = "pgs_id") %>%
  pull(facet_label)

# ---- Plot data ----
barplot_data <- tbl_bar_source_DISPLAY %>%
  filter(pgs_id %in% pgs_for_barplot) %>%
  left_join(pgs_meta, by = "pgs_id") %>%
  mutate(
    ancestry_display = fct_relevel(ancestry_display, display_levels),
    n_label          = lab_si(n_eval),
    facet_label      = factor(facet_label, levels = facet_levels)
  ) %>%
  arrange(facet_label, ancestry_display)

if (nrow(barplot_data) == 0) stop("No rows for barplot after selection.")

# ---- Legend control (force show, stable order) ----
levels_train <- levels(fct_infreq(barplot_data$train_bucket))
fill_vals    <- setNames(pal_okabe_vec[seq_len(length(levels_train))], levels_train)

# ---- Plot ----
p_barplots <- ggplot(
  barplot_data,
  aes(x = ancestry_display, y = auc_mean_w, fill = train_bucket)
) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey75") +
  geom_col(width = 0.62, color = "grey25", linewidth = 0.15) +
  geom_text(
    aes(label = n_label, y = pmin(0.995, auc_mean_w + 0.012)),
    color = "black", vjust = 0, size = 3.0, fontface = "bold", show.legend = FALSE
  ) +
  facet_wrap(~ facet_label, scales = "free_x") +
  coord_cartesian(ylim = c(0.5, 1.0)) +
  scale_y_continuous(breaks = seq(0.5, 1.0, 0.05)) +
  scale_fill_manual(
    name   = "Training bucket",
    breaks = levels_train,
    values = fill_vals,
    drop   = FALSE
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      nrow = 2,
      byrow = TRUE,
      override.aes = list(alpha = 1)
    )
  ) +
  labs(
    title    = "AUC by evaluation ancestry",
    subtitle = "Facet label shows PGS ID and phenotype(s); numbers above bars = evaluation N (short scale)",
    x = "Evaluation ancestry",
    y = "AUC"
  ) +
  theme_portability +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold")
  )

# ---- Save ----
outfile_pdf <- file.path(out_dir, paste0("barplots_topPGS_", stamp, "_with_phenotype.pdf"))
ggsave(outfile_pdf, p_barplots, width = 12, height = 11, device = cairo_pdf)
message("Saved bar plots → ", outfile_pdf)
