#!/usr/bin/env Rscript

# ============================================================
# Script: 05_pgs_auc_ci_audit.R
# Project: pgscatalog_meta
# Purpose: Harmonize AUC values and confidence intervals, audit CI availability, and produce the final evaluation table.
# Inputs: Evaluation-level PGS performance tables generated upstream
# Outputs: results/pgs_auc_ci_audit_<DATE>/eval_df_final_auc_ci.csv
# Run after: 04_PGS_systemic_portability_unique_pss.R
# Run before: 06_meta_ivw.R
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(janitor); library(purrr); library(forcats)
  library(ggplot2); library(scales); library(patchwork)
})

# =========================
# 0) CONFIG
# =========================
cache_dir   <- "data/pgs_cache"
results_dir <- "results"
stamp       <- format(Sys.Date(), "%Y%m%d")
out_dir     <- file.path(results_dir, paste0("pgs_auc_ci_audit_", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Palettes (adjust if you like)
bucket_pal <- c(
  "Train: European"    = "#3B82F6",
  "Train: African"     = "#EF4444",
  "Train: East Asian"  = "#10B981",
  "Train: South Asian" = "#F59E0B",
  "Train: Admixed"     = "#8B5CF6",
  "Train: Unknown"     = "grey60"
)
anc_pal <- c(
  "European" = "#3B82F6",
  "African"  = "#EF4444",
  "East Asian" = "#10B981",
  "South Asian"= "#F59E0B",
  "Hispanic/LatAm" = "#8B5CF6",
  "African American/Afro-Caribbean" = "#6366F1",
  "Other" = "#A3A3A3",
  "Unknown" = "grey70"
)

# =========================
# Helpers
# =========================
message_glue <- function(fmt, ...) cat(sprintf(paste0(fmt, "\n"), ...))

latest_of <- function(pattern, dir = cache_dir, recursive = FALSE, fail_if_missing = TRUE) {
  fs <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = recursive)
  if (!length(fs)) {
    if (fail_if_missing) stop(sprintf("No files matching '%s' under %s", pattern, normalizePath(dir, mustWork = FALSE)), call. = FALSE)
    return(NA_character_)
  }
  fs[order(file.info(fs)$mtime, decreasing = TRUE)][1]
}
latest_any <- function(patterns, dir = cache_dir, recursive = FALSE, fail_if_missing = TRUE) {
  for (pat in patterns) {
    fp <- latest_of(pat, dir = dir, recursive = recursive, fail_if_missing = FALSE)
    if (!is.na(fp)) return(fp)
  }
  if (fail_if_missing) stop(sprintf("No files matching any of: %s under %s", paste(patterns, collapse=" | "), normalizePath(dir, mustWork = FALSE)), call. = FALSE)
  NA_character_
}
read_csv_safe <- function(fp, ...) {
  if (is.na(fp)) return(tibble(.empty=TRUE)[,0])
  suppressMessages(readr::read_csv(fp, show_col_types = FALSE, progress = FALSE, ...))
}
coerce01 <- function(x) {
  xnum <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(xnum) & xnum > 1 & xnum <= 100, xnum/100, xnum)
}
theme_clean <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
}
canon_anc <- function(x){
  x <- as.character(x)
  dplyr::case_when(
    is.na(x) | x == "" | grepl("Unknown", x, ignore.case = TRUE) ~ "Unknown",
    x %in% c("European") ~ "European",
    x %in% c("East Asian") ~ "East Asian",
    x %in% c("South Asian") ~ "South Asian",
    x %in% c("African") ~ "African",
    x %in% c("African American or Afro-Caribbean") ~ "African American/Afro-Caribbean",
    x %in% c("Hispanic or Latin American") ~ "Hispanic/LatAm",
    grepl("Middle|Greater Middle", x, ignore.case = TRUE) ~ "Other",
    grepl("Central|South American|Oceania|Native|Admixed|unspecified|Other", x, ignore.case = TRUE) ~ "Other",
    TRUE ~ x
  )
}

# -------- Vector-safe robust parser for AUROC/C-index strings --------
parse_value_ci <- function(x) {
  parse_one <- function(x1) {
    s <- as.character(x1)
    if (is.na(s) || !nzchar(s)) return(c(NA_real_, NA_real_, NA_real_))
    s <- stringr::str_replace_all(s, "[\u2013\u2014\u2212]", "-")
    s <- stringr::str_squish(s)
    s <- stringr::str_replace_all(s, "(?i)\\b(AUROC|AUC|C-?index|C-stat(istic)?)\\s*[:=]?", "")
    s <- stringr::str_squish(s)
    
    num   <- "[-+]?\\d*[\\.,]?\\d+(?:[eE][-+]?\\d+)?"
    ci_any<- paste0("(?:\\[|\\()\\s*(", num, ")\\s*(?:,|;|\\-|to)\\s*(", num, ")\\s*(?:\\]|\\))")
    
    val1 <- stringr::str_match(s, paste0("\\b(", num, ")\\b"))[,2]
    ci_m <- stringr::str_match(s, ci_any)
    lo1  <- ci_m[,2]; hi1 <- ci_m[,3]
    if (is.na(lo1) || is.na(hi1)) {
      ci2 <- stringr::str_match(s, paste0("(?i)(?:ci|95%\\s*ci)\\s*(", num, ")\\s*(?:\\-|to)\\s*(", num, ")"))
      lo1 <- if (!is.na(ci2[,1])) ci2[,2] else lo1
      hi1 <- if (!is.na(ci2[,1])) ci2[,3] else hi1
    }
    
    norm_num <- function(z) suppressWarnings(as.numeric(gsub(",", ".", z, fixed = FALSE)))
    val_num <- norm_num(val1); lo_num <- norm_num(lo1); hi_num <- norm_num(hi1)
    
    clamp01 <- function(z) dplyr::case_when(
      is.na(z) ~ as.numeric(NA),
      z <= 1 ~ z,
      z > 1 & z <= 100 ~ z/100,
      TRUE ~ z
    )
    c(
      estimate_value    = clamp01(val_num),
      estimate_ci_lower = clamp01(lo_num),
      estimate_ci_upper = clamp01(hi_num)
    )
  }
  m <- vapply(x, parse_one, FUN.VALUE = c(estimate_value=0, estimate_ci_lower=0, estimate_ci_upper=0))
  tibble::as_tibble(t(m))
}

# =========================
# 1) Locate inputs
# =========================
message_glue(">> Scanning cache dir: %s", cache_dir)

fp_ppm <- latest_of("^bulk_performance_metrics_\\d{8}\\.csv(\\.gz)?$")
fp_pss <- latest_any(c("^bulk_evaluation_sample_sets_\\d{8}\\.csv(\\.gz)?$",
                       "^bulk_sample_sets_\\d{8}\\.csv(\\.gz)?$"))

cand_eval <- list.files(results_dir, pattern = "(^|_)eval_df(.*)\\.(csv|parquet)$",
                        full.names = TRUE, recursive = TRUE)
fp_eval_df <- if (length(cand_eval)) cand_eval[order(file.info(cand_eval)$mtime, decreasing = TRUE)][1] else NA_character_

cand_bucket <- c(
  list.files(results_dir, pattern = "bucket.*\\.csv(\\.gz)?$", full.names = TRUE, recursive = TRUE),
  list.files(cache_dir,   pattern = "bucket.*\\.csv(\\.gz)?$", full.names = TRUE, recursive = TRUE),
  list.files(results_dir, pattern = "train.*bucket.*\\.csv(\\.gz)?$", full.names = TRUE, recursive = TRUE)
)
fp_bucket <- if (length(cand_bucket)) cand_bucket[order(file.info(cand_bucket)$mtime, decreasing = TRUE)][1] else NA_character_

message_glue("  * PPM: %s", basename(fp_ppm))
message_glue("  * PSS: %s", basename(fp_pss))
message_glue("  * eval_df (optional): %s", ifelse(is.na(fp_eval_df), "<not found>", basename(fp_eval_df)))
message_glue("  * bucket map (optional): %s", ifelse(is.na(fp_bucket), "<not found>", basename(fp_bucket)))

# =========================
# 2) Read PPM & extract AUROC single + parsed CI
# =========================
ppm_raw <- read_csv_safe(fp_ppm) %>% clean_names()
col_perf <- intersect(names(ppm_raw), c("pgs_performance_metric_ppm_id","performance_id","ppm_id"))[1]
col_pgs  <- intersect(names(ppm_raw), c("evaluated_score","pgs_id"))[1]
col_pss  <- intersect(names(ppm_raw), c("pgs_sample_set_pss","sampleset_id","pss_id"))[1]
stopifnot(!is.na(col_perf), !is.na(col_pgs), !is.na(col_pss))

col_auroc <- intersect(names(ppm_raw),
                       c("area_under_the_receiver_operating_characteristic_curve_auroc",
                         "auroc","auc","area_under_the_roc_curve"))[1]
if (is.na(col_auroc)) stop("No AUROC/AUC column found in PPM CSV.")

ppm_auc_raw <- ppm_raw %>%
  transmute(
    performance_id = .data[[col_perf]],
    pgs_id         = .data[[col_pgs]],
    sampleset_id   = .data[[col_pss]],
    reported_trait = ppm_raw$reported_trait %||% NA_character_,
    auroc_raw      = .data[[col_auroc]]
  )

auc_parsed <- parse_value_ci(ppm_auc_raw$auroc_raw)

ppm_auc_norm <- bind_cols(ppm_auc_raw, auc_parsed) %>%
  transmute(
    performance_id, pgs_id, sampleset_id, reported_trait,
    auc_parsed_single   = estimate_value,
    auc_parsed_ci_lower = estimate_ci_lower,
    auc_parsed_ci_upper = estimate_ci_upper
  )

# =========================
# 3) Pull native CI from PPM/eval_df (if they ever exist)
#    and reconcile with priority: native PPM > native eval_df > parsed
# =========================
find_native_auc_cols <- function(nms) {
  # nothing numeric in current dumps, but keep hook for future
  grep("^estimate_value$|(^|_)auc_value$|(^|_)auroc_value$", nms, value = TRUE, ignore.case = TRUE)
}
find_native_ci_cols <- function(nms) {
  list(
    lo = grep("^estimate_ci_lower$|(^|_)auc_ci_lower$|(^|_)auroc_ci_lower$", nms, value = TRUE, ignore.case = TRUE),
    hi = grep("^estimate_ci_upper$|(^|_)auc_ci_upper$|(^|_)auroc_ci_upper$", nms, value = TRUE, ignore.case = TRUE)
  )
}

# 3A) native from PPM
ppm_native <- {
  auc_candidates <- find_native_auc_cols(names(ppm_raw))
  ci_candidates  <- find_native_ci_cols(names(ppm_raw))
  tibble(
    performance_id = ppm_raw[[col_perf]],
    auc_native     = if (length(auc_candidates)) coerce01(ppm_raw[[auc_candidates[1]]]) else NA_real_,
    ci_native_lower= if (length(ci_candidates$lo)) coerce01(ppm_raw[[ci_candidates$lo[1]]]) else NA_real_,
    ci_native_upper= if (length(ci_candidates$hi)) coerce01(ppm_raw[[ci_candidates$hi[1]]]) else NA_real_
  )
}

# 3B) optional native from eval_df (harmonized output you may already export)
eval_native <- if (!is.na(fp_eval_df)) {
  edf <- read_csv_safe(fp_eval_df) %>% clean_names()
  keys <- intersect(names(edf), c("performance_id","pgs_id","sampleset_id"))
  if (!("performance_id" %in% keys)) tibble(performance_id = character(), auc_native_evaldf=double(),
                                            ci_native_lower_evaldf=double(), ci_native_upper_evaldf=double())
  else edf %>%
    transmute(
      performance_id,
      auc_native_evaldf        = coerce01(edf$auc %||% NA_real_),
      ci_native_lower_evaldf   = coerce01(edf$estimate_ci_lower %||% NA_real_),
      ci_native_upper_evaldf   = coerce01(edf$estimate_ci_upper %||% NA_real_)
    )
} else tibble(performance_id = character(), auc_native_evaldf=double(),
              ci_native_lower_evaldf=double(), ci_native_upper_evaldf=double())

# 3C) Reconcile to single AUC + reported CI (priority native > parsed)
ci_records <- ppm_auc_norm %>%
  left_join(ppm_native,  by = "performance_id") %>%
  left_join(eval_native, by = "performance_id") %>%
  mutate(
    auc_native_any      = coalesce(auc_native, auc_native_evaldf),
    ci_native_lower_any = coalesce(ci_native_lower, ci_native_lower_evaldf),
    ci_native_upper_any = coalesce(ci_native_upper, ci_native_upper_evaldf),
    
    auc_single          = coalesce(auc_native_any, auc_parsed_single),
    reported_ci_lower   = coalesce(ci_native_lower_any, auc_parsed_ci_lower),
    reported_ci_upper   = coalesce(ci_native_upper_any, auc_parsed_ci_upper),
    
    # For compatibility with prior plots/steps:
    auc                 = auc_single,
    estimate_ci_lower   = reported_ci_lower,
    estimate_ci_upper   = reported_ci_upper
  )

# =========================
# 4) Join ancestry + optional train bucket
# =========================
pss <- read_csv_safe(fp_pss) %>% clean_names()
col_anc <- intersect(names(pss),
                     c("ancestry_category","ancestry_display","ancestry_prediction","ancestry_broad",
                       "broad_ancestry_category","ancestry"))
id_pss <- if ("pgs_sample_set_pss" %in% names(pss)) "pgs_sample_set_pss" else
  if ("sampleset_id" %in% names(pss)) "sampleset_id" else NA_character_

pss_pick <- if (is.na(id_pss) || !length(col_anc)) {
  warning("PSS ancestry columns not found. ancestry_eval will be NA.")
  tibble(sampleset_id = unique(ci_records$sampleset_id), ancestry_eval = NA_character_)
} else {
  pss %>%
    transmute(sampleset_id = .data[[id_pss]], ancestry_eval = .data[[col_anc[1]]]) %>%
    mutate(ancestry_eval = if_else(is.na(ancestry_eval) | ancestry_eval == "", "Unknown", ancestry_eval))
}

ci_records <- ci_records %>% left_join(pss_pick, by = "sampleset_id")

train_bucket <- if (!is.na(fp_bucket)) {
  tmp <- read_csv_safe(fp_bucket) %>% clean_names()
  id_col <- intersect(names(tmp), c("pgs_id","pgs","pgsid","evaluated_score"))[1]
  b_col  <- intersect(names(tmp), c("train_bucket","bucket","training_bucket","trainbucket"))[1]
  if (is.na(id_col) || is.na(b_col)) {
    warning("Bucket file found but could not detect 'pgs_id'/'train_bucket' columns. Skipping.")
    tibble(pgs_id=character(), train_bucket=character())
  } else tmp %>% transmute(pgs_id = .data[[id_col]], train_bucket = .data[[b_col]])
} else tibble(pgs_id=character(), train_bucket=character())

ci_records <- ci_records %>% left_join(train_bucket, by = "pgs_id")

# Export raw + augmented
write_csv(ci_records, file.path(out_dir, "ci_records_reconciled.csv"))

# =========================
# 5) Dedup baseline (ONE row per performance_id)
# =========================
ci_auc_base <- ci_records %>% distinct(performance_id, .keep_all = TRUE)
write_csv(
  ci_records %>% anti_join(ci_auc_base %>% select(performance_id), by = "performance_id"),
  file.path(out_dir, "dedup_dropped_rows.csv")
)

# =========================
# 6) Reported-only summaries & plots (baseline)
# =========================
safe_share <- function(x) mean(x, na.rm = TRUE)

summ_overall_auc <- ci_auc_base %>% summarise(
  total_rows      = n(),
  with_auc_single = sum(is.finite(auc_single)),
  with_ci_reported= sum(is.finite(estimate_ci_lower) & is.finite(estimate_ci_upper)),
  share_auc       = with_auc_single/total_rows,
  share_ci        = with_ci_reported/total_rows
)
write_csv(summ_overall_auc, file.path(out_dir, "auc_ci_summary_overall.csv"))

min_n_group <- 30; top_k_anc <- 6; width_q_cap <- 0.99

ci_auc2 <- ci_auc_base %>%
  mutate(
    ancestry_eval = canon_anc(ancestry_eval),
    train_bucket  = if_else(is.na(train_bucket) | train_bucket == "", "Train: Unknown", train_bucket)
  )

# Lollipop: share with CI (reported only)
share_by_anc <- ci_auc2 %>%
  count(ancestry_eval, has_ci = is.finite(estimate_ci_lower) & is.finite(estimate_ci_upper), name = "n") %>%
  group_by(ancestry_eval) %>%
  summarise(n_total = sum(n), n_with = sum(n[has_ci], na.rm=TRUE),
            share = n_with / n_total, .groups="drop") %>%
  filter(n_total >= min_n_group) %>%
  arrange(desc(share), desc(n_total)) %>%
  mutate(ancestry_eval = fct_inorder(ancestry_eval))

p_share_lollipop <- ggplot(share_by_anc, aes(x = share, y = ancestry_eval)) +
  geom_segment(aes(x = 0, xend = share, y = ancestry_eval, yend = ancestry_eval),
               linewidth = 0.7, color = "grey60") +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%d/%d", n_with, n_total), x = pmax(share, 0.02)),
            hjust = 0, vjust = 0.5, size = 3) +
  scale_x_continuous(labels = percent_format(accuracy=1), limits = c(0,1)) +
  labs(title = "AUC evaluations reporting CIs (reported-only, dedup baseline)",
       x = "Share with CI", y = NULL) + theme_clean()
ggsave(file.path(out_dir, "p_share_with_ci_by_ancestry_lollipop.png"), p_share_lollipop, width = 8.5, height = 6, dpi = 300)
ggsave(file.path(out_dir, "p_share_with_ci_by_ancestry_lollipop.pdf"),  p_share_lollipop, width = 8.5, height = 6)

# =========================
# 7) Hanley–McNeil inference (only when AUC present and CI missing AND we have counts)
# =========================
pss_cc <- read_csv_safe(fp_pss) %>% clean_names()
cand_cases    <- intersect(names(pss_cc), c("number_of_cases","n_cases","cases","num_cases","case_n","n_case","cases_n"))
cand_controls <- intersect(names(pss_cc), c("number_of_controls","n_controls","controls","num_controls","control_n","n_control","controls_n"))
idcol <- intersect(names(pss_cc), c("pgs_sample_set_pss","sampleset_id","pss_id"))[1]

if (is.na(idcol) || !length(cand_cases) || !length(cand_controls)) {
  warning("Case/control count columns not found in PSS. Skipping Hanley–McNeil inference.")
  hm_ci <- tibble()
  ci_auc_filled <- ci_auc_base %>%
    mutate(
      ci_lower_inferred = NA_real_,
      ci_upper_inferred = NA_real_,
      ci_lower_final    = estimate_ci_lower,
      ci_upper_final    = estimate_ci_upper,
      ci_source_final   = ifelse(is.finite(estimate_ci_lower)&is.finite(estimate_ci_upper),
                                 "reported_native_or_parsed","none"),
      has_ci_final      = !is.na(ci_lower_final) & !is.na(ci_upper_final)
    )
} else {
  pss_counts <- pss_cc %>%
    transmute(
      sampleset_id = .data[[idcol]],
      n_cases      = suppressWarnings(as.numeric(.data[[cand_cases[1]]])),
      n_controls   = suppressWarnings(as.numeric(.data[[cand_controls[1]]])),
      n_total      = suppressWarnings(as.numeric(dplyr::coalesce(.data[["number_of_individuals"]], NA_real_)))
    ) %>%
    filter(is.finite(n_cases), is.finite(n_controls), n_cases > 1, n_controls > 1) %>%
    arrange(sampleset_id) %>%
    distinct(sampleset_id, .keep_all = TRUE)
  
  auc_for_inf <- ci_auc_base %>%
    left_join(pss_counts, by = "sampleset_id") %>%
    mutate(
      has_auc   = is.finite(auc_single),
      has_ci    = is.finite(estimate_ci_lower) & is.finite(estimate_ci_upper),
      needs_ci  = has_auc & !has_ci
    ) %>%
    filter(needs_ci, is.finite(n_cases), is.finite(n_controls))
  
  hm_ci <- auc_for_inf %>%
    mutate(
      A  = auc_single,
      Q1 = A / (2 - A),
      Q2 = 2*A*A / (1 + A),
      varA = ( A*(1 - A) + (n_cases - 1)*(Q1 - A*A) + (n_controls - 1)*(Q2 - A*A) ) / (n_cases * n_controls),
      seA  = sqrt(pmax(varA, 0)),
      ci_lo_norm  = pmax(0, pmin(1, A - 1.96*seA)),
      ci_hi_norm  = pmax(0, pmin(1, A + 1.96*seA)),
      logitA      = qlogis(pmin(pmax(A, 1e-6), 1 - 1e-6)),
      se_logit    = seA / (A*(1 - A)),
      ci_lo_logit = plogis(logitA - 1.96*se_logit),
      ci_hi_logit = plogis(logitA + 1.96*se_logit)
    ) %>%
    select(performance_id, pgs_id, sampleset_id, reported_trait, ancestry_eval, train_bucket,
           auc_single, n_cases, n_controls, n_total, ci_lo_logit, ci_hi_logit)
  
  # Deduplicate per performance_id deterministically
  hm_dup_report <- hm_ci %>% count(performance_id) %>% filter(n > 1)
  if (nrow(hm_dup_report)) {
    message_glue("[INFO] Resolving %d duplicated performance_id(s) in hm_ci.", nrow(hm_dup_report))
    write_csv(hm_ci %>% semi_join(hm_dup_report, by = "performance_id") %>% arrange(performance_id),
              file.path(out_dir, "hm_ci_duplicates_report.csv"))
  }
  
  hm_ci_unique <- hm_ci %>%
    mutate(
      eff_n    = as.numeric(n_cases) + as.numeric(n_controls),
      var_proxy= (ci_hi_logit - ci_lo_logit)^2
    ) %>%
    group_by(performance_id) %>%
    arrange(desc(eff_n), var_proxy, .by_group = TRUE) %>%
    slice(1) %>% ungroup()
  
  if (nrow(hm_dup_report)) {
    kept_ids <- hm_ci_unique %>% mutate(kept = TRUE) %>%
      select(performance_id, kept, eff_n, var_proxy, n_cases, n_controls, auc_single, ci_lo_logit, ci_hi_logit)
    dropped <- hm_ci %>%
      anti_join(hm_ci_unique, by = "performance_id") %>%
      mutate(kept = FALSE,
             eff_n = as.numeric(n_cases) + as.numeric(n_controls),
             var_proxy = (ci_hi_logit - ci_lo_logit)^2) %>%
      select(performance_id, kept, eff_n, var_proxy, n_cases, n_controls, auc_single, ci_lo_logit, ci_hi_logit)
    write_csv(bind_rows(kept_ids, dropped) %>% arrange(performance_id, desc(kept)),
              file.path(out_dir, "hm_ci_dedup_decisions.csv"))
  }
  
  # Merge back to baseline
  ci_auc_filled <- ci_auc_base %>%
    left_join(hm_ci_unique %>% select(performance_id, ci_lo_logit, ci_hi_logit),
              by = "performance_id") %>%
    mutate(
      ci_lower_inferred = if_else(!is.finite(estimate_ci_lower) & is.finite(ci_lo_logit), ci_lo_logit, NA_real_),
      ci_upper_inferred = if_else(!is.finite(estimate_ci_upper) & is.finite(ci_hi_logit), ci_hi_logit, NA_real_),
      
      ci_lower_final  = coalesce(estimate_ci_lower, ci_lower_inferred),
      ci_upper_final  = coalesce(estimate_ci_upper, ci_upper_inferred),
      has_ci_final    = is.finite(ci_lower_final) & is.finite(ci_upper_final),
      ci_source_final = dplyr::case_when(
        is.finite(estimate_ci_lower) & is.finite(estimate_ci_upper) ~ "reported_native_or_parsed",
        is.finite(ci_lower_inferred) & is.finite(ci_upper_inferred) ~ "inferred_hm",
        TRUE ~ "none"
      )
    )
  stopifnot(nrow(ci_auc_filled) == nrow(ci_auc_base))
}

if (nrow(hm_ci)) write_csv(hm_ci, file.path(out_dir, "ci_inferred_hanley_mcneil_rows.csv"))
write_csv(ci_auc_filled, file.path(out_dir, "ci_auc_with_final_ci.csv"))

# =========================
# 8) Apples-to-apples lollipop (reported vs final)
# =========================
rescue_raw <- ci_auc_base %>%
  mutate(ancestry_eval = canon_anc(ancestry_eval)) %>%
  group_by(ancestry_eval) %>%
  summarise(n_total = n(),
            n_with  = sum(is.finite(estimate_ci_lower) & is.finite(estimate_ci_upper)),
            share   = n_with / n_total,
            src     = "Reported (baseline)",
            .groups = "drop")

rescue_final <- ci_auc_filled %>%
  mutate(ancestry_eval = canon_anc(ancestry_eval)) %>%
  group_by(ancestry_eval) %>%
  summarise(n_total = n(),
            n_with  = sum(has_ci_final),
            share   = n_with / n_total,
            src     = "Final (reported+inferred)",
            .groups = "drop")

rescue_by_anc <- bind_rows(rescue_raw, rescue_final) %>%
  group_by(ancestry_eval) %>%
  summarise(
    n_total        = max(n_total, na.rm = TRUE),
    share_reported = max(ifelse(src=="Reported (baseline)", share, NA_real_), na.rm = TRUE),
    share_with_ci  = max(ifelse(src=="Final (reported+inferred)", share, NA_real_), na.rm = TRUE),
    delta_share    = share_with_ci - share_reported,
    .groups="drop"
  ) %>%
  filter(n_total >= min_n_group) %>%
  arrange(desc(delta_share), desc(n_total))

write_csv(rescue_by_anc, file.path(out_dir, "ci_rescue_by_ancestry.csv"))

mk_lolli <- function(df, title){
  ggplot(df, aes(x = share, y = ancestry_eval)) +
    geom_segment(aes(x = 0, xend = share, y = ancestry_eval, yend = ancestry_eval),
                 linewidth = 0.7, color = "grey60") +
    geom_point(size = 3) +
    geom_text(aes(label = sprintf("%d/%d", n_with, n_total), x = pmax(share, 0.02)),
              hjust = 0, vjust = 0.5, size = 3) +
    scale_x_continuous(labels = percent_format(accuracy=1), limits = c(0,1)) +
    labs(title = title, x = "Share with CI", y = NULL) + theme_clean()
}
anc_levels <- rescue_by_anc %>% arrange(desc(share_with_ci)) %>% pull(ancestry_eval) %>% unique()

p_lolli_raw  <- mk_lolli(rescue_raw  %>% filter(ancestry_eval %in% anc_levels), "Reported CIs (baseline)")
p_lolli_fill <- mk_lolli(rescue_final%>% filter(ancestry_eval %in% anc_levels), "Final CIs (same baseline)")
p_side <- p_lolli_raw + p_lolli_fill + patchwork::plot_layout(ncol = 2, widths = c(1,1))
ggsave(file.path(out_dir, "p_lollipop_apples_to_apples.png"), p_side, width = 15, height = 6, dpi = 300)
ggsave(file.path(out_dir, "p_lollipop_apples_to_apples.pdf"),  p_side, width = 15, height = 6)

# =========================
# 9) AUC-aware funnel & audit
# =========================
pss_counts_any <- pss %>% clean_names() %>%
  transmute(
    sampleset_id = .data[["pgs_sample_set_pss"]] %||% .data[["sampleset_id"]],
    n_cases      = suppressWarnings(as.numeric(.data[["number_of_cases"]] %||% NA_real_)),
    n_controls   = suppressWarnings(as.numeric(.data[["number_of_controls"]] %||% NA_real_))
  ) %>%
  distinct(sampleset_id, .keep_all = TRUE)

elig <- ci_auc_base %>%
  left_join(pss_counts_any, by = "sampleset_id") %>%
  mutate(
    has_auc     = is.finite(auc_single),
    has_counts  = is.finite(n_cases) & is.finite(n_controls) & n_cases > 1 & n_controls > 1,
    has_ci_rep  = is.finite(estimate_ci_lower) & is.finite(estimate_ci_upper),
    needs_ci    = has_auc & !has_ci_rep,
    eligible    = needs_ci & has_counts
  )

funnel_auc_aware <- tibble(
  total_rows             = nrow(ci_auc_base),
  with_auc               = sum(elig$has_auc),
  reported_ci            = sum(elig$has_ci_rep),
  missing_ci_among_auc   = sum(elig$needs_ci),
  with_counts_among_auc  = sum(elig$has_counts & elig$has_auc),
  eligible_for_infer     = sum(elig$eligible)
) %>%
  mutate(
    pct_with_auc              = with_auc / total_rows,
    pct_reported_ci_over_total= reported_ci / total_rows,
    pct_missing_ci_over_auc   = ifelse(with_auc>0, missing_ci_among_auc / with_auc, NA_real_),
    pct_eligible_over_auc     = ifelse(with_auc>0, eligible_for_infer / with_auc, NA_real_)
  )
write_csv(funnel_auc_aware, file.path(out_dir, "ci_recovery_funnel_auc_aware.csv"))

reported_before <- sum(is.finite(ci_auc_base$estimate_ci_lower) & is.finite(ci_auc_base$estimate_ci_upper))
reported_now    <- sum(is.finite(ci_auc_filled$ci_lower_final)   & is.finite(ci_auc_filled$ci_upper_final))
rescued         <- reported_now - reported_before
ceiling         <- funnel_auc_aware$eligible_for_infer[1]
realized_frac   <- ifelse(ceiling > 0, rescued / ceiling, NA_real_)
tibble(total_auc = nrow(ci_auc_base),
       reported_before, reported_now, rescued, ceiling, realized_fraction = realized_frac) %>%
  write_csv(file.path(out_dir, "ci_recovery_realized_vs_ceiling.csv"))

# Audit: why eligible rows weren’t rescued (if any)
eligible_universe <- ci_auc_base %>%
  left_join(pss_counts_any, by = "sampleset_id") %>%
  mutate(
    has_auc  = is.finite(auc_single),
    has_ci   = is.finite(estimate_ci_lower) & is.finite(estimate_ci_upper),
    needs_ci = has_auc & !has_ci,
    has_counts = is.finite(n_cases) & is.finite(n_controls) & n_cases > 1 & n_controls > 1,
    eligible = needs_ci & has_counts
  ) %>% filter(eligible)

inferred_all  <- if (exists("hm_ci")) hm_ci else tibble()
inferred_kept <- if (exists("hm_ci_unique")) hm_ci_unique else tibble()

audit <- eligible_universe %>%
  transmute(performance_id, pgs_id, sampleset_id, auc_single, n_cases, n_controls) %>%
  mutate(has_auc = is.finite(auc_single)) %>%
  left_join(inferred_all  %>% select(performance_id) %>% mutate(inferred_any = TRUE), by = "performance_id") %>%
  left_join(inferred_kept %>% select(performance_id) %>% mutate(inferred_kept = TRUE), by = "performance_id") %>%
  mutate(
    inferred_any  = coalesce(inferred_any, FALSE),
    inferred_kept = coalesce(inferred_kept, FALSE),
    reason = dplyr::case_when(
      !has_auc                           ~ "AUC_missing_or_unparsed",
      has_auc & !inferred_any            ~ "not_computed_upstream",
      inferred_any & !inferred_kept      ~ "dedup_dropped_tie_or_lower_N",
      inferred_kept                      ~ "rescued",
      TRUE                               ~ "other"
    )
  )

audit_summary <- audit %>% count(reason, name = "n") %>% mutate(share = n / sum(n))
write_csv(audit,        file.path(out_dir, "ci_inference_gap_audit_rows.csv"))
write_csv(audit_summary,file.path(out_dir, "ci_inference_gap_audit_summary.csv"))

# =========================
# 10) FINAL: eval_df-like table with correct AUC & CI
#     (plus provenance columns)
# =========================
eval_df_final <- ci_auc_filled %>%
  transmute(
    performance_id,
    pgs_id,
    sampleset_id,
    reported_trait,
    ancestry_eval,
    train_bucket,
    auc                     = auc_single,             # the reconciled AUC single
    estimate_ci_lower       = ci_lower_final,         # final CI (reported or inferred)
    estimate_ci_upper       = ci_upper_final,
    ci_source               = ci_source_final,        # "reported_native_or_parsed" | "inferred_hm" | "none"
    # provenance & sanity helpers
    auc_from                = case_when(
      is.finite(auc_native) | is.finite(auc_native_evaldf) ~ "native",
      is.finite(auc_parsed_single)                          ~ "parsed_text",
      TRUE                                                  ~ "none"
    ),
    reported_ci_from        = case_when(
      is.finite(ci_native_lower_any) & is.finite(ci_native_upper_any) ~ "native",
      is.finite(auc_parsed_ci_lower) & is.finite(auc_parsed_ci_upper) ~ "parsed_text",
      TRUE                                                            ~ "none"
    )
  )

write_csv(eval_df_final, file.path(out_dir, "eval_df_final_auc_ci.csv"))

# =========================
# 11) Logging / session info
# =========================
log_lines <- c(
  sprintf("PPM: %s", fp_ppm),
  sprintf("PSS: %s", fp_pss),
  sprintf("eval_df(input): %s", ifelse(is.na(fp_eval_df), "<not found>", fp_eval_df)),
  sprintf("bucket map: %s", ifelse(is.na(fp_bucket), "<not found>", fp_bucket)),
  "",
  capture.output(sessionInfo())
)
writeLines(log_lines, file.path(out_dir, "session_info.txt"))

message_glue("==> Done. Outputs in: %s", normalizePath(out_dir, mustWork = FALSE))
