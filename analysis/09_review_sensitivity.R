#!/usr/bin/env Rscript

# ============================================================
# Script: 09_review_sensitivity.R
# Project: pgscatalog_meta
# Purpose: Review-sensitivity analyses on the same Catalog stamp as Figure 1C.
# Inputs: frozen bulk metadata, the Stage-1 tables, and the script 08 forest.
# Outputs: results/review_sensitivity_<STAMP>/tables and plots
# Run after: 08_generate_publication_figures.R
# ============================================================
# Rebuilds the recovery-pack figures from this download and the single
# assigned evaluation-ancestry rule. Titles and captions do not name
# referees. The publication Figure 1C remains the script 08 forest.
# Catalog evaluations whose cohort is All of Us are flagged; they are
# not the manuscript's own All of Us analysis.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(glue)
  library(stringr)
  library(forcats)
  library(scales)
  library(purrr)
  library(janitor)
  library(patchwork)
})

if (!file.exists("R/load.R")) {
  stop("R/load.R not found. Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

stamp <- pipeline_stamp
root_out   <- Sys.getenv("REVIEW_CHECK_OUT", unset = file.path("results", paste0("review_sensitivity_", stamp)))
dir_tables <- file.path(root_out, "tables")
dir_plots  <- file.path(root_out, "plots")
invisible(lapply(c(dir_tables, dir_plots), dir.create, recursive = TRUE, showWarnings = FALSE))

pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf

save_plot <- function(filename, plot, width, height, ...) {
  ggsave(
    file.path(dir_plots, paste0(filename, ".pdf")),
    plot, width = width, height = height, device = pdf_device, ...
  )
  ggsave(
    file.path(dir_plots, paste0(filename, ".png")),
    plot, width = width, height = height, dpi = 300, bg = "white", ...
  )
}

write_tbl <- function(x, name) write_csv(x, file.path(dir_tables, paste0(name, ".csv")))

# ---------------------------------------------------------------------------
# Shared Stage-1 helpers come from R/. Stage-2 rebuilds below keep the
# extra SE scales and the random-effects model used only in this script.

pool_stage1 <- function(df, by_class = TRUE) {
  if (isTRUE(by_class) && !"covariate_class" %in% names(df)) {
    stop("Stage 1 needs covariate_class.", call. = FALSE)
  }
  groups <- c("pgs_id", "trait_label", "ancestry_display")
  if (isTRUE(by_class)) groups <- c(groups, "covariate_class")
  df %>%
    mutate(tmp = pmap(list(auc, estimate_ci_lower, estimate_ci_upper), auc_ci_to_logit)) %>%
    tidyr::unnest_wider(tmp) %>%
    filter(is.finite(eta), is.finite(se), se > 0) %>%
    group_by(across(all_of(groups))) %>%
    summarise(pooled = list(ivw_pool_logit(eta, se)), .groups = "drop") %>%
    tidyr::unnest_wider(pooled) %>%
    mutate(
      auc_pooled = plogis(eta),
      lo_pooled  = plogis(eta - 1.96 * se),
      hi_pooled  = plogis(eta + 1.96 * se)
    )
}

i2_filter <- function(stage1_tbl) {
  stage1_tbl %>%
    mutate(flag_high_I2 = !is.na(I2) & k_eval >= 2 & I2 > 80) %>%
    filter(!flag_high_I2) %>%
    select(-flag_high_I2)
}

attach_train_bucket <- function(stage1_tbl, eval_df) {
  bucket_map <- eval_df %>%
    filter(!is.na(pgs_id), !is.na(train_bucket), train_bucket != "") %>%
    distinct(pgs_id, train_bucket) %>%
    group_by(pgs_id) %>%
    slice(1) %>%
    ungroup() %>%
    rename(trained_bucket = train_bucket)
  stage1_tbl %>%
    mutate(pgs_id = as.character(pgs_id)) %>%
    left_join(bucket_map, by = "pgs_id")
}

keep_bucket_pattern <- c("European-only", "Multi incl")

# ---------------------------------------------------------------------------
# Stage 2 rebuild, parameterised by SE scale and pooling model
#   se_scale:
#     "logit_as_implemented"  submitted Stage-2 standard error (logit scale)
#     "delta_method"          se_auc = se_logit * AUC(1-AUC)          (proposed fix)
#     "from_ci"               se_auc = (hi_pooled - lo_pooled)/(2*1.96) (cross-check)
#   model:
#     "fixed"                 inverse-variance fixed effect (as submitted)
#     "random_dl"             DerSimonian-Laird random effects
# ---------------------------------------------------------------------------

se_scale_levels <- c("logit_as_implemented", "delta_method", "from_ci")

stage1_to_scatter <- function(stage1_tbl, keep_traits) {
  if (!"covariate_class" %in% names(stage1_tbl)) {
    stage1_tbl$covariate_class <- NA_character_
  }
  stage1_tbl %>%
    transmute(
      trait_label,
      ancestry_display = as.character(ancestry_display),
      pgs_id,
      auc_pooled, lo_pooled, hi_pooled,
      se_logit = se, k_eval, trained_bucket,
      covariate_class = as.character(covariate_class)
    ) %>%
    filter(
      ancestry_display != "Not reported",
      grepl(paste(keep_bucket_pattern, collapse = "|"), trained_bucket, ignore.case = TRUE),
      trait_label %in% keep_traits
    ) %>%
    mutate(
      se_delta_method = se_logit * auc_pooled * (1 - auc_pooled),
      se_from_ci      = (hi_pooled - lo_pooled) / (2 * 1.96)
    )
}

pair_stage1 <- function(scat_df, se_scale = "logit_as_implemented") {
  se_scale <- match.arg(se_scale, se_scale_levels)
  se_col <- c(logit_as_implemented = "se_logit", delta_method = "se_delta_method", from_ci = "se_from_ci")[[se_scale]]
  if (!"covariate_class" %in% names(scat_df) || anyNA(scat_df$covariate_class)) {
    stop("Stage-1 rows are missing covariate_class.", call. = FALSE)
  }
  pair_stage1_by_class(scat_df, se_col = se_col) %>%
    pick_one_class_per_score() %>%
    mutate(
      delta    = tgt_auc - eur_auc,
      delta_se = sqrt(tgt_se^2 + eur_se^2),
      se_scale = se_scale
    ) %>%
    filter(is.finite(delta), is.finite(delta_se), delta_se > 0)
}

# The pairing used before covariate class was part of the join. Kept only
# to record which deltas change under the class rule.
pair_stage1_ignoring_class <- function(scat_df, se_col) {
  eur_df <- scat_df %>%
    filter(ancestry_display == "European") %>%
    transmute(
      trait_label, pgs_id, trained_bucket,
      eur_auc = auc_pooled, eur_se = .data[[se_col]], eur_k = k_eval
    )
  tgt_df <- scat_df %>%
    filter(ancestry_display != "European") %>%
    transmute(
      trait_label, pgs_id, trained_bucket,
      target_ancestry = ancestry_display,
      tgt_auc = auc_pooled, tgt_se = .data[[se_col]], tgt_k = k_eval
    )
  tgt_df %>%
    inner_join(eur_df, by = c("trait_label", "pgs_id", "trained_bucket")) %>%
    mutate(
      delta = tgt_auc - eur_auc,
      delta_se = sqrt(tgt_se^2 + eur_se^2)
    ) %>%
    filter(is.finite(delta), is.finite(delta_se), delta_se > 0)
}

pool_delta <- function(delta, se, model = "fixed") {
  ok <- is.finite(delta) & is.finite(se) & se > 0
  delta <- delta[ok]
  se    <- se[ok]
  k <- length(delta)
  if (k == 0) {
    return(tibble(delta_hat = NA_real_, se_hat = NA_real_, tau2 = NA_real_, Q_stage2 = NA_real_, I2_stage2 = NA_real_))
  }
  w  <- 1 / se^2
  dh <- sum(w * delta) / sum(w)
  Q  <- sum(w * (delta - dh)^2)
  df <- k - 1
  tau2 <- 0
  if (model == "random_dl" && k >= 2) {
    C <- sum(w) - sum(w^2) / sum(w)
    tau2 <- if (is.finite(C) && C > 0) max(0, (Q - df) / C) else 0
  }
  ws <- 1 / (se^2 + tau2)
  tibble(
    delta_hat = sum(ws * delta) / sum(ws),
    se_hat    = sqrt(1 / sum(ws)),
    tau2      = tau2,
    Q_stage2  = if (k >= 2) Q else NA_real_,
    I2_stage2 = if (k >= 2 && Q > 0) max(0, (Q - df) / Q) * 100 else NA_real_
  )
}

pool_stage2 <- function(paired_df, model = "fixed") {
  paired_df %>%
    group_by(trait_label, trained_bucket, target_ancestry) %>%
    summarise(
      n_pairs        = n(),
      all_not_pooled = all(!is_pooled, na.rm = TRUE),
      k_eval_eur_min = min(eur_k, na.rm = TRUE),
      k_eval_tgt_min = min(tgt_k, na.rm = TRUE),
      res            = list(pool_delta(delta, delta_se, model)),
      .groups        = "drop"
    ) %>%
    tidyr::unnest_wider(res) %>%
    mutate(
      lo = delta_hat - 1.96 * se_hat,
      hi = delta_hat + 1.96 * se_hat,
      ci_excludes_zero = (lo > 0) | (hi < 0),
      model = model
    )
}

rebuild_paired_delta <- function(stage1_tbl, keep_traits, se_scale = "logit_as_implemented", model = "fixed") {
  paired <- pair_stage1(stage1_to_scatter(stage1_tbl, keep_traits), se_scale = se_scale)
  list(paired = paired, cells = pool_stage2(paired, model = model))
}

# ---------------------------------------------------------------------------
# Plot conventions
# ---------------------------------------------------------------------------

theme_reviewer <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

# Colour-safe (Okabe-Ito) palette: editor's colour-vision-deficiency request
pal_ancestry <- c(
  "European"                          = "#0072B2",
  "African"                           = "#D55E00",
  "East Asian"                        = "#56B4E9",
  "South Asian"                       = "#009E73",
  "Hispanic or Latin American"        = "#E69F00",
  "Middle Eastern or North African"   = "#CC79A7",
  "Other/Mixed"                       = "#999999",
  "Not reported"                      = "#000000",
  "Multi-ancestry including European" = "#F0E442",
  "Multi-ancestry excluding European" = "#8E6718"
)

display_levels <- names(pal_ancestry)

anc_keep <- c("European", "African", "East Asian", "South Asian",
              "Hispanic or Latin American", "Middle Eastern or North African")

trait_display <- c(
  "Incident chronic obstructive pulmonary disease before age 50 years" = "Incident COPD before age 50",
  "Coronary heart disease (incident and prevalent)" = "Coronary heart disease (incident + prevalent)",
  "Incident coronary heart disease" = "Incident coronary heart disease",
  "Stroke excluding subarachnoid hemorrhage" = "Stroke (excl. subarachnoid hemorrhage)",
  "Epithelial non-mucinous ovarian cancer" = "Ovarian cancer (epithelial, non-mucinous)",
  "Apparent Treatment-Resistant Hypertension" = "Treatment-resistant hypertension"
)
wrap_trait <- function(x) {
  x <- as.character(x)
  x <- ifelse(x %in% names(trait_display), trait_display[x], x)
  str_wrap(x, 24)
}

fmt_cases <- function(x) {
  ifelse(
    !is.finite(x), NA_character_,
    ifelse(x >= 1000, paste0(formatC(x / 1000, format = "f", digits = 1), "k"),
           formatC(round(x), format = "d", big.mark = ","))
  )
}

pal_bucket <- c("Train: European-only" = "#0072B2", "Train: Multi incl. EUR" = "#D55E00")
pooled_levels <- c("Pooled Stage-1 input", "Single evaluation per side")
dodge_w <- 0.65
x_lim_forest <- c(-0.35, 0.30)
x_lab_pos <- 0.31
x_breaks_forest <- seq(-0.3, 0.3, by = 0.1)
x_labels_forest <- label_number(accuracy = 0.1, style_negative = "minus")

# Generic two-or-more-version forest (used by Sections 6, 8, 9, 10)
forest_compare_plot <- function(long_df, pal, shapes, title, subtitle, caption,
                                lost_df = NULL, lost_label = "not estimable",
                                trait_levels, x_lim = x_lim_forest, x_lab = x_lab_pos,
                                x_breaks = x_breaks_forest, x_accuracy = 0.1,
                                label_mode = c("all", "first_k_others_star"),
                                show_k = TRUE) {
  label_mode <- match.arg(label_mode)
  x_labels <- label_number(accuracy = x_accuracy, style_negative = "minus")
  long_df <- long_df %>%
    mutate(
      trait_wrapped = factor(wrap_trait(trait_label), levels = trait_levels),
      target_ancestry = factor(as.character(target_ancestry), levels = rev(anc_keep)),
      trained_bucket = factor(trained_bucket, levels = names(pal_bucket)),
      sig = (lo > 0) | (hi < 0),
      k_lab = paste0("k = ", n_pairs, ifelse(sig, " *", ""))
    )
  if (label_mode == "first_k_others_star") {
    first_src <- levels(long_df$source)[1]
    long_df <- long_df %>%
      mutate(k_lab = ifelse(source == first_src, k_lab, ifelse(sig, "*", "")))
  }
  p <- ggplot(
    long_df,
    aes(x = delta_hat, y = target_ancestry, colour = source, shape = source, group = source)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_errorbar(
      aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.5,
      position = position_dodge(width = dodge_w), orientation = "y"
    ) +
    geom_point(size = 2.3, position = position_dodge(width = dodge_w))
  if (show_k) {
    p <- p + geom_text(
      aes(x = x_lab, label = k_lab, fontface = ifelse(sig, "bold", "plain")),
      size = 2.5, hjust = 0, colour = "grey20",
      position = position_dodge(width = dodge_w), show.legend = FALSE
    )
  }
  if (!is.null(lost_df) && nrow(lost_df) > 0) {
    lost_df <- lost_df %>%
      mutate(
        trait_wrapped = factor(wrap_trait(trait_label), levels = trait_levels),
        target_ancestry = factor(as.character(target_ancestry), levels = rev(anc_keep)),
        trained_bucket = factor(trained_bucket, levels = names(pal_bucket))
      )
    p <- p + geom_text(
      data = lost_df,
      aes(x = x_lab, y = target_ancestry, label = lost_label),
      inherit.aes = FALSE, size = 2.5, hjust = 0, colour = "#E15759", fontface = "italic",
      nudge_y = 0.33
    )
  }
  p +
    scale_colour_manual(values = pal, name = NULL) +
    scale_shape_manual(values = shapes, name = NULL) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    coord_cartesian(xlim = x_lim, clip = "off") +
    facet_grid(trait_wrapped ~ trained_bucket, scales = "free_y", space = "free_y", switch = "y") +
    labs(title = title, subtitle = subtitle, x = "ΔAUC (target − European) with 95% CI", y = NULL, caption = caption) +
    theme_reviewer(base_size = 9) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8, face = "bold", lineheight = 0.85),
      strip.text.x = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.spacing.y = unit(0.35, "lines"),
      panel.spacing.x = unit(1.6, "lines"),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(5.5, 45, 5.5, 5.5),
      plot.caption = element_text(hjust = 0, size = 7.5),
      plot.title = element_text(size = 10.5)
    )
}

# ---------------------------------------------------------------------------
# 0) Load snapshot files
# ---------------------------------------------------------------------------

message(glue(">> Loading Catalog stamp {stamp}"))

fp_eval   <- file.path("results", paste0("pgs_auc_ci_audit_", stamp), "eval_df_final_auc_ci.csv")
fp_funnel <- file.path("results", paste0("pgs_auc_ci_audit_", stamp), "ci_recovery_funnel_auc_aware.csv")
fp_pss    <- catalog_bulk_file("evaluation_sample_sets", stamp)
fp_ppm    <- catalog_bulk_file("performance_metrics", stamp)
fp_scores <- catalog_bulk_file("scores", stamp)
fp_s1     <- file.path(catalog_results_dir(stamp), "stage1", paste0("stage1_pooled_cells_", stamp, ".csv"))
fp_s1f    <- file.path(catalog_results_dir(stamp), "stage1", paste0("stage1_pooled_cells_I2filtered_", stamp, ".csv"))
fp_forest <- file.path(catalog_results_dir(stamp), "stage2", paste0("deltaAUC_forest_faceted_data_", stamp, ".csv"))

stopifnot(
  file.exists(fp_eval), file.exists(fp_pss), file.exists(fp_ppm), file.exists(fp_scores),
  file.exists(fp_s1), file.exists(fp_s1f), file.exists(fp_forest)
)

eval_raw   <- read_csv(fp_eval, show_col_types = FALSE)
funnel_in  <- if (file.exists(fp_funnel)) read_csv(fp_funnel, show_col_types = FALSE) else NULL
stage1     <- read_csv(fp_s1, show_col_types = FALSE)
stage1_f   <- read_csv(fp_s1f, show_col_types = FALSE)
forest_raw <- read_csv(fp_forest, show_col_types = FALSE)

pss_raw <- read_csv(fp_pss, show_col_types = FALSE, col_types = cols(.default = col_character()))
pss_field_names <- names(pss_raw)
pss <- pss_raw %>%
  janitor::clean_names() %>%
  transmute(
    sampleset_id  = as.character(pgs_sample_set_pss),
    n_individuals = suppressWarnings(as.numeric(number_of_individuals)),
    n_cases       = suppressWarnings(as.numeric(number_of_cases)),
    n_controls    = suppressWarnings(as.numeric(number_of_controls)),
    pct_male      = suppressWarnings(as.numeric(percent_of_participants_who_are_male)),
    cohort_txt    = if ("cohort_s" %in% names(.)) as.character(cohort_s) else NA_character_
  ) %>%
  distinct(sampleset_id, .keep_all = TRUE)

# Per-PPM metadata (covariates, OR, HR) for Sections 10, 12, 13
ppm_meta <- read_csv(fp_ppm, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
  janitor::clean_names() %>%
  transmute(
    performance_id = as.character(pgs_performance_metric_ppm_id),
    covariates_txt = covariates_included_in_the_model,
    other_metric_txt = other_metric_s,
    other_info_txt = pgs_performance_other_relevant_information,
    or_txt = odds_ratio_or,
    hr_txt = hazard_ratio_hr
  ) %>%
  distinct(performance_id, .keep_all = TRUE)

# Per-PGS development metadata for Section 11
scores_meta <- read_csv(fp_scores, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
  janitor::clean_names() %>%
  transmute(
    pgs_id = as.character(polygenic_score_pgs_id),
    pgs_name = pgs_name,
    dev_method = pgs_development_method,
    n_variants = suppressWarnings(as.numeric(number_of_variants)),
    anc_gwas_txt = ancestry_distribution_percent_source_of_variant_associations_gwas,
    anc_dev_txt  = ancestry_distribution_percent_score_development_training,
    release_date = release_date
  ) %>%
  distinct(pgs_id, .keep_all = TRUE)

eval_df <- eval_raw %>%
  mutate(
    pgs_id = as.character(pgs_id),
    sampleset_id = as.character(sampleset_id),
    performance_id = as.character(performance_id),
    reported_trait = as.character(reported_trait),
    train_bucket = as.character(train_bucket),
    ci_source = as.character(ci_source),
    ancestry_display = as.character(ancestry_eval),
    raw_label = dplyr::coalesce(na_if(trimws(reported_trait), ""), NA_character_),
    trait_label = harmonize_trait_label(NA_character_, raw_label),
    has_auc = is.finite(auc),
    has_ci  = is.finite(estimate_ci_lower) & is.finite(estimate_ci_upper),
    ci_class = case_when(
      ci_source == "reported_native_or_parsed" ~ "Reported in Catalog",
      ci_source == "inferred_hm" ~ "Hanley-McNeil recovery",
      TRUE ~ "No usable CI"
    )
  ) %>%
  left_join(pss, by = "sampleset_id") %>%
  left_join(ppm_meta, by = "performance_id") %>%
  mutate(
    catalog_all_of_us = grepl("all of us|allofus|\\baou\\b", coalesce(cohort_txt, ""), ignore.case = TRUE)
  )

write_tbl(
  eval_df %>%
    filter(catalog_all_of_us) %>%
    transmute(
      performance_id, pgs_id, sampleset_id, reported_trait, ancestry_display,
      train_bucket, auc, ci_source, cohort_txt,
      note = "Catalog evaluation whose cohort is All of Us; not the manuscript All of Us analysis"
    ),
  "00_catalog_all_of_us_evaluations"
)
write_tbl(
  tibble(
    n_catalog_all_of_us_evaluations = sum(eval_df$catalog_all_of_us, na.rm = TRUE),
    n_catalog_all_of_us_with_usable_ci = sum(eval_df$catalog_all_of_us & eval_df$has_auc & eval_df$has_ci, na.rm = TRUE),
    n_catalog_all_of_us_sample_sets = n_distinct(eval_df$sampleset_id[eval_df$catalog_all_of_us]),
    note = "These published Catalog evaluations are not the manuscript's own All of Us analysis."
  ),
  "00_catalog_all_of_us_summary"
)

if (!"covariate_class" %in% names(stage1) || !"covariate_class" %in% names(stage1_f)) {
  stop(
    "Stage-1 tables have no covariate_class. Re-run 06_meta_ivw.R and 07_i_square_filter.R.",
    call. = FALSE
  )
}
stage1   <- attach_train_bucket(stage1, eval_df)
stage1_f <- attach_train_bucket(stage1_f, eval_df)

forest_cells_pub <- forest_raw %>%
  distinct(trait_label, trained_bucket, target_ancestry, n_pairs, delta_hat, se_hat, lo, hi, all_not_pooled)

keep_traits <- sort(unique(forest_cells_pub$trait_label))
message(glue(">> Figure 1C traits: {length(keep_traits)}; unique cells in published file: {nrow(forest_cells_pub)}"))

# Rebuild paired ΔAUC from the I²-filtered Stage-1 table with the submitted
# logit-scale SEs and with the published delta-method SEs.
rebuilt     <- rebuild_paired_delta(stage1_f, keep_traits, se_scale = "logit_as_implemented")
paired_main <- rebuilt$paired
delta_main  <- rebuilt$cells

rebuilt_dm  <- rebuild_paired_delta(stage1_f, keep_traits, se_scale = "delta_method")
paired_dm   <- rebuilt_dm$paired
delta_dm    <- rebuilt_dm$cells

rebuilt_ci  <- rebuild_paired_delta(stage1_f, keep_traits, se_scale = "from_ci")
delta_ci    <- rebuilt_ci$cells

fp_kept <- file.path(catalog_results_dir(stamp), "stage1", paste0("stage1_kept_evaluations_", stamp, ".csv"))
if (!file.exists(fp_kept)) {
  stop("Kept-evaluation table not found: ", fp_kept, ". Re-run 06_meta_ivw.R.", call. = FALSE)
}
kept_eval <- read_csv(fp_kept, show_col_types = FALSE)
usable_keys <- kept_eval %>%
  transmute(
    performance_id = as.character(performance_id),
    pgs_id = as.character(pgs_id),
    trait_label = as.character(trait_label),
    ancestry_display = as.character(ancestry_display),
    sampleset_id = as.character(sampleset_id),
    covariate_class = as.character(covariate_class),
    covariates = as.character(covariates)
  ) %>%
  left_join(eval_df %>% distinct(performance_id, ci_source), by = "performance_id")
if (anyNA(usable_keys$ci_source)) {
  stop("A kept evaluation has no CI source in the audit table.", call. = FALSE)
}
kept_n <- usable_keys %>%
  count(pgs_id, trait_label, ancestry_display, covariate_class, name = "n_kept")
k_check <- stage1 %>%
  mutate(
    pgs_id = as.character(pgs_id),
    trait_label = as.character(trait_label),
    ancestry_display = as.character(ancestry_display),
    covariate_class = as.character(covariate_class)
  ) %>%
  left_join(kept_n, by = c("pgs_id", "trait_label", "ancestry_display", "covariate_class")) %>%
  mutate(n_kept = coalesce(n_kept, 0L))
if (any(k_check$k_eval != k_check$n_kept)) {
  bad <- k_check %>%
    filter(k_eval != n_kept) %>%
    select(pgs_id, trait_label, ancestry_display, covariate_class, k_eval, n_kept)
  stop(
    "Kept evaluations do not match the Stage-1 average.\n",
    paste(capture.output(print(head(bad, 20))), collapse = "\n"),
    call. = FALSE
  )
}

# ===========================================================================
# 1) Pipeline funnel
# ===========================================================================

message(">> [1] Pipeline funnel")

n_ppm       <- nrow(eval_df)
n_auc       <- sum(eval_df$has_auc)
n_reported  <- sum(eval_df$ci_source == "reported_native_or_parsed", na.rm = TRUE)
n_inferred  <- sum(eval_df$ci_source == "inferred_hm", na.rm = TRUE)
n_usable    <- sum(eval_df$has_auc & eval_df$has_ci)
n_auc_no_ci <- n_auc - n_usable
n_s1        <- nrow(stage1)
n_s1f       <- nrow(stage1_f)
n_i2_drop   <- n_s1 - n_s1f
n_forest_traits <- n_distinct(forest_cells_pub$trait_label)
n_forest_cells  <- nrow(forest_cells_pub)

funnel <- tibble(
  step = c(
    "1. All PGS Catalog performance records (PPM)",
    "2. Records with a finite AUC",
    "3a. AUC with a published/parsed 95% CI",
    "3b. AUC with Hanley-McNeil recovered CI",
    "3c. AUC dropped (no published CI and not recoverable)",
    "4. Evaluations entering IVW (usable 95% CI)",
    "5. Stage-1 averages (trait x score x ancestry x covariate class)",
    "6. Stage-1 cells after I2 > 80% filter",
    "7. Figure 1C traits / trait-bucket-ancestry cells"
  ),
  unit = c(
    "evaluations", "evaluations", "evaluations", "evaluations", "evaluations",
    "evaluations", "cells", "cells", "traits / cells"
  ),
  n = c(n_ppm, n_auc, n_reported, n_inferred, n_auc_no_ci, n_usable, n_s1, n_s1f, NA_real_),
  n_label = c(
    as.character(n_ppm), as.character(n_auc), as.character(n_reported),
    as.character(n_inferred), as.character(n_auc_no_ci), as.character(n_usable),
    as.character(n_s1), as.character(n_s1f),
    paste0(n_forest_traits, " / ", n_forest_cells)
  ),
  catalog_snapshot = stamp
)

write_tbl(funnel, "01_pipeline_funnel")

funnel_plot_df <- tibble(
  step = factor(
    c("All PPM", "Finite AUC", "Reported CI", "Hanley-McNeil", "Usable CI",
      "Stage-1 cells", "After I2 filter", "Figure 1C cells"),
    levels = c("All PPM", "Finite AUC", "Reported CI", "Hanley-McNeil", "Usable CI",
               "Stage-1 cells", "After I2 filter", "Figure 1C cells")
  ),
  n = c(n_ppm, n_auc, n_reported, n_inferred, n_usable, n_s1, n_s1f, n_forest_cells),
  group = c("Evaluations", "Evaluations", "Evaluations", "Evaluations", "Evaluations",
            "Pooled cells", "Pooled cells", "Pooled cells")
)

p_funnel <- ggplot(funnel_plot_df, aes(x = n, y = fct_rev(step), fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(n)), hjust = -0.1, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)), labels = comma) +
  scale_fill_manual(values = c("Evaluations" = "#4E79A7", "Pooled cells" = "#F28E2B"), name = NULL) +
  labs(
    title = glue("PGS Catalog analysis funnel (stamp {stamp})"),
    subtitle = "What entered the IVW meta-analysis and Figure 1C",
    x = "Count", y = NULL
  ) +
  theme_reviewer()
save_plot("01_pipeline_funnel", p_funnel, width = 9, height = 5)

# ===========================================================================
# 2) CI source by ancestry
# ===========================================================================

message(">> [2] CI source by ancestry")

ci_by_anc <- eval_df %>%
  filter(has_auc) %>%
  mutate(
    ancestry_display = factor(
      ancestry_display,
      levels = unique(c(display_levels, sort(unique(as.character(ancestry_display)))))
    ),
    ci_class = factor(ci_class, levels = c("Reported in Catalog", "Hanley-McNeil recovery", "No usable CI"))
  ) %>%
  count(ancestry_display, ci_class, name = "n") %>%
  group_by(ancestry_display) %>%
  mutate(n_auc = sum(n), share = n / n_auc) %>%
  ungroup()

write_tbl(ci_by_anc, "02_ci_source_by_ancestry")

ci_plot_df <- ci_by_anc %>%
  mutate(anc_plot = ifelse(n_auc < 10, "Other ancestry labels", as.character(ancestry_display))) %>%
  group_by(anc_plot, ci_class) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(anc_plot) %>%
  mutate(n_auc = sum(n), share = n / n_auc) %>%
  ungroup() %>%
  mutate(
    anc_plot = factor(
      anc_plot,
      levels = c(intersect(display_levels, unique(anc_plot)),
                 setdiff(sort(unique(anc_plot)), c(display_levels, "Other ancestry labels")),
                 "Other ancestry labels")
    ),
    anc_label = glue("{anc_plot}  (n = {comma(n_auc)})"),
    anc_label = fct_reorder(anc_label, as.integer(anc_plot))
  )

p_ci <- ggplot(ci_plot_df, aes(x = share, y = fct_rev(anc_label), fill = ci_class)) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = ifelse(share >= 0.08, percent(share, accuracy = 1), "")),
    position = position_stack(vjust = 0.5), size = 3, colour = "white"
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(
    values = c("Reported in Catalog" = "#4E79A7", "Hanley-McNeil recovery" = "#F28E2B", "No usable CI" = "#BAB0AC"),
    name = "CI source"
  ) +
  labs(
    title = "Where Catalog AUC variances came from",
    subtitle = glue(
      "Among {comma(n_auc)} finite AUCs: {comma(n_reported)} reported CIs, ",
      "{comma(n_inferred)} Hanley-McNeil recoveries, {comma(n_auc_no_ci)} dropped"
    ),
    x = "Share of evaluations with a finite AUC", y = NULL,
    caption = "n = evaluations with a finite AUC per evaluation-ancestry label. 'Other ancestry labels' pools labels with fewer than 10 evaluations."
  ) +
  theme_reviewer()
save_plot("02_ci_source_by_ancestry", p_ci, width = 9, height = 5.2)

# ===========================================================================
# 3) Sample sizes
# ===========================================================================

message(">> [3] Sample sizes")

usable <- eval_df %>% filter(has_auc, has_ci)

sample_overall <- usable %>%
  mutate(ancestry_display = ifelse(is.na(ancestry_display) | ancestry_display == "", "Not reported", ancestry_display)) %>%
  group_by(ancestry_display) %>%
  summarise(
    n_evaluations = n(),
    n_distinct_pgs = n_distinct(pgs_id),
    n_distinct_pss = n_distinct(sampleset_id),
    n_with_case_control = sum(is.finite(n_cases) & is.finite(n_controls) & n_cases > 0 & n_controls > 0),
    median_n_cases = median(n_cases, na.rm = TRUE),
    iqr_n_cases_lo = quantile(n_cases, 0.25, na.rm = TRUE),
    iqr_n_cases_hi = quantile(n_cases, 0.75, na.rm = TRUE),
    min_n_cases = min(n_cases, na.rm = TRUE),
    median_n_controls = median(n_controls, na.rm = TRUE),
    min_n_controls = min(n_controls, na.rm = TRUE),
    median_n_individuals = median(n_individuals, na.rm = TRUE),
    median_case_fraction = median(n_cases / (n_cases + n_controls), na.rm = TRUE),
    .groups = "drop"
  )

write_tbl(sample_overall, "03_sample_sizes_overall_by_ancestry")

sample_forest <- usable %>%
  filter(trait_label %in% keep_traits) %>%
  group_by(trait_label, ancestry_display) %>%
  summarise(
    n_evaluations = n(),
    n_distinct_pgs = n_distinct(pgs_id),
    n_reported_ci = sum(ci_source == "reported_native_or_parsed"),
    n_inferred_ci = sum(ci_source == "inferred_hm"),
    n_with_case_control = sum(is.finite(n_cases) & is.finite(n_controls) & n_cases > 0 & n_controls > 0),
    median_n_cases = median(n_cases, na.rm = TRUE),
    iqr_n_cases_lo = as.numeric(quantile(n_cases, 0.25, na.rm = TRUE)),
    iqr_n_cases_hi = as.numeric(quantile(n_cases, 0.75, na.rm = TRUE)),
    min_n_cases = suppressWarnings(min(n_cases, na.rm = TRUE)),
    median_n_controls = median(n_controls, na.rm = TRUE),
    min_n_controls = suppressWarnings(min(n_controls, na.rm = TRUE)),
    median_n_individuals = median(n_individuals, na.rm = TRUE),
    median_case_fraction = median(n_cases / (n_cases + n_controls), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    min_n_cases = ifelse(is.infinite(min_n_cases), NA_real_, min_n_cases),
    min_n_controls = ifelse(is.infinite(min_n_controls), NA_real_, min_n_controls)
  ) %>%
  arrange(trait_label, ancestry_display)

write_tbl(sample_forest, "03_sample_sizes_forest_traits")

trait_order <- sample_forest %>%
  filter(ancestry_display %in% anc_keep) %>%
  group_by(trait_label) %>%
  summarise(n_total = sum(n_evaluations), .groups = "drop") %>%
  arrange(desc(n_total)) %>%
  pull(trait_label)
trait_levels_wrapped <- wrap_trait(trait_order)

heat_df <- sample_forest %>%
  filter(ancestry_display %in% anc_keep) %>%
  mutate(
    ancestry_display = factor(ancestry_display, levels = anc_keep),
    trait_wrapped = factor(wrap_trait(trait_label), levels = rev(trait_levels_wrapped)),
    cases_lab = coalesce(fmt_cases(median_n_cases), "n/a"),
    n_dark = n_evaluations > 60,
    cases_dark = is.finite(median_n_cases) & median_n_cases > 3000
  )

text_pal <- c(`TRUE` = "white", `FALSE` = "black")

p_n_a <- ggplot(heat_df, aes(x = ancestry_display, y = trait_wrapped, fill = n_evaluations)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = n_evaluations, colour = n_dark), size = 3.1) +
  scale_fill_gradient(low = "#EFF3FF", high = "#08519C", name = "Evaluations", breaks = c(1, 50, 100, 150)) +
  scale_colour_manual(values = text_pal, guide = "none") +
  scale_x_discrete(labels = function(x) str_wrap(x, 12), position = "top") +
  labs(title = "A  Evaluations with a usable 95% CI", x = NULL, y = NULL) +
  theme_reviewer(base_size = 9) +
  theme(legend.position = "right", panel.grid = element_blank(),
        axis.text.x.top = element_text(size = 8), axis.text.y = element_text(size = 8, lineheight = 0.85))

p_n_b <- ggplot(heat_df, aes(x = ancestry_display, y = trait_wrapped, fill = median_n_cases)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = cases_lab, colour = cases_dark), size = 3.1) +
  scale_fill_gradient(
    low = "#FEE6CE", high = "#A63603", trans = "log10",
    breaks = c(100, 300, 1000, 3000), labels = c("100", "300", "1k", "3k"),
    na.value = "grey90", name = "Median cases"
  ) +
  scale_colour_manual(values = text_pal, guide = "none") +
  scale_x_discrete(labels = function(x) str_wrap(x, 12), position = "top") +
  labs(title = "B  Median number of cases per evaluation", x = NULL, y = NULL) +
  theme_reviewer(base_size = 9) +
  theme(legend.position = "right", panel.grid = element_blank(),
        axis.text.x.top = element_text(size = 8), axis.text.y = element_blank())

p_n <- (p_n_a | p_n_b) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Figure 1C traits: how much Catalog evidence sits behind each trait x ancestry",
    caption = paste0(
      "Blank tile = no evaluation with a usable CI. ",
      "n/a = evaluations exist but the Catalog sample-set table has no case/control counts. ",
      "Traits ordered by total evaluations."
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 11), plot.caption = element_text(size = 8, hjust = 0))
  )
save_plot("03_sample_sizes_forest_traits", p_n, width = 12, height = 6.2)

# ===========================================================================
# 4) Independence between the two sides of a pair
# ===========================================================================

message(">> [4] Pair independence")

eval_keys <- usable_keys

independence_rows <- paired_main %>%
  mutate(pair_id = paste(trait_label, trained_bucket, target_ancestry, pgs_id, sep = "||")) %>%
  rowwise() %>%
  mutate(
    eur_perf = list(eval_keys$performance_id[eval_keys$pgs_id == pgs_id & eval_keys$trait_label == trait_label & eval_keys$ancestry_display == "European" & eval_keys$covariate_class == covariate_class]),
    tgt_perf = list(eval_keys$performance_id[eval_keys$pgs_id == pgs_id & eval_keys$trait_label == trait_label & eval_keys$ancestry_display == target_ancestry & eval_keys$covariate_class == covariate_class]),
    eur_pss  = list(eval_keys$sampleset_id[eval_keys$pgs_id == pgs_id & eval_keys$trait_label == trait_label & eval_keys$ancestry_display == "European" & eval_keys$covariate_class == covariate_class]),
    tgt_pss  = list(eval_keys$sampleset_id[eval_keys$pgs_id == pgs_id & eval_keys$trait_label == trait_label & eval_keys$ancestry_display == target_ancestry & eval_keys$covariate_class == covariate_class])
  ) %>%
  ungroup() %>%
  mutate(
    n_eur_eval = lengths(eur_perf),
    n_tgt_eval = lengths(tgt_perf),
    n_overlap_performance_id = map2_int(eur_perf, tgt_perf, ~ length(intersect(.x, .y))),
    n_overlap_sampleset_id   = map2_int(eur_pss, tgt_pss, ~ length(intersect(.x, .y))),
    target_is_multi_incl_eur = grepl("Multi-ancestry including European", target_ancestry, ignore.case = TRUE)
  ) %>%
  select(
    trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class,
    eur_k, tgt_k, n_eur_eval, n_tgt_eval,
    n_overlap_performance_id, n_overlap_sampleset_id,
    target_is_multi_incl_eur, delta, delta_se
  )

write_tbl(independence_rows, "04_independence_pairs")

independence_summary <- tibble(
  n_pgs_pairs = nrow(independence_rows),
  n_pairs_with_overlapping_performance_id = sum(independence_rows$n_overlap_performance_id > 0),
  n_pairs_with_overlapping_sampleset_id = sum(independence_rows$n_overlap_sampleset_id > 0),
  n_pairs_target_multi_incl_eur = sum(independence_rows$target_is_multi_incl_eur),
  n_figure1c_cells = n_forest_cells,
  n_figure1c_cells_multi_incl_eur = sum(forest_cells_pub$target_ancestry == "Multi-ancestry including European")
)

write_tbl(independence_summary, "04_independence_summary")

# ===========================================================================
# 5) What Figure 1C actually pools (published file)
#    Published intervals use the delta-method SE; Section 8 compares
#    the submitted logit-scale recalculation.
# ===========================================================================

message(">> [5] Figure 1C cells, recovery-pack layout")

fig1c <- forest_cells_pub %>%
  mutate(
    ci_excludes_zero = (lo > 0) | (hi < 0),
    pooling_status = ifelse(all_not_pooled, "All single-eval cells", "Includes pooled Stage-1 cells"),
    what_is_pooled = paste0(n_pairs, " PGS-level ΔAUC values (same PGS in EUR vs ", target_ancestry, ")"),
    se_scale_note = "Published Figure 1C uses the delta-method SE; see 08_scale_fix_cells for the submitted logit-scale recalculation"
  ) %>%
  arrange(trait_label, trained_bucket, target_ancestry)

cmp_rebuild <- fig1c %>%
  inner_join(
    delta_dm %>% select(trait_label, trained_bucket, target_ancestry, delta_hat_rebuild = delta_hat, se_hat_rebuild = se_hat),
    by = c("trait_label", "trained_bucket", "target_ancestry")
  ) %>%
  mutate(abs_diff = abs(delta_hat - delta_hat_rebuild), abs_diff_se = abs(se_hat - se_hat_rebuild))

write_tbl(fig1c, "05_figure1c_cells")
write_tbl(cmp_rebuild, "05_rebuild_vs_published_delta")

prep_forest <- function(df) {
  df %>%
    mutate(
      trait_wrapped = factor(wrap_trait(trait_label), levels = trait_levels_wrapped),
      target_ancestry = factor(as.character(target_ancestry), levels = rev(anc_keep)),
      trained_bucket = factor(trained_bucket, levels = names(pal_bucket)),
      pooled = factor(ifelse(all_not_pooled, pooled_levels[2], pooled_levels[1]), levels = pooled_levels),
      sig = (lo > 0) | (hi < 0),
      k_lab = paste0("k = ", n_pairs, ifelse(sig, " *", ""))
    )
}

fig1c_plot_df <- prep_forest(fig1c)

p_fig1c <- ggplot(fig1c_plot_df, aes(x = delta_hat, y = target_ancestry, colour = trained_bucket, group = trained_bucket)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.55, position = position_dodge(width = dodge_w), orientation = "y") +
  geom_point(aes(shape = pooled), size = 2.6, fill = "white", stroke = 0.9, position = position_dodge(width = dodge_w)) +
  geom_text(aes(x = x_lab_pos, label = k_lab, fontface = ifelse(sig, "bold", "plain")),
            size = 2.7, hjust = 0, colour = "grey20", position = position_dodge(width = dodge_w), show.legend = FALSE) +
  scale_colour_manual(values = pal_bucket, name = "PGS training weights") +
  scale_shape_manual(values = setNames(c(16, 21), pooled_levels), name = "Point fill") +
  scale_x_continuous(breaks = x_breaks_forest, labels = x_labels_forest) +
  coord_cartesian(xlim = x_lim_forest, clip = "off") +
  facet_grid(trait_wrapped ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title = "Figure 1C cells in the recovery-pack layout: within-PGS ΔAUC (target − European)",
    subtitle = "Intervals are the published delta-method standard errors. The script 08 forest remains the publication figure.",
    x = "ΔAUC (target − European) with 95% CI", y = NULL,
    caption = paste0(
      "k = number of PGS with a paired European/target estimate;  * = 95% CI excludes zero.\n",
      "Hollow points: every PGS in the cell has a single evaluation on each side (no Stage-1 pooling)."
    )
  ) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 15, size = 3.5)), shape = guide_legend(order = 2)) +
  theme_reviewer(base_size = 9) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8, face = "bold", lineheight = 0.85),
    strip.background = element_rect(fill = "grey95", colour = NA),
    panel.spacing.y = unit(0.35, "lines"),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(5.5, 50, 5.5, 5.5),
    plot.caption = element_text(hjust = 0, size = 7.5),
    plot.title = element_text(size = 11),
    legend.box = "horizontal",
    legend.spacing.x = unit(1.2, "lines")
  )
n_row_plot <- n_distinct(paste(fig1c$trait_label, fig1c$target_ancestry))
save_plot("05_figure1c_cells", p_fig1c, width = 10, height = max(6.5, 0.30 * n_row_plot + 2.8), limitsize = FALSE)

# ===========================================================================
# 6) Sensitivity: reported-CI-only (drop Hanley-McNeil recoveries)
#    Run twice: logit-scale SE as submitted (06_) and delta-method SE (06b_).
# ===========================================================================

message(">> [6] Hanley-McNeil sensitivity on both standard-error scales")

eval_reported <- eval_df %>%
  filter(has_auc, has_ci, ci_source == "reported_native_or_parsed") %>%
  filter(is.finite(auc), is.finite(estimate_ci_lower), is.finite(estimate_ci_upper)) %>%
  mutate(covariate_class = classify_covariates(covariates_txt)) %>%
  keep_one_evaluation(
    n_col = "n_individuals",
    cases_col = "n_cases",
    controls_col = "n_controls",
    covariates_col = "covariates_txt",
    ancestry_col = "ancestry_display",
    class_col = "covariate_class"
  )

stage1_rep   <- pool_stage1(eval_reported)
stage1_rep_f <- i2_filter(stage1_rep) %>% attach_train_bucket(eval_df)

run_hm_sensitivity <- function(baseline_cells, se_scale, prefix, scale_label) {
  rebuilt_rep <- rebuild_paired_delta(stage1_rep_f, keep_traits, se_scale = se_scale)

  sens <- baseline_cells %>%
    select(trait_label, trained_bucket, target_ancestry, n_pairs, delta_hat, se_hat, lo, hi, ci_excludes_zero) %>%
    rename(n_pairs_all = n_pairs, delta_all = delta_hat, se_all = se_hat, lo_all = lo, hi_all = hi, sig_all = ci_excludes_zero) %>%
    left_join(
      rebuilt_rep$cells %>%
        select(trait_label, trained_bucket, target_ancestry,
               n_pairs_reported = n_pairs, delta_reported = delta_hat, se_reported = se_hat,
               lo_reported = lo, hi_reported = hi, sig_reported = ci_excludes_zero),
      by = c("trait_label", "trained_bucket", "target_ancestry")
    ) %>%
    mutate(
      se_scale = se_scale,
      delta_diff = delta_reported - delta_all,
      sign_flip = !is.na(delta_reported) & (sign(delta_reported) != sign(delta_all)) & abs(delta_all) > 1e-8 & abs(delta_reported) > 1e-8,
      inference_flip = !is.na(sig_reported) & (sig_all != sig_reported)
    )

  write_tbl(rebuilt_rep$cells %>% mutate(se_scale = se_scale), paste0(prefix, "_reported_ci_only_delta"))
  write_tbl(sens, paste0(prefix, "_sensitivity_comparison"))

  sens_summary <- tibble(
    se_scale = se_scale,
    n_figure1c_cells = nrow(sens),
    n_cells_recovered_without_hanley_mcneil = sum(!is.na(sens$delta_reported)),
    n_cells_lost_without_hanley_mcneil = sum(is.na(sens$delta_reported)),
    n_sign_flips = sum(sens$sign_flip, na.rm = TRUE),
    n_inference_flips = sum(sens$inference_flip, na.rm = TRUE),
    n_sig_all = sum(sens$sig_all, na.rm = TRUE),
    n_sig_reported_only = sum(sens$sig_reported, na.rm = TRUE),
    max_abs_delta_diff = max(abs(sens$delta_diff), na.rm = TRUE)
  )
  write_tbl(sens_summary, paste0(prefix, "_sensitivity_summary"))

  sens_long <- bind_rows(
    baseline_cells %>% transmute(trait_label, trained_bucket, target_ancestry, n_pairs, delta_hat, lo, hi,
                                 source = "All usable CIs (reported + Hanley-McNeil)"),
    rebuilt_rep$cells %>% transmute(trait_label, trained_bucket, target_ancestry, n_pairs, delta_hat, lo, hi,
                                    source = "Reported CIs only")
  ) %>%
    mutate(source = factor(source, levels = c("All usable CIs (reported + Hanley-McNeil)", "Reported CIs only")))

  lost_cells <- sens %>% filter(is.na(delta_reported)) %>% select(trait_label, trained_bucket, target_ancestry)

  p_sens <- forest_compare_plot(
    sens_long,
    pal = c("#4E79A7", "#E15759"), shapes = c(16, 17),
    title = glue("Sensitivity: ΔAUC with vs without Hanley-McNeil recovered CIs ({scale_label})"),
    subtitle = glue(
      "{sens_summary$n_cells_recovered_without_hanley_mcneil} of {sens_summary$n_figure1c_cells} cells remain estimable with reported CIs only; ",
      "{sens_summary$n_sign_flips} sign flips; {sens_summary$n_inference_flips} changes in whether the 95% CI excludes zero."
    ),
    caption = paste0(
      "k = number of PGS with a paired European/target estimate;  * = 95% CI excludes zero.\n",
      "'not estimable' = every PGS in the cell depended on a Hanley-McNeil recovered CI."
    ),
    lost_df = lost_cells, trait_levels = trait_levels_wrapped
  )
  save_plot(paste0(prefix, "_sensitivity_forest"), p_sens, width = 12, height = max(6.5, 0.30 * n_row_plot + 2.6), limitsize = FALSE)

  p_sc <- ggplot(sens %>% filter(!is.na(delta_reported)), aes(x = delta_all, y = delta_reported, color = target_ancestry)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_point(size = 2) +
    scale_color_manual(values = pal_ancestry, name = "Target ancestry") +
    labs(
      title = glue("Reported-CI-only vs all-usable-CI ΔAUC ({scale_label})"),
      subtitle = glue("{sens_summary$n_inference_flips} / {sens_summary$n_cells_recovered_without_hanley_mcneil} overlapping cells change whether the 95% CI excludes zero"),
      x = "ΔAUC using all usable CIs", y = "ΔAUC using reported CIs only"
    ) +
    guides(color = guide_legend(nrow = 2)) +
    theme_reviewer() +
    coord_equal()
  save_plot(paste0(prefix, "_sensitivity_scatter"), p_sc, width = 7.5, height = 6.6)

  list(sens = sens, summary = sens_summary)
}

hm_impl <- run_hm_sensitivity(delta_main, "logit_as_implemented", "06",  "logit-scale SE as submitted")
hm_dm   <- run_hm_sensitivity(delta_dm,   "delta_method",         "06b", "delta-method SE")

# ===========================================================================
# 7) Within-cell dependence: do the k PGS pooled in a cell share sample sets?
# ===========================================================================

message(">> [7] Within-cell sample-set sharing")

pair_sides <- paired_main %>%
  select(trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class)

sample_set_counts <- function(keys_side) {
  keys_side %>%
    group_by(trait_label, trained_bucket, target_ancestry, sampleset_id) %>%
    summarise(n_pgs = n_distinct(pgs_id), .groups = "drop") %>%
    group_by(trait_label, trained_bucket, target_ancestry) %>%
    summarise(
      n_sample_sets = n(),
      max_pgs_sharing_one_sample_set = max(n_pgs),
      .groups = "drop"
    )
}

eur_sets <- pair_sides %>%
  inner_join(
    usable_keys %>%
      filter(ancestry_display == "European") %>%
      select(pgs_id, trait_label, covariate_class, sampleset_id),
    by = c("pgs_id", "trait_label", "covariate_class"),
    relationship = "many-to-many"
  ) %>%
  sample_set_counts() %>%
  rename(
    eur_n_sample_sets = n_sample_sets,
    eur_max_pgs_sharing_one_sample_set = max_pgs_sharing_one_sample_set
  )

tgt_sets <- pair_sides %>%
  inner_join(
    usable_keys %>%
      select(pgs_id, trait_label, ancestry_display, covariate_class, sampleset_id),
    by = c("pgs_id", "trait_label", "covariate_class", "target_ancestry" = "ancestry_display")
  ) %>%
  sample_set_counts() %>%
  rename(
    tgt_n_sample_sets = n_sample_sets,
    tgt_max_pgs_sharing_one_sample_set = max_pgs_sharing_one_sample_set
  )

cell_pss_sharing <- paired_main %>%
  group_by(trait_label, trained_bucket, target_ancestry) %>%
  summarise(n_pairs = n_distinct(pgs_id), .groups = "drop") %>%
  left_join(eur_sets, by = c("trait_label", "trained_bucket", "target_ancestry")) %>%
  left_join(tgt_sets, by = c("trait_label", "trained_bucket", "target_ancestry")) %>%
  mutate(
    eur_n_sample_sets = coalesce(eur_n_sample_sets, 0L),
    eur_max_pgs_sharing_one_sample_set = coalesce(eur_max_pgs_sharing_one_sample_set, 0L),
    tgt_n_sample_sets = coalesce(tgt_n_sample_sets, 0L),
    tgt_max_pgs_sharing_one_sample_set = coalesce(tgt_max_pgs_sharing_one_sample_set, 0L),
    all_pgs_share_one_target_sample_set = n_pairs >= 2 & tgt_max_pgs_sharing_one_sample_set == n_pairs,
    all_pgs_share_one_eur_sample_set    = n_pairs >= 2 & eur_max_pgs_sharing_one_sample_set == n_pairs
  ) %>%
  arrange(desc(n_pairs), trait_label, trained_bucket, target_ancestry)

write_tbl(cell_pss_sharing, "07_within_cell_sampleset_sharing")

sharing_summary <- tibble(
  n_cells_total = nrow(cell_pss_sharing),
  n_cells_k_ge_2 = sum(cell_pss_sharing$n_pairs >= 2),
  n_cells_k_ge_2_all_pgs_same_target_pss = sum(cell_pss_sharing$all_pgs_share_one_target_sample_set),
  n_cells_k_ge_2_all_pgs_same_eur_pss = sum(cell_pss_sharing$all_pgs_share_one_eur_sample_set),
  n_cells_k_ge_2_all_pgs_same_pss_both_sides = sum(cell_pss_sharing$all_pgs_share_one_target_sample_set & cell_pss_sharing$all_pgs_share_one_eur_sample_set),
  n_cells_k_ge_2_target_from_single_pss = sum(cell_pss_sharing$n_pairs >= 2 & cell_pss_sharing$tgt_n_sample_sets == 1)
)
write_tbl(sharing_summary, "07_within_cell_sampleset_sharing_summary")

# ===========================================================================
# 8) Stage-2 scale comparison
#    (a) PGS000004 African-vs-European interval check
#    (b) rebuild all cells with three SE conventions
#    (c) mixed-scale per-PGS bounds used on the submitted scatter
# ===========================================================================

message(">> [8] Stage-2 scale comparison")

# (a) Interval-formula check on the published PGS000004 inputs, then this
#     download's Stage-1 rows (the European Breast cancer cell is I2-dropped).
ex_afr_eta <- 0.291508
ex_eur_eta <- 0.405465
ex_afr_se  <- 0.01630624
ex_eur_se  <- 0.02126195
ex_afr_auc <- plogis(ex_afr_eta)
ex_eur_auc <- plogis(ex_eur_eta)
ex_d       <- ex_afr_auc - ex_eur_auc
ex_se_impl <- sqrt(ex_afr_se^2 + ex_eur_se^2)
ex_se_dm   <- sqrt((ex_afr_se * ex_afr_auc * (1 - ex_afr_auc))^2 +
                     (ex_eur_se * ex_eur_auc * (1 - ex_eur_auc))^2)
example_check <- tibble(
  quantity = c(
    "African pooled logit(AUC)", "European pooled logit(AUC)",
    "African logit-scale SE", "European logit-scale SE",
    "African AUC", "European AUC", "delta AUC",
    "SE(delta) as implemented (logit SEs)", "95% CI lo as implemented", "95% CI hi as implemented",
    "SE(delta) delta method", "95% CI lo delta method", "95% CI hi delta method"
  ),
  recovered = c(
    ex_afr_eta, ex_eur_eta, ex_afr_se, ex_eur_se, ex_afr_auc, ex_eur_auc, ex_d,
    ex_se_impl, ex_d - 1.96 * ex_se_impl, ex_d + 1.96 * ex_se_impl,
    ex_se_dm, ex_d - 1.96 * ex_se_dm, ex_d + 1.96 * ex_se_dm
  ),
  expected = c(
    0.291508, 0.405465, 0.01630624, 0.02126195, NA, NA, -0.027635,
    NA, -0.080153, 0.024883,
    NA, -0.040332, -0.014937
  )
) %>%
  mutate(abs_diff = abs(recovered - expected), matches_expected = ifelse(is.na(expected), NA, abs_diff < 1e-5))
write_tbl(example_check, "08_scale_fix_PGS000004_example")
stopifnot(all(example_check$matches_expected, na.rm = TRUE))

ex_this_download <- bind_rows(
  stage1 %>%
    filter(pgs_id == "PGS000004", trait_label == "Breast cancer", ancestry_display %in% c("African", "European")) %>%
    mutate(stage1_table = "unfiltered"),
  stage1_f %>%
    filter(pgs_id == "PGS000004", trait_label == "Breast cancer", ancestry_display %in% c("African", "European")) %>%
    mutate(stage1_table = "I2_filtered")
) %>%
  transmute(
    catalog_snapshot = stamp, stage1_table, pgs_id, trait_label, ancestry_display,
    eta, se, I2, k_eval, auc_pooled, lo_pooled, hi_pooled,
    note = "On this download the European Breast cancer cell has I2 > 80 and is dropped before Figure 1C."
  )
write_tbl(ex_this_download, "08_scale_fix_PGS000004_this_download")

# (b) All cells under the three SE conventions
scale_cells <- delta_main %>%
  transmute(trait_label, trained_bucket, target_ancestry, n_pairs, all_not_pooled,
            delta_impl = delta_hat, se_impl = se_hat, lo_impl = lo, hi_impl = hi, sig_impl = ci_excludes_zero) %>%
  inner_join(
    delta_dm %>% transmute(trait_label, trained_bucket, target_ancestry,
                           delta_dm = delta_hat, se_dm = se_hat, lo_dm = lo, hi_dm = hi, sig_dm = ci_excludes_zero),
    by = c("trait_label", "trained_bucket", "target_ancestry")
  ) %>%
  inner_join(
    delta_ci %>% transmute(trait_label, trained_bucket, target_ancestry,
                           delta_ci = delta_hat, se_ci = se_hat, lo_ci = lo, hi_ci = hi, sig_ci = ci_excludes_zero),
    by = c("trait_label", "trained_bucket", "target_ancestry")
  ) %>%
  mutate(
    se_inflation_factor = se_impl / se_dm,
    abs_delta_shift = abs(delta_impl - delta_dm),
    inference_change = case_when(
      !sig_impl & sig_dm & delta_dm < 0 ~ "gains significance (deficit)",
      !sig_impl & sig_dm & delta_dm > 0 ~ "gains significance (advantage)",
      sig_impl & !sig_dm ~ "loses significance",
      TRUE ~ "unchanged"
    ),
    dm_vs_ci_agree = sig_dm == sig_ci
  ) %>%
  arrange(delta_dm)

write_tbl(scale_cells, "08_scale_fix_cells")

scale_summary <- tibble(
  n_cells = nrow(scale_cells),
  n_sig_as_implemented = sum(scale_cells$sig_impl),
  n_sig_delta_method = sum(scale_cells$sig_dm),
  n_sig_delta_method_negative = sum(scale_cells$sig_dm & scale_cells$delta_dm < 0),
  n_sig_delta_method_positive = sum(scale_cells$sig_dm & scale_cells$delta_dm > 0),
  n_sig_delta_method_k1 = sum(scale_cells$sig_dm & scale_cells$n_pairs == 1),
  n_sig_from_ci = sum(scale_cells$sig_ci),
  n_cells_dm_ci_disagree = sum(!scale_cells$dm_vs_ci_agree),
  se_inflation_min = min(scale_cells$se_inflation_factor),
  se_inflation_max = max(scale_cells$se_inflation_factor),
  max_abs_delta_shift = max(scale_cells$abs_delta_shift)
)
write_tbl(scale_summary, "08_scale_fix_summary")

scale_long <- bind_rows(
  delta_main %>% transmute(trait_label, trained_bucket, target_ancestry, n_pairs, delta_hat, lo, hi, source = "As submitted (logit-scale SE)"),
  delta_dm   %>% transmute(trait_label, trained_bucket, target_ancestry, n_pairs, delta_hat, lo, hi, source = "Corrected (delta-method SE)")
) %>%
  mutate(source = factor(source, levels = c("As submitted (logit-scale SE)", "Corrected (delta-method SE)")))

p_scale <- forest_compare_plot(
  scale_long,
  pal = c("#BAB0AC", "#0072B2"), shapes = c(1, 16),
  title = "Figure 1C before and after putting Stage-2 standard errors on the AUC scale",
  subtitle = glue(
    "Cells with 95% CI excluding zero: {scale_summary$n_sig_as_implemented} as submitted -> {scale_summary$n_sig_delta_method} corrected ",
    "({scale_summary$n_sig_delta_method_negative} deficits, {scale_summary$n_sig_delta_method_positive} advantages). ",
    "SE inflation as submitted: {round(scale_summary$se_inflation_min, 1)}-{round(scale_summary$se_inflation_max, 1)}x."
  ),
  caption = paste0(
    "Point estimates differ by at most ", signif(scale_summary$max_abs_delta_shift, 2), " (weights change with the SE). ",
    "k = number of PGS per cell; * = 95% CI excludes zero (bold label refers to the corrected version).\n",
    "Corrected SE_AUC = SE_logit x AUC(1-AUC); identical inference when SE is derived from the back-transformed CI (", scale_summary$n_cells_dm_ci_disagree, " disagreements)."
  ),
  trait_levels = trait_levels_wrapped
)
save_plot("08_scale_fix_forest", p_scale, width = 12, height = max(6.5, 0.30 * n_row_plot + 2.6), limitsize = FALSE)

# (c) Mixed-scale per-PGS bounds used for scatter error bars (lines 285/297)
mixed_bounds <- paired_main %>%
  transmute(
    trait_label, trained_bucket, target_ancestry, pgs_id,
    eur_auc, eur_lo_mixed = eur_auc - 1.96 * eur_se_logit, eur_hi_mixed = eur_auc + 1.96 * eur_se_logit, eur_lo_correct = eur_lo, eur_hi_correct = eur_hi,
    tgt_auc, tgt_lo_mixed = tgt_auc - 1.96 * tgt_se_logit, tgt_hi_mixed = tgt_auc + 1.96 * tgt_se_logit, tgt_lo_correct = tgt_lo, tgt_hi_correct = tgt_hi
  ) %>%
  mutate(
    any_mixed_bound_outside_0_1 = eur_lo_mixed < 0 | eur_hi_mixed > 1 | tgt_lo_mixed < 0 | tgt_hi_mixed > 1,
    max_abs_bound_error = pmax(abs(eur_lo_mixed - eur_lo_correct), abs(eur_hi_mixed - eur_hi_correct),
                               abs(tgt_lo_mixed - tgt_lo_correct), abs(tgt_hi_mixed - tgt_hi_correct))
  )
write_tbl(mixed_bounds, "08_scale_fix_mixed_scale_pgs_bounds")

# ===========================================================================
# 9) Dependence sensitivities on the delta-method SEs
#    Fixed effect of every PGS in the cell.
#    When k > 1 but every score shares one target sample set, the sensitivity
#    is that single retained PGS. DerSimonian-Laird is not defined for one ΔAUC.
#    When at least two target sample sets remain, random effects are fit to
#    those independent ΔAUC values (one PGS per target sample set, lowest PGS ID).
#    The fixed effect of that reduced set is drawn only when it drops a duplicate.
# ===========================================================================

message(">> [9] Dependence sensitivities")

pair_class <- paired_dm %>%
  select(trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class)

pss_key_for <- function(pairs, ancestry_match) {
  pairs %>%
    inner_join(
      ancestry_match,
      by = c("pgs_id", "trait_label", "covariate_class"),
      relationship = "many-to-many"
    ) %>%
    group_by(trait_label, trained_bucket, target_ancestry, pgs_id) %>%
    summarise(pss_key = paste(sort(unique(sampleset_id)), collapse = ";"), .groups = "drop")
}

eur_pss_key <- pss_key_for(
  pair_class,
  usable_keys %>%
    filter(ancestry_display == "European") %>%
    select(pgs_id, trait_label, covariate_class, sampleset_id)
) %>%
  rename(eur_pss_key = pss_key)

tgt_pss_key <- pair_class %>%
  inner_join(
    usable_keys %>% select(pgs_id, trait_label, ancestry_display, covariate_class, sampleset_id),
    by = c("pgs_id", "trait_label", "covariate_class", "target_ancestry" = "ancestry_display")
  ) %>%
  group_by(trait_label, trained_bucket, target_ancestry, pgs_id) %>%
  summarise(tgt_pss_key = paste(sort(unique(sampleset_id)), collapse = ";"), .groups = "drop")

pair_tgt_pss <- eur_pss_key %>%
  inner_join(tgt_pss_key, by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id"))

paired_one_per_pss <- paired_dm %>%
  inner_join(pair_tgt_pss, by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")) %>%
  group_by(trait_label, trained_bucket, target_ancestry, tgt_pss_key) %>%
  arrange(pgs_id, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

cells_one_per_pss <- pool_stage2(paired_one_per_pss, model = "fixed") %>%
  mutate(sensitivity = "one PGS per target sample set")
# Random effects only on the independent ΔAUCs, and only when at least two remain.
# A cell with one ΔAUC (k = 1, or several PGS on one target sample set) has no τ² to estimate.
cells_random <- pool_stage2(paired_one_per_pss, model = "random_dl") %>%
  filter(n_pairs >= 2) %>%
  mutate(sensitivity = "random effects (DL)")

src_fixed   <- "IVW fixed effect"
src_single  <- "Single PGS (shared target sample set)"
src_reduced <- "Fixed effect, one PGS per sample set"
src_random  <- "Random effects (DL)"

dependence_cells <- delta_dm %>%
  transmute(trait_label, trained_bucket, target_ancestry, n_pairs,
            delta_fixed = delta_hat, se_fixed = se_hat, lo_fixed = lo, hi_fixed = hi, sig_fixed = ci_excludes_zero) %>%
  left_join(
    cells_one_per_pss %>% transmute(trait_label, trained_bucket, target_ancestry, n_pairs_one_per_pss = n_pairs,
                                    delta_one_per_pss = delta_hat, se_one_per_pss = se_hat, lo_one_per_pss = lo, hi_one_per_pss = hi, sig_one_per_pss = ci_excludes_zero),
    by = c("trait_label", "trained_bucket", "target_ancestry")
  ) %>%
  left_join(
    cells_random %>% transmute(trait_label, trained_bucket, target_ancestry, tau2_random = tau2, I2_stage2 = I2_stage2,
                               delta_random = delta_hat, se_random = se_hat, lo_random = lo, hi_random = hi, sig_random = ci_excludes_zero),
    by = c("trait_label", "trained_bucket", "target_ancestry")
  ) %>%
  left_join(cell_pss_sharing %>% select(trait_label, trained_bucket, target_ancestry, eur_n_sample_sets, tgt_n_sample_sets,
                                        all_pgs_share_one_target_sample_set, all_pgs_share_one_eur_sample_set),
            by = c("trait_label", "trained_bucket", "target_ancestry")) %>%
  mutate(
    re_applicable = n_pairs_one_per_pss >= 2,
    n_delta_random = if_else(re_applicable, n_pairs_one_per_pss, NA_integer_),
    single_pgs_sensitivity = n_pairs > 1 & n_pairs_one_per_pss == 1,
    reduced_fixed_effect = n_pairs_one_per_pss >= 2 & n_pairs_one_per_pss < n_pairs,
    robust_fixed_and_random = re_applicable & sig_fixed & coalesce(sig_random, FALSE),
    robustness = case_when(
      n_pairs == 1 & sig_fixed ~ "significant; single ΔAUC (k = 1; random effects not applicable)",
      n_pairs == 1 & !sig_fixed ~ "not significant; single ΔAUC (k = 1; random effects not applicable)",
      single_pgs_sensitivity & sig_fixed & sig_one_per_pss ~ "significant in the fixed effect and in the single PGS for that target sample set",
      single_pgs_sensitivity & sig_fixed & !sig_one_per_pss ~ "significant in the fixed effect only; the single PGS for that target sample set is not",
      single_pgs_sensitivity & !sig_fixed ~ "not significant (corrected fixed effect); sensitivity is the single PGS, not random effects",
      re_applicable & n_pairs_one_per_pss == n_pairs & sig_fixed & sig_random ~ "significant in the fixed effect and in random effects",
      re_applicable & n_pairs_one_per_pss == n_pairs & sig_fixed & !sig_random ~ "significant in the fixed effect only; not in random effects",
      re_applicable & sig_fixed & sig_one_per_pss & sig_random ~ "significant in the fixed effect, the one-per-sample-set fixed effect, and random effects",
      re_applicable & sig_fixed & sig_one_per_pss & !sig_random ~ "significant in both fixed effects; not in random effects",
      re_applicable & sig_fixed & !sig_one_per_pss & sig_random ~ "significant in the full fixed effect and in random effects; not after one PGS per sample set",
      re_applicable & sig_fixed ~ "significant in the full fixed effect only; sensitive to dependence",
      TRUE ~ "not significant (corrected fixed effect)"
    )
  ) %>%
  arrange(delta_fixed)

write_tbl(dependence_cells, "09_dependence_cells")

dependence_summary <- tibble(
  n_cells = nrow(dependence_cells),
  n_sig_fixed_corrected = sum(dependence_cells$sig_fixed),
  n_cells_single_delta = sum(dependence_cells$n_pairs == 1),
  n_sig_single_delta = sum(dependence_cells$n_pairs == 1 & dependence_cells$sig_fixed),
  n_cells_single_pgs_sensitivity = sum(dependence_cells$single_pgs_sensitivity),
  n_sig_single_pgs_sensitivity = sum(dependence_cells$single_pgs_sensitivity & dependence_cells$sig_one_per_pss, na.rm = TRUE),
  n_cells_random_effects = sum(dependence_cells$re_applicable),
  n_sig_random_dl = sum(dependence_cells$sig_random, na.rm = TRUE),
  n_sig_fixed_and_random = sum(dependence_cells$robust_fixed_and_random, na.rm = TRUE),
  n_cells_reduced_fixed_effect = sum(dependence_cells$reduced_fixed_effect),
  n_cells_k_reduced_by_one_per_pss = sum(dependence_cells$n_pairs_one_per_pss < dependence_cells$n_pairs, na.rm = TRUE)
)
write_tbl(dependence_summary, "09_dependence_summary")

# One row in the figure gets the fixed effect, plus only the sensitivities that
# are a different estimand. Random effects and the single retained PGS never
# share a cell: the former needs ≥2 independent ΔAUCs, the latter is the case
# in which only one ΔAUC remains.
dependence_long <- bind_rows(
  dependence_cells %>%
    transmute(trait_label, trained_bucket, target_ancestry, n_pairs,
              delta_hat = delta_fixed, lo = lo_fixed, hi = hi_fixed, source = src_fixed),
  dependence_cells %>%
    filter(single_pgs_sensitivity) %>%
    transmute(trait_label, trained_bucket, target_ancestry, n_pairs = n_pairs_one_per_pss,
              delta_hat = delta_one_per_pss, lo = lo_one_per_pss, hi = hi_one_per_pss, source = src_single),
  dependence_cells %>%
    filter(reduced_fixed_effect) %>%
    transmute(trait_label, trained_bucket, target_ancestry, n_pairs = n_pairs_one_per_pss,
              delta_hat = delta_one_per_pss, lo = lo_one_per_pss, hi = hi_one_per_pss, source = src_reduced),
  dependence_cells %>%
    filter(re_applicable) %>%
    transmute(trait_label, trained_bucket, target_ancestry, n_pairs = n_delta_random,
              delta_hat = delta_random, lo = lo_random, hi = hi_random, source = src_random)
) %>%
  mutate(source = factor(source, levels = c(src_fixed, src_single, src_reduced, src_random)))

p_dep <- forest_compare_plot(
  dependence_long,
  # Outside the ancestry Okabe-Ito palette (no blue, orange, sky, green, yellow, pink, brown, grey, black).
  pal = c("#332288", "#E41A1C", "#6A3D9A", "#01665E"),
  shapes = c(16, 17, 18, 15),
  title = "Dependence sensitivities on the delta-method ΔAUC",
  subtitle = glue(
    "95% CI excludes zero for {dependence_summary$n_sig_fixed_corrected} IVW fixed effects. ",
    "Single-PGS sensitivity (k > 1, one target sample set): {dependence_summary$n_sig_single_pgs_sensitivity} of {dependence_summary$n_cells_single_pgs_sensitivity}.",
    "\nRandom effects, only where ≥2 independent ΔAUC remain: {dependence_summary$n_sig_random_dl} of {dependence_summary$n_cells_random_effects} exclude zero. ",
    "{dependence_summary$n_cells_single_delta} cells are one ΔAUC and are drawn once."
  ),
  caption = paste0(
    "IVW fixed effect: inverse-variance fixed effect of every PGS in the cell. ",
    "Single PGS: those scores share one target sample set; the point is the retained score (lowest PGS ID). Random effects are not fit to one ΔAUC.\n",
    "Fixed effect, one PGS per sample set: shown only when ≥2 scores remain and a duplicate is dropped. ",
    "Random effects (DL): DerSimonian-Laird on those independent ΔAUC values.\n",
    "k is the number of PGS in the IVW fixed effect. * = its 95% CI excludes zero."
  ),
  trait_levels = trait_levels_wrapped,
  show_k = FALSE
)

dependence_row_lab <- dependence_cells %>%
  mutate(
    trait_wrapped = factor(wrap_trait(trait_label), levels = trait_levels_wrapped),
    target_ancestry = factor(as.character(target_ancestry), levels = rev(anc_keep)),
    trained_bucket = factor(trained_bucket, levels = names(pal_bucket)),
    k_lab = paste0("k = ", n_pairs, ifelse(sig_fixed, " *", ""))
  )

p_dep <- p_dep +
  geom_text(
    data = dependence_row_lab,
    aes(x = x_lab_pos, y = target_ancestry, label = k_lab, fontface = ifelse(sig_fixed, "bold", "plain")),
    inherit.aes = FALSE, size = 2.5, hjust = 0, vjust = 0.5, colour = "grey20",
    show.legend = FALSE
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(
    panel.spacing.x = unit(1.15, "cm"),
    plot.margin = margin(6, 50, 6, 6)
  )
save_plot("09_dependence_forest", p_dep, width = 13.2, height = max(8.2, 0.44 * n_row_plot + 3.8), limitsize = FALSE)

# ===========================================================================
# 10) Covariate class
#     Every delta is one score and one class on both sides. Colour says
#     whether the comparison's scores are score-only, covariate-adjusted,
#     or both kinds.
# ===========================================================================

message(">> [10] Covariate class")

write_tbl(
  usable_keys %>%
    count(covariate_class, name = "n_evaluations") %>%
    mutate(share = n_evaluations / sum(n_evaluations)),
  "10_covariates_class_overall"
)

class_pairs <- pair_stage1_by_class(
  stage1_to_scatter(stage1_f, keep_traits),
  se_col = "se_delta_method"
) %>%
  mutate(
    delta = tgt_auc - eur_auc,
    delta_se = sqrt(tgt_se^2 + eur_se^2)
  ) %>%
  filter(is.finite(delta), is.finite(delta_se), delta_se > 0)

not_chosen <- class_pairs %>%
  anti_join(
    pick_one_class_per_score(class_pairs) %>%
      select(trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class),
    by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id", "covariate_class")
  ) %>%
  group_by(trait_label, trained_bucket, target_ancestry, pgs_id) %>%
  summarise(
    covariate_class_not_chosen = paste(covariate_class, collapse = " | "),
    delta_not_chosen = paste(sprintf("%.6f", delta), collapse = " | "),
    .groups = "drop"
  )

cov_txt <- usable_keys %>%
  group_by(pgs_id, trait_label, ancestry_display, covariate_class) %>%
  summarise(
    covariates_txt = paste(unique(coalesce(covariates, "")), collapse = " || "),
    .groups = "drop"
  )

pair_cov <- paired_dm %>%
  select(trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class, delta, delta_se, is_pooled, eur_k, tgt_k) %>%
  left_join(
    cov_txt %>% filter(ancestry_display == "European") %>%
      select(pgs_id, trait_label, covariate_class, eur_covariates_txt = covariates_txt),
    by = c("pgs_id", "trait_label", "covariate_class")
  ) %>%
  left_join(
    cov_txt %>%
      select(pgs_id, trait_label, ancestry_display, covariate_class, tgt_covariates_txt = covariates_txt),
    by = c("pgs_id", "trait_label", "covariate_class", "target_ancestry" = "ancestry_display")
  ) %>%
  left_join(not_chosen, by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id"))

if (anyDuplicated(pair_cov[c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")])) {
  stop("A score contributed more than one delta to a comparison.", call. = FALSE)
}
if (anyNA(pair_cov$covariate_class) || any(!pair_cov$covariate_class %in% covariate_class_levels)) {
  stop("A paired delta is missing its covariate class.", call. = FALSE)
}
n_class_mismatch <- 0L
chosen_vs_dm <- pair_cov %>%
  select(trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class, delta) %>%
  full_join(
    paired_dm %>% select(trait_label, trained_bucket, target_ancestry, pgs_id, class_dm = covariate_class, delta_dm = delta),
    by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")
  )
if (
  anyNA(chosen_vs_dm$covariate_class) || anyNA(chosen_vs_dm$class_dm) ||
    any(chosen_vs_dm$covariate_class != chosen_vs_dm$class_dm) ||
    any(abs(chosen_vs_dm$delta - chosen_vs_dm$delta_dm) > 1e-8)
) {
  stop("Paired deltas do not match the class chosen for each score.", call. = FALSE)
}

write_tbl(pair_cov, "10_covariates_pairs")

pre_input <- eval_df %>%
  filter(
    has_auc, has_ci,
    ancestry_display %in% evaluation_display_levels,
    is.finite(auc), is.finite(estimate_ci_lower), is.finite(estimate_ci_upper)
  )
pre_kept <- keep_one_evaluation(
  pre_input,
  n_col = "n_individuals",
  cases_col = "n_cases",
  controls_col = "n_controls",
  covariates_col = "covariates_txt",
  ancestry_col = "ancestry_display"
)
pre_paired <- pair_stage1_ignoring_class(
  stage1_to_scatter(
    i2_filter(pool_stage1(pre_kept, by_class = FALSE)) %>% attach_train_bucket(eval_df),
    keep_traits
  ),
  se_col = "se_delta_method"
)

rule_changes <- pair_cov %>%
  select(trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class, delta_new = delta, delta_se_new = delta_se) %>%
  full_join(
    pre_paired %>%
      select(trait_label, trained_bucket, target_ancestry, pgs_id, delta_pre_rule = delta, delta_se_pre_rule = delta_se),
    by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")
  ) %>%
  mutate(delta_diff = delta_new - delta_pre_rule) %>%
  filter(!is.finite(delta_diff) | abs(delta_diff) > 1e-6) %>%
  arrange(trait_label, trained_bucket, target_ancestry, pgs_id)

rule_changes <- rule_changes %>%
  mutate(change_type = case_when(
    is.finite(delta_pre_rule) & is.finite(delta_new) ~ "value",
    is.finite(delta_new) & !is.finite(delta_pre_rule) ~ "admitted",
    is.finite(delta_pre_rule) & !is.finite(delta_new) ~ "dropped",
    TRUE ~ "other"
  ))

expected_value_changes <- tibble(
  trait_label = c(
    rep("Coronary heart disease (incident and prevalent)", 4),
    "Type 2 diabetes"
  ),
  trained_bucket = "Train: Multi incl. EUR",
  target_ancestry = c(
    "African", "Hispanic or Latin American",
    "African", "Hispanic or Latin American",
    "Middle Eastern or North African"
  ),
  pgs_id = c("PGS000013", "PGS000013", "PGS000018", "PGS000018", "PGS000804")
)
got_value_changes <- rule_changes %>%
  filter(change_type == "value") %>%
  select(trait_label, trained_bucket, target_ancestry, pgs_id) %>%
  arrange(trait_label, trained_bucket, target_ancestry, pgs_id)
exp_value_changes <- expected_value_changes %>%
  arrange(trait_label, trained_bucket, target_ancestry, pgs_id)
if (!isTRUE(all.equal(got_value_changes, exp_value_changes, check.attributes = FALSE))) {
  print(as.data.frame(rule_changes), row.names = FALSE)
  stop("Value changes are not the five expected scores.", call. = FALSE)
}
if (any(rule_changes$change_type == "dropped")) {
  print(as.data.frame(rule_changes %>% filter(change_type == "dropped")), row.names = FALSE)
  stop("The class rule dropped a pair.", call. = FALSE)
}
write_tbl(rule_changes, "10_covariates_rule_changes")

cov_colour_levels <- covariate_colour_levels
pal_cov_colour <- c(
  "Uses covariates" = "#332288",
  "Does not use covariates" = "#882255",
  "Both kinds of score" = "#999933"
)

cov_colour_tbl <- pair_cov %>%
  group_by(trait_label, trained_bucket, target_ancestry) %>%
  summarise(
    n_pairs = n(),
    n_uses = sum(covariate_class != "None (PRS only)"),
    n_none = sum(covariate_class == "None (PRS only)"),
    covariate_colour = comparison_covariate_colour(covariate_class),
    .groups = "drop"
  ) %>%
  mutate(covariate_colour = factor(covariate_colour, levels = cov_colour_levels))

write_tbl(cov_colour_tbl, "10_covariates_comparison_colour")
old_matched <- file.path(dir_tables, "10_covariates_cells_matched_sensitivity.csv")
if (file.exists(old_matched)) unlink(old_matched)

# 120 pairs become 127 because seven scores were previously one high-I2 average
# of two covariate classes. Each class is now its own side, so those scores pair.
if (nrow(pair_cov) != 127L || nrow(cov_colour_tbl) != 32L) {
  stop(
    glue("Expected 127 pairs and 32 comparisons; got {nrow(pair_cov)} pairs and {nrow(cov_colour_tbl)} comparisons."),
    call. = FALSE
  )
}

cov_summary <- tibble(
  n_pairs = nrow(pair_cov),
  n_pairs_score_only = sum(pair_cov$covariate_class == "None (PRS only)"),
  n_pairs_age_sex_pcs = sum(pair_cov$covariate_class == "Age/sex/PCs/technical"),
  n_pairs_clinical = sum(pair_cov$covariate_class == "Clinical or other covariates"),
  n_pairs_class_mismatch = n_class_mismatch,
  n_comparisons = nrow(cov_colour_tbl),
  n_comparisons_uses_covariates = sum(cov_colour_tbl$covariate_colour == "Uses covariates"),
  n_comparisons_score_only = sum(cov_colour_tbl$covariate_colour == "Does not use covariates"),
  n_comparisons_both_kinds = sum(cov_colour_tbl$covariate_colour == "Both kinds of score"),
  n_pairs_rule_changed = nrow(rule_changes)
)
write_tbl(cov_summary, "10_covariates_summary")

cov_plot_df <- delta_dm %>%
  transmute(trait_label, trained_bucket, target_ancestry, n_pairs, delta_hat, lo, hi) %>%
  left_join(
    cov_colour_tbl %>%
      mutate(covariate_colour = as.character(covariate_colour)) %>%
      select(trait_label, trained_bucket, target_ancestry, covariate_colour),
    by = c("trait_label", "trained_bucket", "target_ancestry")
  ) %>%
  mutate(
    trait_wrapped = factor(wrap_trait(trait_label), levels = trait_levels_wrapped),
    target_ancestry = factor(as.character(target_ancestry), levels = rev(anc_keep)),
    trained_bucket = factor(trained_bucket, levels = names(pal_bucket)),
    covariate_colour = factor(covariate_colour, levels = cov_colour_levels),
    sig = (lo > 0) | (hi < 0)
  )

n_cov_uses <- sum(cov_colour_tbl$covariate_colour == "Uses covariates")
n_cov_none <- sum(cov_colour_tbl$covariate_colour == "Does not use covariates")
n_cov_both <- sum(cov_colour_tbl$covariate_colour == "Both kinds of score")

p_cov <- ggplot(
  cov_plot_df,
  aes(x = delta_hat, y = target_ancestry, colour = covariate_colour)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.5, orientation = "y") +
  geom_point(size = 2.3, shape = 16) +
  scale_colour_manual(values = pal_cov_colour, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = x_breaks_forest, labels = x_labels_forest) +
  coord_cartesian(xlim = x_lim_forest, clip = "off") +
  facet_grid(trait_wrapped ~ trained_bucket, scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title = "ΔAUC by covariate content of the scores in each comparison",
    subtitle = str_wrap(glue(
      "{n_cov_uses} comparisons use covariates on every score. ",
      "{n_cov_none} use none. ",
      "{n_cov_both} include both kinds of score. ",
      "Every delta is one covariate class on the European side and on the non-European side. ",
      "The accompanying figure shows each score in the comparisons that include both kinds."
    ), width = 120),
    x = "ΔAUC (non-European − European) with 95% CI",
    y = NULL,
    caption = str_wrap(paste0(
      "Each point is the average of all scores in the comparison. ",
      "Colour is the covariate class of those scores: covariates on every score, ",
      "no covariates on any score, or both kinds of score. ",
      "A class present on one side only is not part of a delta. ",
      "When one score has more than one class on both sides, the score-only delta is kept, ",
      "then age, sex, and principal components, then clinical covariates. ",
      "Age, sex, and principal components count as covariates."
    ), width = 150)
  ) +
  theme_reviewer(base_size = 9) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8, face = "bold", lineheight = 0.85),
    strip.text.x = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey95", colour = NA),
    panel.spacing.y = unit(0.35, "lines"),
    panel.spacing.x = unit(1.6, "lines"),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(5.5, 20, 5.5, 5.5),
    plot.caption = element_text(hjust = 0, size = 7.5),
    plot.title = element_text(size = 10.5)
  )
save_plot("10_covariates_forest", p_cov, width = 12.5, height = max(7.2, 0.32 * n_row_plot + 2.8), limitsize = FALSE)

training_label <- c(
  "Train: European-only" = "European-only training",
  "Train: Multi incl. EUR" = "Multi-ancestry training"
)
focus_class_label <- c(
  "None (PRS only)" = "Does not use covariates",
  "Age/sex/PCs/technical" = "Age, sex, and principal components",
  "Clinical or other covariates" = "Clinical covariates"
)
pal_focus_class <- c(
  "Does not use covariates" = "#882255",
  "Age, sex, and principal components" = "#332288",
  "Clinical covariates" = "#E41A1C"
)

# One row per covariate class. Several scores become their IVW.
# A class with one score stays that score.
focus_class_estimate <- function(d) {
  if (nrow(d) == 1L) {
    return(tibble(
      estimate = "One score",
      row_label = d$pgs_id[[1]],
      n_scores = 1L,
      delta = d$delta[[1]],
      delta_se = d$delta_se[[1]]
    ))
  }
  averaged <- pool_delta(d$delta, d$delta_se, model = "fixed")
  tibble(
    estimate = "IVW",
    row_label = paste0("IVW of ", nrow(d), " scores"),
    n_scores = nrow(d),
    delta = averaged$delta_hat,
    delta_se = averaged$se_hat
  )
}

focus_deltas <- pair_cov %>%
  inner_join(
    cov_colour_tbl %>%
      filter(covariate_colour == "Both kinds of score") %>%
      select(trait_label, trained_bucket, target_ancestry),
    by = c("trait_label", "trained_bucket", "target_ancestry")
  ) %>%
  group_by(trait_label, trained_bucket, target_ancestry, covariate_class) %>%
  group_modify(~ focus_class_estimate(.x)) %>%
  ungroup() %>%
  mutate(
    class_label = unname(focus_class_label[covariate_class]),
    class_label = factor(class_label, levels = names(pal_focus_class)),
    estimate = factor(estimate, levels = c("IVW", "One score")),
    training = unname(training_label[as.character(trained_bucket)]),
    strip = paste0(wrap_trait(trait_label), "\n", training, "\n", target_ancestry),
    lo = delta - 1.96 * delta_se,
    hi = delta + 1.96 * delta_se
  )

strip_levels <- focus_deltas %>%
  distinct(trait_label, trained_bucket, target_ancestry, strip) %>%
  mutate(
    trait_label = factor(trait_label, levels = trait_order),
    trained_bucket = factor(trained_bucket, levels = names(pal_bucket)),
    target_ancestry = factor(target_ancestry, levels = anc_keep)
  ) %>%
  arrange(trait_label, trained_bucket, target_ancestry) %>%
  pull(strip)

focus_deltas <- focus_deltas %>%
  mutate(strip = factor(strip, levels = strip_levels)) %>%
  arrange(strip, delta) %>%
  mutate(row_id = factor(
    paste(row_label, strip, sep = " | "),
    levels = unique(paste(row_label, strip, sep = " | "))
  ))

if (any(focus_deltas %>% count(strip) %>% pull(n) < 2L)) {
  stop("A comparison on the covariate figure has only one estimate.", call. = FALSE)
}
focus_check <- focus_deltas %>%
  mutate(strip_key = paste(trait_label, trained_bucket, target_ancestry, sep = " | "))
stopifnot(
  any(focus_check$strip_key == "Breast cancer | Train: European-only | South Asian" &
        focus_check$covariate_class == "None (PRS only)" & focus_check$n_scores == 10L),
  any(focus_check$row_label == "PGS001336"),
  any(focus_check$strip_key == "Type 2 diabetes | Train: European-only | South Asian" &
        focus_check$covariate_class == "Age/sex/PCs/technical" & focus_check$n_scores == 12L),
  any(focus_check$row_label == "PGS004106"),
  sum(focus_check$strip_key == "Type 2 diabetes | Train: Multi incl. EUR | African" &
        focus_check$estimate == "One score") == 2L,
  any(focus_check$strip_key == "Type 2 diabetes | Train: Multi incl. EUR | Hispanic or Latin American" &
        focus_check$covariate_class == "Age/sex/PCs/technical" & focus_check$n_scores == 4L),
  all(c("PGS000804", "PGS005242") %in% focus_check$row_label)
)

write_tbl(
  focus_deltas %>%
    transmute(
      trait_label, trained_bucket, target_ancestry, covariate_class, class_label,
      estimate, row_label, n_scores, delta, delta_se, lo, hi
    ),
  "10b_covariates_class_estimates"
)

focus_span <- range(c(focus_deltas$lo, focus_deltas$hi), na.rm = TRUE)
focus_pad <- max(0.01, 0.08 * diff(focus_span))

p_focus <- ggplot(focus_deltas, aes(x = delta, y = row_id, colour = class_label, shape = estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.55, orientation = "y") +
  geom_point(size = 3) +
  scale_colour_manual(values = pal_focus_class, name = NULL, drop = TRUE) +
  scale_shape_manual(values = c("IVW" = 18, "One score" = 16), name = NULL) +
  scale_y_discrete(labels = function(x) sub(" \\|.*", "", x)) +
  scale_x_continuous(labels = label_number(accuracy = 0.02, style_negative = "minus")) +
  coord_cartesian(xlim = c(focus_span[1] - focus_pad, focus_span[2] + focus_pad)) +
  facet_grid(strip ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title = "IVW of the larger covariate class, and the single score of each other class",
    subtitle = str_wrap(paste0(
      "Comparisons that include more than one covariate class. ",
      "A class with several scores is their IVW. ",
      "A class with one score stays that score."
    ), width = 105),
    x = "ΔAUC (non-European − European) with 95% CI",
    y = NULL,
    caption = str_wrap(paste0(
      "The European side and the non-European side of every estimate use the same covariate class. ",
      "Does not use covariates: None (PRS only). ",
      "Age, sex, and principal components are separate from clinical covariates. ",
      "Breast cancer, South Asian participants: IVW of 10 scores with no covariates, and PGS001336. ",
      "Type 2 diabetes, South Asian participants, European-only training: IVW of 12 scores with age, sex, and principal components, and PGS004106. ",
      "African participants: PGS000804 and PGS004859, one score each. ",
      "Hispanic or Latin American participants: IVW of 4 scores with age, sex, and principal components, PGS000804, and PGS005242."
    ), width = 120)
  ) +
  theme_reviewer(base_size = 9) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8, face = "bold", lineheight = 0.9),
    strip.background = element_rect(fill = "grey95", colour = NA),
    panel.spacing.y = unit(0.7, "lines"),
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0, size = 7.5),
    plot.title = element_text(size = 10.5),
    legend.box = "vertical"
  )
save_plot(
  "10b_covariates_by_phenotype",
  p_focus,
  width = 8.8,
  height = max(6.4, 0.42 * nrow(focus_deltas) + 1.15 * nlevels(focus_deltas$strip) + 2.4),
  limitsize = FALSE
)
old_chd <- file.path(dir_plots, c("10b_covariates_chd.png", "10b_covariates_chd.pdf"))
unlink(old_chd[file.exists(old_chd)])


# ===========================================================================
# 11) Training ancestry
#     Parse GWAS and development ancestry distributions per PGS; flag whether
#     the target ancestry was represented; cross-check the training bucket.
# ===========================================================================

message(">> [11] Training ancestry")

map_dev_label <- function(lbl) {
  lbl <- trimws(lbl)
  case_when(
    grepl("Multi-ancestry \\(including European\\)", lbl, ignore.case = TRUE) ~ "Multi-ancestry including European",
    grepl("Multi-ancestry \\(excluding European\\)", lbl, ignore.case = TRUE) ~ "Multi-ancestry excluding European",
    grepl("^European", lbl, ignore.case = TRUE) ~ "European",
    grepl("African|Sub-Saharan|Afro", lbl, ignore.case = TRUE) ~ "African",
    grepl("East Asian", lbl, ignore.case = TRUE) ~ "East Asian",
    grepl("South Asian", lbl, ignore.case = TRUE) ~ "South Asian",
    grepl("Hispanic|Latin", lbl, ignore.case = TRUE) ~ "Hispanic or Latin American",
    grepl("Middle Eastern|North African", lbl, ignore.case = TRUE) ~ "Middle Eastern or North African",
    grepl("Not reported|NR", lbl, ignore.case = TRUE) ~ "Not reported",
    TRUE ~ "Other/Mixed"
  )
}

parse_anc_dist <- function(txt) {
  if (is.na(txt) || txt == "") return(tibble(group = character(), pct = numeric()))
  parts <- str_split(txt, "\\|")[[1]]
  tibble(
    group = map_dev_label(str_replace(parts, ":[^:]*$", "")),
    pct = suppressWarnings(as.numeric(str_extract(parts, "[0-9.]+$")))
  ) %>% group_by(group) %>% summarise(pct = sum(pct, na.rm = TRUE), .groups = "drop")
}

pgs_in_pairs <- sort(unique(paired_dm$pgs_id))

pgs_dev <- scores_meta %>%
  filter(pgs_id %in% pgs_in_pairs) %>%
  mutate(
    gwas = map(anc_gwas_txt, parse_anc_dist),
    dev  = map(anc_dev_txt, parse_anc_dist)
  )

pgs_dev_long <- bind_rows(
  pgs_dev %>% select(pgs_id, dist = gwas) %>% unnest(dist) %>% mutate(stage = "GWAS (variant associations)"),
  pgs_dev %>% select(pgs_id, dist = dev) %>% unnest(dist) %>% mutate(stage = "Score development/training")
)
write_tbl(pgs_dev_long, "11_training_ancestry_pgs_long")

pct_of <- function(pgs, stage_tbl_col, grp) {
  tbl <- pgs_dev[[stage_tbl_col]][[match(pgs, pgs_dev$pgs_id)]]
  if (is.null(tbl) || !nrow(tbl)) return(NA_real_)
  v <- tbl$pct[tbl$group == grp]
  if (length(v)) v[1] else 0
}

pair_train <- paired_dm %>%
  select(trait_label, trained_bucket, target_ancestry, pgs_id, delta, delta_se) %>%
  rowwise() %>%
  mutate(
    gwas_pct_european = pct_of(pgs_id, "gwas", "European"),
    gwas_pct_target = pct_of(pgs_id, "gwas", target_ancestry),
    gwas_pct_multi_incl_eur = pct_of(pgs_id, "gwas", "Multi-ancestry including European"),
    dev_pct_european = pct_of(pgs_id, "dev", "European"),
    dev_pct_target = pct_of(pgs_id, "dev", target_ancestry),
    dev_pct_multi_incl_eur = pct_of(pgs_id, "dev", "Multi-ancestry including European"),
    dev_reported = !is.na(dev_pct_european) | !is.na(dev_pct_target) | !is.na(dev_pct_multi_incl_eur)
  ) %>%
  ungroup() %>%
  mutate(
    target_in_gwas = coalesce(gwas_pct_target, 0) > 0,
    target_in_dev = coalesce(dev_pct_target, 0) > 0,
    multi_incl_eur_in_gwas_or_dev = coalesce(gwas_pct_multi_incl_eur, 0) > 0 | coalesce(dev_pct_multi_incl_eur, 0) > 0,
    target_possibly_represented = target_in_gwas | target_in_dev | multi_incl_eur_in_gwas_or_dev,
    bucket_consistent = case_when(
      grepl("European-only", trained_bucket) ~ !target_possibly_represented,
      grepl("Multi incl", trained_bucket) ~ TRUE,
      TRUE ~ NA
    )
  ) %>%
  left_join(scores_meta %>% select(pgs_id, pgs_name, dev_method, n_variants, anc_gwas_txt, anc_dev_txt), by = "pgs_id")

write_tbl(pair_train, "11_training_ancestry_pairs")

train_cells <- pair_train %>%
  group_by(trait_label, trained_bucket, target_ancestry) %>%
  summarise(
    n_pairs = n(),
    n_pairs_target_in_gwas = sum(target_in_gwas),
    n_pairs_target_in_dev = sum(target_in_dev),
    n_pairs_multi_incl_eur_ambiguous = sum(multi_incl_eur_in_gwas_or_dev & !target_in_gwas & !target_in_dev),
    n_pairs_target_possibly_represented = sum(target_possibly_represented),
    n_pairs_bucket_inconsistent = sum(!coalesce(bucket_consistent, TRUE)),
    median_gwas_pct_european = median(gwas_pct_european, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(delta_dm %>% select(trait_label, trained_bucket, target_ancestry, delta_hat, lo, hi, ci_excludes_zero),
            by = c("trait_label", "trained_bucket", "target_ancestry")) %>%
  mutate(strict_portability_cell = n_pairs_target_possibly_represented == 0)

write_tbl(train_cells, "11_training_ancestry_cells")

train_summary <- tibble(
  n_pairs = nrow(pair_train),
  n_pairs_target_in_gwas = sum(pair_train$target_in_gwas),
  n_pairs_target_in_dev = sum(pair_train$target_in_dev),
  n_pairs_multi_incl_eur_ambiguous = sum(pair_train$multi_incl_eur_in_gwas_or_dev & !pair_train$target_in_gwas & !pair_train$target_in_dev),
  n_pairs_eur_only_bucket_but_target_represented = sum(grepl("European-only", pair_train$trained_bucket) & pair_train$target_possibly_represented),
  n_cells = nrow(train_cells),
  n_cells_strict_portability = sum(train_cells$strict_portability_cell),
  n_cells_strict_and_significant = sum(train_cells$strict_portability_cell & train_cells$ci_excludes_zero, na.rm = TRUE)
)
write_tbl(train_summary, "11_training_ancestry_summary")

# ===========================================================================
# 12) Complementary metric: OR / HR per SD
#     Same two-stage design on log(OR) (HR used when OR absent).
#     The Catalog does not guarantee per-SD scaling; only pairs where both
#     sides report the same metric are used.
# ===========================================================================

message(">> [12] OR/HR per SD and case fractions")

parse_ratio <- function(txt) {
  m <- str_match(coalesce(txt, ""), "^\\s*([0-9.]+)\\s*\\[\\s*([0-9.]+)\\s*,\\s*([0-9.]+)\\s*\\]")
  tibble(est = as.numeric(m[, 2]), lo = as.numeric(m[, 3]), hi = as.numeric(m[, 4]))
}

or_df <- eval_df %>%
  select(performance_id, pgs_id, trait_label, ancestry_display, train_bucket, sampleset_id, or_txt, hr_txt) %>%
  mutate(
    metric = case_when(!is.na(or_txt) & or_txt != "" ~ "OR", !is.na(hr_txt) & hr_txt != "" ~ "HR", TRUE ~ NA_character_),
    ratio_txt = ifelse(metric == "OR", or_txt, hr_txt)
  ) %>%
  filter(!is.na(metric)) %>%
  bind_cols(parse_ratio(.$ratio_txt)) %>%
  filter(is.finite(est), is.finite(lo), is.finite(hi), est > 0, lo > 0, hi > lo) %>%
  mutate(log_est = log(est), log_se = (log(hi) - log(lo)) / (2 * 1.96)) %>%
  filter(is.finite(log_se), log_se > 0, ancestry_display != "Not reported")

stage1_or <- or_df %>%
  group_by(pgs_id, trait_label, ancestry_display, metric) %>%
  summarise(pooled = list(ivw_pool_logit(log_est, log_se)), .groups = "drop") %>%
  tidyr::unnest_wider(pooled) %>%
  rename(log_ratio = eta, log_se = se) %>%
  mutate(flag_high_I2 = !is.na(I2) & k_eval >= 2 & I2 > 80) %>%
  filter(!flag_high_I2) %>%
  select(-flag_high_I2) %>%
  attach_train_bucket(eval_df) %>%
  filter(trait_label %in% keep_traits,
         grepl(paste(keep_bucket_pattern, collapse = "|"), trained_bucket, ignore.case = TRUE))

or_eur <- stage1_or %>% filter(ancestry_display == "European") %>%
  transmute(trait_label, pgs_id, trained_bucket, metric, eur_log = log_ratio, eur_se = log_se, eur_k = k_eval)
or_tgt <- stage1_or %>% filter(ancestry_display != "European") %>%
  transmute(trait_label, pgs_id, trained_bucket, metric, target_ancestry = ancestry_display, tgt_log = log_ratio, tgt_se = log_se, tgt_k = k_eval)

or_pairs <- or_tgt %>%
  inner_join(or_eur, by = c("trait_label", "pgs_id", "trained_bucket", "metric")) %>%
  mutate(is_pooled = eur_k >= 2 & tgt_k >= 2, delta = tgt_log - eur_log, delta_se = sqrt(tgt_se^2 + eur_se^2)) %>%
  filter(is.finite(delta), is.finite(delta_se), delta_se > 0)

write_tbl(or_pairs, "12_or_per_sd_pairs")

or_cells <- pool_stage2(or_pairs, model = "fixed") %>%
  rename(delta_log_ratio = delta_hat, se_log_ratio = se_hat, lo_log_ratio = lo, hi_log_ratio = hi, sig_log_ratio = ci_excludes_zero) %>%
  mutate(ratio_of_ratios = exp(delta_log_ratio), ratio_lo = exp(lo_log_ratio), ratio_hi = exp(hi_log_ratio)) %>%
  left_join(or_pairs %>% group_by(trait_label, trained_bucket, target_ancestry) %>% summarise(metrics_used = paste(sort(unique(metric)), collapse = "+"), .groups = "drop"),
            by = c("trait_label", "trained_bucket", "target_ancestry")) %>%
  full_join(delta_dm %>% transmute(trait_label, trained_bucket, target_ancestry, n_pairs_auc = n_pairs, delta_auc = delta_hat, lo_auc = lo, hi_auc = hi, sig_auc = ci_excludes_zero),
            by = c("trait_label", "trained_bucket", "target_ancestry")) %>%
  mutate(
    in_figure1c = !is.na(delta_auc),
    direction_agrees = sign(delta_log_ratio) == sign(delta_auc),
    both_significant_same_direction = sig_log_ratio & sig_auc & direction_agrees,
    sig_opposite_direction = sig_log_ratio & sig_auc & !direction_agrees
  ) %>%
  arrange(desc(in_figure1c), trait_label, trained_bucket, target_ancestry)

write_tbl(or_cells, "12_or_per_sd_cells")

or_cells_1c <- or_cells %>% filter(in_figure1c)

or_summary <- tibble(
  n_ppm_with_or_or_hr = nrow(or_df),
  n_stage1_or_cells = nrow(stage1_or),
  n_or_pairs = nrow(or_pairs),
  n_figure1c_cells = nrow(delta_dm),
  n_figure1c_cells_with_or_estimate = sum(!is.na(or_cells_1c$delta_log_ratio)),
  n_cells_direction_agrees = sum(or_cells_1c$direction_agrees, na.rm = TRUE),
  n_cells_direction_disagrees = sum(!or_cells_1c$direction_agrees, na.rm = TRUE),
  n_cells_both_sig_same_direction = sum(or_cells_1c$both_significant_same_direction, na.rm = TRUE),
  n_cells_both_sig_opposite_direction = sum(or_cells_1c$sig_opposite_direction, na.rm = TRUE),
  n_cells_sig_log_ratio = sum(or_cells_1c$sig_log_ratio, na.rm = TRUE),
  n_extra_or_cells_not_in_figure1c = sum(!or_cells$in_figure1c)
)
write_tbl(or_summary, "12_or_per_sd_summary")

or_plot_df <- or_cells %>% filter(!is.na(delta_log_ratio), !is.na(delta_auc)) %>%
  mutate(target_ancestry = factor(as.character(target_ancestry), levels = anc_keep))
if (nrow(or_plot_df) > 0) {
  p_or <- ggplot(or_plot_df, aes(x = delta_auc, y = delta_log_ratio, colour = target_ancestry, shape = trained_bucket)) +
    geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
    geom_errorbar(aes(ymin = lo_log_ratio, ymax = hi_log_ratio), width = 0, alpha = 0.5) +
    geom_errorbar(aes(xmin = lo_auc, xmax = hi_auc), width = 0, alpha = 0.5, orientation = "y") +
    geom_point(size = 2.5) +
    scale_colour_manual(values = pal_ancestry, name = "Target ancestry") +
    scale_shape_manual(values = c(16, 17), name = "Training") +
    labs(
      title = "ΔAUC vs Δlog(OR or HR per SD), same cells",
      subtitle = glue(
        "{or_summary$n_figure1c_cells_with_or_estimate} of {or_summary$n_figure1c_cells} Figure 1C cells have a paired OR/HR estimate; ",
        "direction agrees in {or_summary$n_cells_direction_agrees}; both significant and opposite in {or_summary$n_cells_both_sig_opposite_direction}."
      ),
      x = "ΔAUC (target − European), corrected 95% CI", y = "Δlog(ratio) (target − European), 95% CI",
      caption = "Catalog OR/HR fields parsed as 'est [lo,hi]'; per-SD scaling not verifiable from the Catalog; only pairs with the same metric on both sides are used."
    ) +
    guides(colour = guide_legend(nrow = 2)) +
    theme_reviewer()
  save_plot("12_or_per_sd_scatter", p_or, width = 8.5, height = 7)
}

# Case fractions per cell side
side_pss <- bind_rows(
  paired_dm %>% transmute(trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class, side = "European", ancestry = "European"),
  paired_dm %>% transmute(trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class, side = "Target", ancestry = target_ancestry)
) %>%
  # many-to-many is intended: each pair side expands to the evaluations in its class
  inner_join(
    usable_keys %>% select(pgs_id, trait_label, ancestry_display, covariate_class, sampleset_id, performance_id, ci_source),
    by = c("pgs_id", "trait_label", "ancestry" = "ancestry_display", "covariate_class"),
    relationship = "many-to-many"
  ) %>%
  left_join(pss, by = "sampleset_id")

case_fraction_cells <- side_pss %>%
  distinct(trait_label, trained_bucket, target_ancestry, side, sampleset_id, n_cases, n_controls, n_individuals) %>%
  group_by(trait_label, trained_bucket, target_ancestry, side) %>%
  summarise(
    n_sample_sets = n_distinct(sampleset_id),
    n_sample_sets_with_counts = sum(is.finite(n_cases) & is.finite(n_controls)),
    total_cases = sum(n_cases, na.rm = TRUE),
    total_controls = sum(n_controls, na.rm = TRUE),
    total_individuals = sum(n_individuals, na.rm = TRUE),
    min_cases_in_a_sample_set = suppressWarnings(min(n_cases, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    case_fraction = ifelse(total_cases + total_controls > 0, total_cases / (total_cases + total_controls), NA_real_),
    min_cases_in_a_sample_set = ifelse(is.infinite(min_cases_in_a_sample_set), NA_real_, min_cases_in_a_sample_set)
  ) %>%
  pivot_wider(id_cols = c(trait_label, trained_bucket, target_ancestry),
              names_from = side, values_from = c(n_sample_sets, n_sample_sets_with_counts, total_cases, total_controls, total_individuals, min_cases_in_a_sample_set, case_fraction),
              names_glue = "{tolower(side)}_{.value}")

write_tbl(case_fraction_cells, "12_case_fraction_by_cell")

# ===========================================================================
# 13) Identifiers behind every cell
#     PGS / PPM / PSS, CI source, covariates, cohort text, per side.
# ===========================================================================

message(">> [13] Cell identifiers")

cell_identifiers <- side_pss %>%
  left_join(ppm_meta %>% select(performance_id, covariates_txt, or_txt, hr_txt), by = "performance_id") %>%
  left_join(eval_df %>% select(performance_id, auc, estimate_ci_lower, estimate_ci_upper), by = "performance_id") %>%
  left_join(scores_meta %>% select(pgs_id, pgs_name, dev_method, n_variants), by = "pgs_id") %>%
  transmute(
    catalog_snapshot = stamp,
    trait_label, trained_bucket, target_ancestry, side, evaluation_ancestry = ancestry,
    pgs_id, pgs_name, dev_method, n_variants,
    performance_id, sampleset_id, auc, estimate_ci_lower, estimate_ci_upper, ci_source,
    n_individuals, n_cases, n_controls, pct_male, covariates_txt, cohort_txt, or_txt, hr_txt,
    catalog_all_of_us = grepl("all of us|allofus|\\baou\\b", coalesce(cohort_txt, ""), ignore.case = TRUE)
  ) %>%
  arrange(trait_label, trained_bucket, target_ancestry, pgs_id, side, performance_id)

write_tbl(cell_identifiers, "13_cell_identifiers_long")

cell_identifiers_summary <- cell_identifiers %>%
  group_by(trait_label, trained_bucket, target_ancestry) %>%
  summarise(
    n_pgs = n_distinct(pgs_id),
    pgs_ids = paste(sort(unique(pgs_id)), collapse = ";"),
    n_ppm_eur = n_distinct(performance_id[side == "European"]),
    n_ppm_target = n_distinct(performance_id[side == "Target"]),
    n_pss_eur = n_distinct(sampleset_id[side == "European"]),
    n_pss_target = n_distinct(sampleset_id[side == "Target"]),
    n_ppm_reported_ci = sum(ci_source == "reported_native_or_parsed"),
    n_ppm_recovered_ci = sum(ci_source == "inferred_hm"),
    ppm_ids_eur = paste(sort(unique(performance_id[side == "European"])), collapse = ";"),
    ppm_ids_target = paste(sort(unique(performance_id[side == "Target"])), collapse = ";"),
    pss_ids_eur = paste(sort(unique(sampleset_id[side == "European"])), collapse = ";"),
    pss_ids_target = paste(sort(unique(sampleset_id[side == "Target"])), collapse = ";"),
    .groups = "drop"
  )
write_tbl(cell_identifiers_summary, "13_cell_identifiers_by_cell")

# ===========================================================================
# 14) I2 range check
# ===========================================================================

message(">> [14] I2 range check")

i2_check <- bind_rows(
  stage1   %>% transmute(table = "stage1_pooled_cells", I2, k_eval),
  stage1_f %>% transmute(table = "stage1_pooled_cells_I2filtered", I2, k_eval),
  stage1_or %>% transmute(table = "stage1_or_hr (this script)", I2, k_eval)
) %>%
  group_by(table) %>%
  summarise(
    n_rows = n(),
    n_I2_na = sum(is.na(I2)),
    n_I2_k1_na_expected = sum(is.na(I2) & k_eval == 1),
    n_I2_below_0 = sum(I2 < 0, na.rm = TRUE),
    n_I2_above_100 = sum(I2 > 100, na.rm = TRUE),
    n_I2_above_80 = sum(I2 > 80, na.rm = TRUE),
    I2_min = suppressWarnings(min(I2, na.rm = TRUE)),
    I2_max = suppressWarnings(max(I2, na.rm = TRUE)),
    unit = "percentage (0-100)",
    .groups = "drop"
  ) %>%
  mutate(
    in_range = n_I2_below_0 == 0 & n_I2_above_100 == 0,
    note = "I2 is stored as a percentage on 0-100."
  )
write_tbl(i2_check, "14_i2_range_check")

# ===========================================================================
# 15) Sex composition of the evaluation sample sets behind Figure 1C,
#     per cell side; sex-specific traits flagged.
# ===========================================================================

message(">> [15] Sex composition")

sex_cells <- side_pss %>%
  distinct(trait_label, trained_bucket, target_ancestry, side, sampleset_id, pct_male) %>%
  group_by(trait_label, trained_bucket, target_ancestry, side) %>%
  summarise(
    n_sample_sets = n_distinct(sampleset_id),
    n_sample_sets_with_sex = sum(is.finite(pct_male)),
    median_pct_male = median(pct_male, na.rm = TRUE),
    min_pct_male = suppressWarnings(min(pct_male, na.rm = TRUE)),
    max_pct_male = suppressWarnings(max(pct_male, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    min_pct_male = ifelse(is.infinite(min_pct_male), NA_real_, min_pct_male),
    max_pct_male = ifelse(is.infinite(max_pct_male), NA_real_, max_pct_male),
    sex_specific_trait = grepl("breast|ovarian|prostate|endometri|cervi", trait_label, ignore.case = TRUE),
    single_sex_sample_sets = is.finite(median_pct_male) & (median_pct_male == 0 | median_pct_male == 100)
  ) %>%
  arrange(trait_label, trained_bucket, target_ancestry, side)

write_tbl(sex_cells, "15_sex_composition_by_cell")

sex_summary <- tibble(
  n_cell_sides = nrow(sex_cells),
  n_cell_sides_with_sex_info = sum(sex_cells$n_sample_sets_with_sex > 0),
  n_cells_sex_specific_trait = n_distinct(paste(sex_cells$trait_label, sex_cells$trained_bucket, sex_cells$target_ancestry)[sex_cells$sex_specific_trait]),
  share_usable_evaluations_with_pct_male = mean(is.finite(usable$pct_male)),
  note = "Catalog AUCs are reported for pooled-sex samples except sex-specific traits; no sex-stratified meta-analysis is possible from Catalog records."
)
write_tbl(sex_summary, "15_sex_composition_summary")

# ===========================================================================
# 16) k = 1 gate
#     The figure keeps a trait if any of its cells has at least two PGS.
# ===========================================================================

message(">> [16] k = 1 gate check")

k1_check <- delta_dm %>%
  transmute(trait_label, trained_bucket, target_ancestry, n_pairs, delta_hat, lo, hi, ci_excludes_zero) %>%
  mutate(
    k1_cell = n_pairs == 1,
    trait_has_any_cell_k_ge_2 = trait_label %in% (delta_dm %>% filter(n_pairs >= 2) %>% pull(trait_label) %>% unique())
  ) %>%
  arrange(desc(k1_cell), trait_label, trained_bucket, target_ancestry)
write_tbl(k1_check, "16_k1_gate_check")

k1_summary <- tibble(
  n_cells = nrow(k1_check),
  n_cells_k1 = sum(k1_check$k1_cell),
  n_cells_k1_significant_corrected = sum(k1_check$k1_cell & k1_check$ci_excludes_zero),
  n_traits = n_distinct(k1_check$trait_label),
  n_traits_all_cells_k1 = n_distinct(k1_check$trait_label[!k1_check$trait_has_any_cell_k_ge_2]),
  trait_level_gate_leaves_k1_cells = n_cells_k1 > 0,
  note = "The figure keeps a trait if any of its cells has at least two PGS; single-PGS cells can remain."
)
write_tbl(k1_summary, "16_k1_gate_summary")

# ===========================================================================
# 17) Can Catalog R² be meta-analyzed?
#     The question is not whether an AUROC row also prints an R².
#     A meta-analysis needs the R² itself, an interval, and more than
#     one score in the same trait. Counts use the bulk performance file.
#     "PGS-specific R²" is the Catalog field "PGS R2 (no covariates)".
#     Other R² intervals are counted separately.
#     Field definitions match tests/q_r2_meta_sanity.R.
# ===========================================================================

message(">> [17] Catalog R2 meta-analysis check")

unlink(file.path(dir_tables, c(
  "17_liability_r2_not_derivable.csv",
  "17_liability_r2_class_counts.csv",
  "17_liability_r2_examples.csv",
  "17_liability_r2_prevalence_sensitivity.csv"
)))

num_pat <- "[-+]?(?:\\d+\\.?\\d*|\\d*\\.\\d+)(?:[eE][-+]?\\d+)?"

parse_metric_field <- function(text, label_regex) {
  pat <- sprintf(
    "%s\\s*=\\s*(%s)(?:\\s*\\[\\s*(%s)\\s*,\\s*(%s)\\s*\\])?",
    label_regex, num_pat, num_pat, num_pat
  )
  m <- str_match(text, regex(pat, ignore_case = TRUE))
  tibble(
    present = !is.na(m[, 2]),
    estimate = suppressWarnings(as.numeric(m[, 2])),
    ci_lower = suppressWarnings(as.numeric(m[, 3])),
    ci_upper = suppressWarnings(as.numeric(m[, 4]))
  )
}

r2_ppm <- read_csv(fp_ppm, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
  transmute(
    performance_id = `PGS Performance Metric (PPM) ID`,
    pgs_id = `Evaluated Score`,
    sampleset_id = `PGS Sample Set (PSS)`,
    pgp_id = `PGS Publication (PGP) ID`,
    reported_trait = `Reported Trait`,
    covariates = `Covariates Included in the Model`,
    other_info = `PGS Performance: Other Relevant Information`,
    pmid = `Publication (PMID)`,
    other_metric = `Other Metric(s)`,
    auc = `Area Under the Receiver-Operating Characteristic Curve (AUROC)`
  )

r2_publications <- read_csv(
  catalog_bulk_file("publications", stamp),
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
) %>%
  transmute(
    pgp_id = `PGS Publication/Study (PGP) ID`,
    first_author = `First Author`,
    title = Title,
    publication_year = substr(`Publication Date`, 1, 4)
  )

# One sample-set id can occupy more than one row. Keep the ancestry
# label only when every row of that id agrees.
r2_samples <- pss_raw %>%
  transmute(
    sampleset_id = `PGS Sample Set (PSS)`,
    n_individuals = suppressWarnings(as.numeric(`Number of Individuals`)),
    n_cases = suppressWarnings(as.numeric(`Number of Cases`)),
    broad_ancestry_raw = `Broad Ancestry Category`,
    cohort = `Cohort(s)`,
    sample_note = `Additional Sample/Cohort Information`
  ) %>%
  group_by(sampleset_id) %>%
  summarise(
    n_individuals = suppressWarnings(max(n_individuals, na.rm = TRUE)),
    any_cases = any(is.finite(n_cases)),
    broad_ancestry = {
      labs <- unique(broad_ancestry_raw[!is.na(broad_ancestry_raw) & nzchar(broad_ancestry_raw)])
      if (length(labs) == 1L) labs else paste(sort(labs), collapse = " | ")
    },
    ancestry_conflict = n_distinct(broad_ancestry_raw[!is.na(broad_ancestry_raw) & nzchar(broad_ancestry_raw)]) > 1,
    cohort = paste(unique(cohort[!is.na(cohort) & nzchar(cohort)]), collapse = " | "),
    sample_note = paste(unique(sample_note[!is.na(sample_note) & nzchar(sample_note)]), collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(
    n_individuals = ifelse(is.infinite(n_individuals), NA_real_, n_individuals),
    ancestry_display = assign_evaluation_ancestry(broad_ancestry),
    is_european = !is.na(broad_ancestry) & broad_ancestry == "European",
    is_non_european = !is.na(broad_ancestry) & nzchar(broad_ancestry) & broad_ancestry != "European"
  )

# The full-model field in this release uses the superscript "R²".
# "PGS R2" and "Incremental R2" use an ASCII "R2".
full_r2 <- parse_metric_field(r2_ppm$other_metric, "(?<![A-Za-z])R²")
pgs_r2 <- parse_metric_field(r2_ppm$other_metric, "PGS R2 \\(no covariates\\)")
incr_r2 <- parse_metric_field(r2_ppm$other_metric, "Incremental R2 \\(full-covars\\)")

r2_metrics <- r2_ppm %>%
  bind_cols(
    full_r2 %>% rename_with(~ paste0("full_", .x)),
    pgs_r2 %>% rename_with(~ paste0("pgs_", .x)),
    incr_r2 %>% rename_with(~ paste0("incr_", .x))
  ) %>%
  mutate(
    pgs_has_ci = pgs_present & is.finite(pgs_ci_lower) & is.finite(pgs_ci_upper),
    full_has_ci = full_present & is.finite(full_ci_lower) & is.finite(full_ci_upper),
    incr_has_ci = incr_present & is.finite(incr_ci_lower) & is.finite(incr_ci_upper),
    has_auc_field = !is.na(auc) & nzchar(auc)
  ) %>%
  left_join(r2_samples, by = "sampleset_id") %>%
  left_join(r2_publications, by = "pgp_id") %>%
  left_join(scores_meta, by = "pgs_id") %>%
  mutate(
    sample_note = coalesce(sample_note, ""),
    is_european = coalesce(is_european, FALSE),
    is_non_european = coalesce(is_non_european, FALSE),
    any_cases = coalesce(any_cases, FALSE)
  )

# -----------------------------------------------------------------------
# Yengo et al. 2022 height. Point R², no interval.
# -----------------------------------------------------------------------

yengo <- r2_metrics %>%
  filter(pgp_id == "PGP000382", reported_trait == "Height") %>%
  arrange(pgs_id, broad_ancestry, sampleset_id)

yengo_multi <- yengo %>% filter(pgs_id == "PGS002802")
yengo_multi_eur <- yengo_multi %>% filter(is_european)
yengo_multi_other <- yengo_multi %>% filter(is_non_european)
yengo_ci_n <- sum(yengo$pgs_has_ci | yengo$full_has_ci | yengo$incr_has_ci)
yengo_outside <- yengo %>%
  filter(pgs_id != "PGS002802", is_non_european, round(full_estimate, 2) < 0.09 | round(full_estimate, 2) > 0.21)

yengo_out <- yengo %>%
  transmute(
    performance_id, pgs_id, pgs_name, n_variants, dev_method,
    gwas_ancestry = anc_gwas_txt,
    reported_trait, sampleset_id, broad_ancestry, ancestry_display,
    n_individuals, cohort, covariates,
    r2 = full_estimate,
    r2_ci_lower = full_ci_lower,
    r2_ci_upper = full_ci_upper,
    score_role = if_else(pgs_id == "PGS002802", "multi-ancestry score", "ancestry-specific score"),
    other_metric, pmid, pgp_id, first_author
  )
write_tbl(yengo_out, "17_r2_yengo_height")

# -----------------------------------------------------------------------
# Named field: PGS R2 (no covariates), with and without an interval.
# -----------------------------------------------------------------------

pgs_ci <- r2_metrics %>% filter(pgs_has_ci)
pgs_no_ci <- r2_metrics %>%
  filter(pgs_present, !pgs_has_ci) %>%
  mutate(
    nagelkerke_note = str_detect(coalesce(other_info, ""), regex("nagelkerke", ignore_case = TRUE)),
    auroc_interval = str_detect(
      coalesce(other_metric, ""),
      regex("PGS AUROC \\(no covariates\\)\\s*=\\s*[-+0-9.eE]+\\s*\\[[^\\]]+\\]", ignore_case = TRUE)
    )
  )

pair_ancestry <- pgs_ci %>%
  group_by(pgs_id, reported_trait) %>%
  summarise(
    n_evaluations = n(),
    n_european = sum(is_european),
    n_non_european = sum(is_non_european),
    ancestries = paste(sort(unique(broad_ancestry)), collapse = "; "),
    manuscript_ancestries = paste(sort(unique(na.omit(ancestry_display))), collapse = "; "),
    pmids = paste(sort(unique(pmid)), collapse = "; "),
    pgp_ids = paste(sort(unique(pgp_id)), collapse = "; "),
    first_authors = paste(sort(unique(first_author)), collapse = "; "),
    publication_years = paste(sort(unique(publication_year)), collapse = "; "),
    incr_present_n = sum(incr_present),
    incr_ci_n = sum(incr_has_ci),
    full_ci_n = sum(full_has_ci),
    any_cases_n = sum(any_cases),
    auc_n = sum(has_auc_field),
    .groups = "drop"
  )

both_ancestry <- pair_ancestry %>%
  filter(n_european > 0L, n_non_european > 0L)

trait_scores <- both_ancestry %>%
  group_by(reported_trait) %>%
  summarise(
    n_scores = n_distinct(pgs_id),
    pmids = paste(sort(unique(unlist(str_split(pmids, "; ")))), collapse = "; "),
    .groups = "drop"
  )
write_tbl(trait_scores, "17_r2_scores_per_trait")

within_paper <- pgs_ci %>%
  semi_join(both_ancestry, by = c("pgs_id", "reported_trait")) %>%
  distinct(pmid, reported_trait, pgs_id) %>%
  count(pmid, reported_trait, name = "n_scores")

eur_notes <- pgs_ci %>%
  semi_join(both_ancestry, by = c("pgs_id", "reported_trait")) %>%
  filter(is_european) %>%
  group_by(pmid, pgs_id, reported_trait) %>%
  summarise(
    n_holdout = sum(str_detect(sample_note, regex("holdout|testing cohort", ignore_case = TRUE))),
    cohorts = paste(sort(unique(cohort)), collapse = "; "),
    .groups = "drop"
  )

non_eur_labels <- pgs_ci %>%
  semi_join(both_ancestry, by = c("pgs_id", "reported_trait")) %>%
  filter(is_non_european) %>%
  count(broad_ancestry, ancestry_display, name = "n_evaluations") %>%
  arrange(desc(n_evaluations))
write_tbl(non_eur_labels, "17_r2_non_european_labels")

pairs_out <- both_ancestry %>%
  left_join(pgs_ci %>% distinct(pgs_id, pgs_name, dev_method), by = "pgs_id") %>%
  arrange(pmids, reported_trait, pgs_id)
write_tbl(pairs_out, "17_r2_pgs_ci_score_trait_pairs")

# An R² interval that is not the named PGS R2 (no covariates) field.
other_ci_pat <- regex(
  "(?i)((?:incremental\\s+variance\\s+explained[^,=]{0,80}|R²|R2))\\s*=\\s*[-+0-9.eE]+\\s*\\[[^\\]]+\\]"
)
other_ci <- r2_metrics %>%
  filter(!pgs_has_ci) %>%
  filter(str_detect(coalesce(other_metric, ""), other_ci_pat))

other_pairs <- other_ci %>%
  group_by(pmid, first_author, publication_year, pgp_id, title, pgs_id, reported_trait) %>%
  summarise(
    n_evaluations = n(),
    n_european = sum(is_european),
    n_non_european = sum(is_non_european),
    ancestries = paste(sort(unique(broad_ancestry)), collapse = "; "),
    covariates = paste(sort(unique(covariates)), collapse = " || "),
    other_info = paste(sort(unique(na.omit(other_info[nzchar(other_info)]))), collapse = " || "),
    .groups = "drop"
  ) %>%
  filter(n_european > 0L, n_non_european > 0L)

paper_from_pairs <- function(pairs, r2_field) {
  pairs %>%
    group_by(pmid, first_author, publication_year, pgp_id, title) %>%
    summarise(
      n_score_trait_pairs = n(),
      n_evaluations = sum(n_evaluations),
      n_traits = n_distinct(reported_trait),
      n_scores = n_distinct(pgs_id),
      max_scores_in_one_trait = max(as.numeric(table(reported_trait))),
      ancestries = paste(sort(unique(unlist(str_split(ancestries, "; ")))), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(r2_field = r2_field, .before = 1)
}

named_pairs_for_papers <- both_ancestry %>%
  transmute(
    pmid = pmids,
    first_author = first_authors,
    publication_year = publication_years,
    pgp_id = pgp_ids,
    title = NA_character_,
    pgs_id, reported_trait,
    n_evaluations, ancestries
  ) %>%
  left_join(
    r2_publications %>% distinct(pgp_id, title_full = title),
    by = "pgp_id"
  ) %>%
  mutate(title = coalesce(title_full, title)) %>%
  select(-title_full)

r2_papers <- bind_rows(
  paper_from_pairs(named_pairs_for_papers, "PGS R2 (no covariates) with interval"),
  paper_from_pairs(other_pairs, "Other R2 with interval")
) %>%
  arrange(r2_field, desc(n_score_trait_pairs))
write_tbl(r2_papers, "17_r2_cross_ancestry_papers")

n_ppm <- nrow(r2_metrics)
n_yengo <- nrow(yengo)
n_pgs_ci <- nrow(pgs_ci)
n_pairs <- nrow(both_ancestry)
n_pairs_2022 <- sum(both_ancestry$pmids == "35324888")
n_pairs_2023 <- sum(both_ancestry$pmids == "37890495")
n_traits <- nrow(trait_scores)
n_traits_one <- sum(trait_scores$n_scores == 1L)
n_traits_two <- sum(trait_scores$n_scores == 2L)
n_within_paper_multi <- sum(within_paper$n_scores > 1L)
n_incr_ci <- sum(pgs_ci$incr_has_ci)
n_full_ci <- sum(pgs_ci$full_has_ci)
n_pgs_ci_auc <- sum(pgs_ci$has_auc_field)
n_holdout_pairs <- sum(eur_notes$pmid == "35324888" & eur_notes$n_holdout == 1L)
yengo_eur_txt <- paste(sprintf("%.3f", sort(yengo_multi_eur$full_estimate)), collapse = " and ")
n_other_papers <- n_distinct(other_pairs$pmid)
n_other_pairs <- nrow(other_pairs)
monti_quant_n <- other_pairs %>%
  filter(pmid == "38908374", str_detect(other_info, fixed("pearson correlation"))) %>%
  nrow()
top_other <- r2_papers %>%
  filter(r2_field == "Other R2 with interval") %>%
  slice_max(n_score_trait_pairs, n = 1, with_ties = FALSE)
other_paper_detail <- if (nrow(top_other) == 0L) {
  "No other publication has an R² interval together with a European and a non-European evaluation."
} else {
  sprintf(
    "Largest other set: %s %s (PMID %s), %d score–trait pairs, up to %d scores in one trait. Of the Monti et al. 2024 pairs, %d are quantitative (the note says the beta is a Pearson correlation) and the covariates field is 0. See 17_r2_cross_ancestry_papers.csv.",
    top_other$first_author, top_other$publication_year, top_other$pmid,
    top_other$n_score_trait_pairs, top_other$max_scores_in_one_trait,
    monti_quant_n
  )
}
multi_eur_lo <- if (nrow(yengo_multi_other) == 0L) NA_real_ else min(yengo_multi_other$full_estimate)
multi_eur_hi <- if (nrow(yengo_multi_other) == 0L) NA_real_ else max(yengo_multi_other$full_estimate)
yengo_min <- if (nrow(yengo) == 0L) NA_real_ else min(yengo$full_estimate)

r2_meta_counts <- tibble(
  question = c(
    "1",
    "1",
    "1",
    "1",
    "2",
    "2",
    "2",
    "2",
    "2",
    "2",
    "2",
    "2"
  ),
  clause = c(
    "Yengo 2022 height evaluations",
    "Yengo 2022 height evaluations with an interval",
    "Named PGS R2 (no covariates) with an interval",
    "Those evaluations that also have an AUROC",
    "Score-trait pairs with a named PGS R2 interval and both ancestries",
    "Those pairs from Tanigawa 2022 (PMID 35324888)",
    "Those pairs from Tanigawa 2023 (PMID 37890495)",
    "Incremental R2 intervals on the named-PGS-R2 evaluations",
    "Full-model R2 intervals on the named-PGS-R2 evaluations",
    "Named PGS R2 reported with no interval",
    "Traits with more than one score inside one paper",
    "Other publications with a different R2 interval and both ancestries"
  ),
  n = c(
    n_yengo,
    yengo_ci_n,
    n_pgs_ci,
    n_pgs_ci_auc,
    n_pairs,
    n_pairs_2022,
    n_pairs_2023,
    n_incr_ci,
    n_full_ci,
    nrow(pgs_no_ci),
    n_within_paper_multi,
    n_other_papers
  ),
  denominator = c(
    n_ppm,
    n_yengo,
    n_ppm,
    n_pgs_ci,
    n_pgs_ci,
    n_pairs,
    n_pairs,
    n_pgs_ci,
    n_pgs_ci,
    n_ppm,
    nrow(within_paper),
    n_other_pairs
  ),
  detail = c(
    sprintf(
      "PMID 36224396. No interval, so no standard error. The multi-ancestry score PGS002802 has European R² %s; its other ancestries run from %.3f to %.3f. Ancestry-specific scores in the same set go down to %.3f. These scores have about a million variants, not the 12,111-variant predictor cited by the reviewer.",
      yengo_eur_txt, multi_eur_lo, multi_eur_hi, yengo_min
    ),
    "Zero intervals. Yengo height cannot be meta-analyzed from the Catalog.",
    "This field is the unadjusted PGS R². It is not an R² attached to an AUROC row.",
    "The named PGS R² interval does not sit on the manuscript's AUROC evaluations.",
    sprintf(
      "Each pair is one score in one trait, with a European evaluation and a non-European evaluation. %d of %d traits have one score; %d traits appear in both papers. See 17_r2_scores_per_trait.csv.",
      n_traits_one, n_traits, n_traits_two
    ),
    sprintf(
      "%d of these pairs label one European sample as a UK Biobank held-out testing cohort. The other European sample is a second UK Biobank partition.",
      n_holdout_pairs
    ),
    "UK Biobank. One score per trait in this paper.",
    "The covariate-adjusted PGS contribution has no interval, so it has no standard error.",
    "The full-model R² includes the covariates. Its interval is not a PGS-only interval.",
    sprintf(
      "All %d are Tanigawa 2022. %d note Nagelkerke's R². %d put the interval on PGS AUROC rather than on the R².",
      nrow(pgs_no_ci),
      sum(pgs_no_ci$nagelkerke_note),
      sum(pgs_no_ci$auroc_interval)
    ),
    sprintf(
      "Within a paper the maximum number of scores in one trait is %d. There is no second score to pool.",
      if (nrow(within_paper) == 0L) NA_integer_ else max(within_paper$n_scores)
    ),
    other_paper_detail
  )
)
write_tbl(r2_meta_counts, "17_r2_meta_counts")

pick_r2_example <- function(df, role) {
  df %>%
    arrange(performance_id) %>%
    slice(1) %>%
    transmute(
      role = role,
      performance_id, pgs_id, pgs_name, reported_trait,
      broad_ancestry, pmid, first_author,
      covariates, other_info, other_metric
    )
}
r2_examples <- bind_rows(
  pick_r2_example(yengo %>% filter(pgs_id == "PGS002802", is_european), "yengo_multi_ancestry_european"),
  pick_r2_example(yengo %>% filter(pgs_id == "PGS002802", is_non_european), "yengo_multi_ancestry_non_european"),
  pick_r2_example(yengo_outside, "yengo_ancestry_specific_outside_0.09_0.21"),
  pick_r2_example(pgs_ci %>% filter(pmid == "35324888"), "tanigawa_2022_pgs_r2_with_interval"),
  pick_r2_example(pgs_ci %>% filter(pmid == "37890495"), "tanigawa_2023_pgs_r2_with_interval"),
  pick_r2_example(pgs_no_ci, "tanigawa_2022_pgs_r2_without_interval"),
  pick_r2_example(other_ci %>% filter(pmid == "38908374"), "monti_2024_other_r2_with_interval")
)
write_tbl(r2_examples, "17_r2_examples")

stopifnot(
  nrow(r2_metrics) == n_ppm,
  n_yengo == 14L,
  yengo_ci_n == 0L,
  all(yengo$full_present),
  nrow(yengo_multi_eur) == 2L,
  n_pgs_ci == 3150L,
  n_pgs_ci_auc == 0L,
  n_pairs == 630L,
  n_pairs_2022 == 570L,
  n_pairs_2023 == 60L,
  n_incr_ci == 0L,
  n_full_ci == n_pgs_ci,
  nrow(pgs_no_ci) == 1189L,
  all(pgs_no_ci$nagelkerke_note),
  all(pgs_no_ci$auroc_interval),
  n_within_paper_multi == 0L,
  n_traits_one == 586L,
  n_traits_two == 22L,
  n_traits_one + n_traits_two == n_traits,
  n_holdout_pairs == n_pairs_2022,
  sum(r2_papers$n_score_trait_pairs[r2_papers$r2_field == "PGS R2 (no covariates) with interval"]) == n_pairs
)

# ===========================================================================
# 00) Sanity check of this run, not the older Catalog snapshot
# ===========================================================================

sanity <- tibble(
  item = c(
    "Catalog stamp",
    "Figure 1C traits",
    "Figure 1C unique cells",
    "Cells excluding zero, logit-scale SE as submitted",
    "Cells excluding zero, delta-method SE",
    "Rebuilt delta-method cells matching published (max |dAUC|)",
    "Rebuilt delta-method cells matching published (max |SE|)",
    "PGS pairs with overlapping performance_id",
    "Figure 1C cells with Multi-ancestry including European",
    "I2 values outside 0-100 (pipeline tables)",
    "Catalog All of Us evaluations",
    "Figure 1C rows whose cohort is All of Us",
    "Pairs whose two sides differ in covariate class"
  ),
  value = c(
    as.numeric(stamp),
    n_forest_traits,
    n_forest_cells,
    scale_summary$n_sig_as_implemented,
    scale_summary$n_sig_delta_method,
    max(cmp_rebuild$abs_diff, na.rm = TRUE),
    max(cmp_rebuild$abs_diff_se, na.rm = TRUE),
    independence_summary$n_pairs_with_overlapping_performance_id,
    independence_summary$n_figure1c_cells_multi_incl_eur,
    sum(i2_check$n_I2_below_0[i2_check$table != "stage1_or_hr (this script)"]) +
      sum(i2_check$n_I2_above_100[i2_check$table != "stage1_or_hr (this script)"]),
    sum(eval_df$catalog_all_of_us, na.rm = TRUE),
    sum(cell_identifiers$catalog_all_of_us, na.rm = TRUE),
    n_class_mismatch
  ),
  check = c(
    stamp,
    as.character(n_forest_traits),
    as.character(n_forest_cells),
    as.character(scale_summary$n_sig_as_implemented),
    as.character(scale_summary$n_sig_delta_method),
    "near zero if script 08 used delta-method SEs",
    "near zero if script 08 used delta-method SEs",
    "should be 0",
    "should be 0 under the single-ancestry rule",
    "should be 0",
    "flag only; not a count lock",
    "flag only; not a count lock",
    "should be 0"
  )
)

write_tbl(sanity, "00_sanity_check_this_run")

message("\n==== Sanity check ====")
print(as.data.frame(sanity), row.names = FALSE)
message("\n==== Scale comparison ====")
print(as.data.frame(scale_summary), row.names = FALSE)
message("\n==== Dependence ====")
print(as.data.frame(dependence_summary), row.names = FALSE)
message("\n==== Covariates ====")
print(as.data.frame(cov_summary), row.names = FALSE)
message("\n==== Training ancestry ====")
print(as.data.frame(train_summary), row.names = FALSE)
message("\n==== OR/HR per SD ====")
print(as.data.frame(or_summary), row.names = FALSE)
message("\n==== Catalog R2 meta-analysis ====")
print(as.data.frame(r2_meta_counts %>% select(question, clause, n, denominator)), row.names = FALSE)
message(glue("\n>> Wrote review-sensitivity pack to {normalizePath(root_out, mustWork = FALSE)}"))
message(">> Done.")
