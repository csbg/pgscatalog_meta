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

pool_stage1 <- function(df) {
  df %>%
    mutate(tmp = pmap(list(auc, estimate_ci_lower, estimate_ci_upper), auc_ci_to_logit)) %>%
    tidyr::unnest_wider(tmp) %>%
    filter(is.finite(eta), is.finite(se), se > 0) %>%
    group_by(pgs_id, trait_label, ancestry_display) %>%
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
  stage1_tbl %>%
    transmute(
      trait_label,
      ancestry_display = as.character(ancestry_display),
      pgs_id,
      auc_pooled, lo_pooled, hi_pooled,
      se_logit = se, k_eval, trained_bucket
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
  scat_df <- scat_df %>% mutate(se_use = .data[[se_col]])

  eur_df <- scat_df %>%
    filter(ancestry_display == "European") %>%
    transmute(
      trait_label, pgs_id, trained_bucket,
      eur_auc = auc_pooled, eur_lo = lo_pooled, eur_hi = hi_pooled,
      eur_se = se_use, eur_se_logit = se_logit, eur_k = k_eval
    )

  tgt_df <- scat_df %>%
    filter(ancestry_display != "European") %>%
    transmute(
      trait_label, pgs_id, trained_bucket,
      target_ancestry = ancestry_display,
      tgt_auc = auc_pooled, tgt_lo = lo_pooled, tgt_hi = hi_pooled,
      tgt_se = se_use, tgt_se_logit = se_logit, tgt_k = k_eval
    )

  tgt_df %>%
    inner_join(eur_df, by = c("trait_label", "pgs_id", "trained_bucket")) %>%
    mutate(
      is_pooled = (eur_k >= 2) & (tgt_k >= 2),
      delta     = tgt_auc - eur_auc,
      delta_se  = sqrt(tgt_se^2 + eur_se^2),
      se_scale  = se_scale
    ) %>%
    filter(is.finite(eur_auc), is.finite(tgt_auc), is.finite(delta), is.finite(delta_se), delta_se > 0)
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
                                label_mode = c("all", "first_k_others_star"),
                                show_k = TRUE) {
  label_mode <- match.arg(label_mode)
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
    scale_x_continuous(breaks = x_breaks_forest, labels = x_labels_forest) +
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

usable_keys <- eval_df %>%
  filter(has_auc, has_ci) %>%
  select(performance_id, pgs_id, trait_label, ancestry_display, sampleset_id, train_bucket, ci_source)

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
    "5. Stage-1 cells (trait x PGS x evaluation ancestry)",
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
    eur_perf = list(eval_keys$performance_id[eval_keys$pgs_id == pgs_id & eval_keys$trait_label == trait_label & eval_keys$ancestry_display == "European"]),
    tgt_perf = list(eval_keys$performance_id[eval_keys$pgs_id == pgs_id & eval_keys$trait_label == trait_label & eval_keys$ancestry_display == target_ancestry]),
    eur_pss  = list(eval_keys$sampleset_id[eval_keys$pgs_id == pgs_id & eval_keys$trait_label == trait_label & eval_keys$ancestry_display == "European"]),
    tgt_pss  = list(eval_keys$sampleset_id[eval_keys$pgs_id == pgs_id & eval_keys$trait_label == trait_label & eval_keys$ancestry_display == target_ancestry])
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
    trait_label, trained_bucket, target_ancestry, pgs_id,
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
  filter(is.finite(auc), is.finite(estimate_ci_lower), is.finite(estimate_ci_upper))

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

cell_pss_sharing <- paired_main %>%
  distinct(trait_label, trained_bucket, target_ancestry) %>%
  mutate(cell_id = row_number()) %>%
  rowwise() %>%
  mutate(
    pgs_in_cell = list(paired_main$pgs_id[
      paired_main$trait_label == trait_label & paired_main$trained_bucket == trained_bucket & paired_main$target_ancestry == target_ancestry
    ]),
    eur_tbl = list(
      usable_keys %>%
        filter(pgs_id %in% pgs_in_cell, trait_label == .env$trait_label, ancestry_display == "European") %>%
        group_by(sampleset_id) %>% summarise(n_pgs = n_distinct(pgs_id), .groups = "drop")
    ),
    tgt_tbl = list(
      usable_keys %>%
        filter(pgs_id %in% pgs_in_cell, trait_label == .env$trait_label, ancestry_display == .env$target_ancestry) %>%
        group_by(sampleset_id) %>% summarise(n_pgs = n_distinct(pgs_id), .groups = "drop")
    )
  ) %>%
  ungroup() %>%
  mutate(
    n_pairs = lengths(pgs_in_cell),
    eur_n_sample_sets = map_int(eur_tbl, nrow),
    eur_max_pgs_sharing_one_sample_set = map_int(eur_tbl, ~ if (nrow(.x)) max(.x$n_pgs) else 0L),
    tgt_n_sample_sets = map_int(tgt_tbl, nrow),
    tgt_max_pgs_sharing_one_sample_set = map_int(tgt_tbl, ~ if (nrow(.x)) max(.x$n_pgs) else 0L),
    all_pgs_share_one_target_sample_set = n_pairs >= 2 & tgt_max_pgs_sharing_one_sample_set == n_pairs,
    all_pgs_share_one_eur_sample_set    = n_pairs >= 2 & eur_max_pgs_sharing_one_sample_set == n_pairs
  ) %>%
  select(-cell_id, -pgs_in_cell, -eur_tbl, -tgt_tbl) %>%
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

pair_tgt_pss <- paired_dm %>%
  select(trait_label, trained_bucket, target_ancestry, pgs_id) %>%
  rowwise() %>%
  mutate(
    tgt_pss_key = paste(sort(unique(usable_keys$sampleset_id[
      usable_keys$pgs_id == pgs_id & usable_keys$trait_label == trait_label & usable_keys$ancestry_display == target_ancestry
    ])), collapse = ";"),
    eur_pss_key = paste(sort(unique(usable_keys$sampleset_id[
      usable_keys$pgs_id == pgs_id & usable_keys$trait_label == trait_label & usable_keys$ancestry_display == "European"
    ])), collapse = ";")
  ) %>%
  ungroup()

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
# 10) Covariate adjustment
#     Classify the 'Covariates Included in the Model' field per PPM; compare
#     EUR vs target sides; re-pool cells restricted to covariate-matched pairs.
# ===========================================================================

message(">> [10] Covariate class")

demo_tech_regex <- paste0(
  "^(age|sex|gender|pc\\d*|pcs?|principal|genetic pc|genotyp|array|batch|cent(er|re)|birth|",
  "deprivation|townsend|assessment|recruitment|ancestry|region|study|cohort|site|platform|chip|country|",
  "year|season|ses\\b|socio|education|income|\\d+\\s*pcs?|pc\\s*1|genetic principal|first \\d+ pc|top \\d+ pc|",
  "ehr|electronic health|duration of|follow-?up|enrol)"
)

classify_covariates <- function(txt) {
  x <- tolower(trimws(coalesce(txt, "")))
  if (x == "" || x %in% c("0", "none", "no", "na", "n/a", "unadjusted", "no covariates", "-")) return("None (PRS only)")
  toks <- str_split(x, "[,;+/]|\\band\\b|\\bwith\\b")[[1]]
  toks <- trimws(toks)
  toks <- toks[toks != ""]
  if (!length(toks)) return("None (PRS only)")
  demo <- grepl(demo_tech_regex, toks, perl = TRUE) |
    grepl("(pc|principal component|genotyp|array|batch|age|sex)", toks, perl = TRUE)
  if (all(demo)) "Age/sex/PCs/technical" else "Clinical or other covariates"
}

ppm_cov <- usable_keys %>%
  left_join(ppm_meta %>% select(performance_id, covariates_txt), by = "performance_id") %>%
  mutate(covariate_class = map_chr(covariates_txt, classify_covariates))

write_tbl(
  ppm_cov %>% count(covariate_class, name = "n_evaluations") %>% mutate(share = n_evaluations / sum(n_evaluations)),
  "10_covariates_class_overall"
)

# Per pair: covariate classes on each side
pair_cov <- paired_dm %>%
  select(trait_label, trained_bucket, target_ancestry, pgs_id, delta, delta_se, is_pooled, eur_k, tgt_k) %>%
  rowwise() %>%
  mutate(
    eur_classes = paste(sort(unique(ppm_cov$covariate_class[ppm_cov$pgs_id == pgs_id & ppm_cov$trait_label == trait_label & ppm_cov$ancestry_display == "European"])), collapse = " | "),
    tgt_classes = paste(sort(unique(ppm_cov$covariate_class[ppm_cov$pgs_id == pgs_id & ppm_cov$trait_label == trait_label & ppm_cov$ancestry_display == target_ancestry])), collapse = " | "),
    eur_covariates_txt = paste(unique(coalesce(ppm_cov$covariates_txt[ppm_cov$pgs_id == pgs_id & ppm_cov$trait_label == trait_label & ppm_cov$ancestry_display == "European"], "")), collapse = " || "),
    tgt_covariates_txt = paste(unique(coalesce(ppm_cov$covariates_txt[ppm_cov$pgs_id == pgs_id & ppm_cov$trait_label == trait_label & ppm_cov$ancestry_display == target_ancestry], "")), collapse = " || ")
  ) %>%
  ungroup() %>%
  mutate(
    covariate_matched = eur_classes == tgt_classes & !grepl("\\|", eur_classes),
    either_side_clinical = grepl("Clinical", eur_classes) | grepl("Clinical", tgt_classes),
    either_side_mixed_classes = grepl("\\|", eur_classes) | grepl("\\|", tgt_classes)
  )

write_tbl(pair_cov, "10_covariates_pairs")

cells_cov_matched <- paired_dm %>%
  inner_join(pair_cov %>% filter(covariate_matched) %>% select(trait_label, trained_bucket, target_ancestry, pgs_id),
             by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")) %>%
  pool_stage2(model = "fixed")

cov_cells <- delta_dm %>%
  transmute(trait_label, trained_bucket, target_ancestry, n_pairs, delta_all = delta_hat, lo_all = lo, hi_all = hi, sig_all = ci_excludes_zero) %>%
  left_join(
    pair_cov %>% group_by(trait_label, trained_bucket, target_ancestry) %>%
      summarise(n_pairs_covariate_matched = sum(covariate_matched), n_pairs_either_clinical = sum(either_side_clinical),
                n_pairs_mixed_classes = sum(either_side_mixed_classes), .groups = "drop"),
    by = c("trait_label", "trained_bucket", "target_ancestry")
  ) %>%
  left_join(
    cells_cov_matched %>% transmute(trait_label, trained_bucket, target_ancestry,
                                    delta_matched = delta_hat, lo_matched = lo, hi_matched = hi, sig_matched = ci_excludes_zero),
    by = c("trait_label", "trained_bucket", "target_ancestry")
  ) %>%
  mutate(inference_flip_matched = !is.na(sig_matched) & sig_matched != sig_all)

write_tbl(cov_cells, "10_covariates_cells_matched_sensitivity")

cov_summary <- tibble(
  n_pairs = nrow(pair_cov),
  n_pairs_covariate_matched = sum(pair_cov$covariate_matched),
  n_pairs_either_side_clinical = sum(pair_cov$either_side_clinical),
  n_pairs_mixed_classes_within_side = sum(pair_cov$either_side_mixed_classes),
  n_cells = nrow(cov_cells),
  n_cells_estimable_matched_only = sum(!is.na(cov_cells$delta_matched)),
  n_cells_inference_flip_matched_only = sum(cov_cells$inference_flip_matched)
)
write_tbl(cov_summary, "10_covariates_summary")

cov_long <- bind_rows(
  delta_dm %>% transmute(trait_label, trained_bucket, target_ancestry, n_pairs, delta_hat, lo, hi, source = "All pairs (corrected)"),
  cells_cov_matched %>% transmute(trait_label, trained_bucket, target_ancestry, n_pairs, delta_hat, lo, hi, source = "Covariate-matched pairs only")
) %>% mutate(source = factor(source, levels = c("All pairs (corrected)", "Covariate-matched pairs only")))

p_cov <- forest_compare_plot(
  cov_long,
  pal = c("#0072B2", "#CC79A7"), shapes = c(16, 17),
  title = "ΔAUC restricted to pairs with the same covariate class on both sides",
  subtitle = glue("{cov_summary$n_pairs_covariate_matched} of {cov_summary$n_pairs} pairs matched; {cov_summary$n_cells_estimable_matched_only} of {cov_summary$n_cells} cells estimable; {cov_summary$n_cells_inference_flip_matched_only} inference changes."),
  caption = "Covariate class from the Catalog 'Covariates Included in the Model' field: None (PRS only) / Age-sex-PCs-technical / Clinical or other. Matched = identical single class on both sides.",
  lost_df = cov_cells %>% filter(is.na(delta_matched)) %>% select(trait_label, trained_bucket, target_ancestry),
  lost_label = "no matched pair", trait_levels = trait_levels_wrapped
)
save_plot("10_covariates_forest", p_cov, width = 12, height = max(6.5, 0.30 * n_row_plot + 2.6), limitsize = FALSE)

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
  paired_dm %>% transmute(trait_label, trained_bucket, target_ancestry, pgs_id, side = "European", ancestry = "European"),
  paired_dm %>% transmute(trait_label, trained_bucket, target_ancestry, pgs_id, side = "Target", ancestry = target_ancestry)
) %>%
  # many-to-many is intended: each pair side expands to all its evaluations
  inner_join(usable_keys %>% select(pgs_id, trait_label, ancestry_display, sampleset_id, performance_id, ci_source),
             by = c("pgs_id", "trait_label", "ancestry" = "ancestry_display"),
             relationship = "many-to-many") %>%
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
# 17) Liability-scale R² is not consistently derivable
#     The Catalog stores an evaluation case fraction and a free-text R².
#     Converting AUC or an OR per SD to liability-scale R² also needs the
#     population prevalence (Wray et al. 2010; Lee et al. 2012).
# ===========================================================================

message(">> [17] Liability-scale R2 check")

# Population-sample mapping from AUC to liability-scale R²
# (Wray, Yang, Goddard & Visscher, PLoS Genet 2010). Ascertained
# case-control samples need the further correction in Lee, Goddard,
# Wray & Visscher, Genet Epidemiol 2012, which requires K and P.
auc_to_liability_r2 <- function(auc, K) {
  r2_to_auc <- function(R2, k) {
    thr <- qnorm(1 - k)
    z <- dnorm(thr)
    mean_c <- z / k
    mean_u <- -z / (1 - k)
    var_c <- 1 + thr * (z / k) - mean_c^2
    var_u <- 1 - thr * (z / (1 - k)) - mean_u^2
    mg_c <- sqrt(R2) * mean_c
    mg_u <- sqrt(R2) * mean_u
    vg_c <- 1 - R2 * (1 - var_c)
    vg_u <- 1 - R2 * (1 - var_u)
    pnorm((mg_c - mg_u) / sqrt(vg_c + vg_u))
  }
  mapply(function(a, k) {
    if (!is.finite(a) || !is.finite(k) || k <= 0 || k >= 1 || a <= 0.5) return(NA_real_)
    lo <- 1e-8
    hi <- 0.95
    if (r2_to_auc(hi, k) < a) return(NA_real_)
    for (.step in seq_len(60)) {
      mid <- 0.5 * (lo + hi)
      if (r2_to_auc(mid, k) < a) lo <- mid else hi <- mid
    }
    0.5 * (lo + hi)
  }, auc, K)
}

r2_classified <- eval_df %>%
  mutate(
    other_metric_txt = coalesce(other_metric_txt, ""),
    other_info_txt = coalesce(other_info_txt, ""),
    r2_blob = tolower(paste(other_metric_txt, other_info_txt, sep = " || ")),
    mentions_liability = str_detect(r2_blob, "liability"),
    states_assumed_prevalence = mentions_liability & str_detect(r2_blob, "assum"),
    r2_class = case_when(
      mentions_liability & str_detect(r2_blob, "nested|incremental") ~ "liability_nested_or_incremental",
      mentions_liability & str_detect(r2_blob, "assum|prev") ~ "liability_with_assumed_prevalence",
      mentions_liability ~ "liability_scale_other",
      str_detect(r2_blob, "nagelkerke") ~ "nagelkerke_pseudo_r2",
      str_detect(r2_blob, "pseudo") ~ "other_pseudo_r2",
      str_detect(r2_blob, "incremental|delta r|nested") ~ "incremental_or_nested_not_called_liability",
      str_detect(r2_blob, "\\br2\\b|r\u00b2|variance explained") ~ "unspecified_r2",
      !str_detect(other_metric_txt, "\\S") &
        !str_detect(other_info_txt, regex("r2|r\u00b2", ignore_case = TRUE)) ~ "no_r2_reported",
      TRUE ~ "other_metric_not_r2"
    )
  )

fig_ppm_ids <- unique(cell_identifiers$performance_id)
auroc_cls <- r2_classified %>% filter(has_auc)
fig_cls <- r2_classified %>%
  filter(performance_id %in% fig_ppm_ids) %>%
  distinct(performance_id, .keep_all = TRUE)

count_class <- function(df, scope_name) {
  df %>%
    count(r2_class, name = "n") %>%
    mutate(scope = scope_name, denominator = sum(n), .before = 1)
}
class_counts <- bind_rows(
  count_class(auroc_cls, "auroc_metrics"),
  count_class(fig_cls, "figure_1c_metrics")
) %>%
  arrange(scope, desc(n))

n_auroc <- nrow(auroc_cls)
n_auroc_liability <- sum(auroc_cls$mentions_liability)
n_auroc_liability_nested <- sum(auroc_cls$r2_class == "liability_nested_or_incremental")
n_auroc_liability_other <- sum(auroc_cls$r2_class == "liability_scale_other")
n_fig_ppm <- nrow(fig_cls)
n_fig_liability <- sum(fig_cls$mentions_liability)
n_ppm_liability <- sum(r2_classified$mentions_liability)
n_assumed_prevalence <- sum(r2_classified$states_assumed_prevalence)
# The sample-set file repeats a PSS identifier, and some repeats disagree
# on cases and controls. Count the file, not the de-duplicated `pss` table.
pss_count_rows <- tibble(
  sampleset_id = as.character(pss_raw[[pss_field_names[1]]]),
  n_cases = suppressWarnings(as.numeric(pss_raw[["Number of Cases"]])),
  n_controls = suppressWarnings(as.numeric(pss_raw[["Number of Controls"]]))
)
pss_case_variants <- pss_count_rows %>%
  filter(is.finite(n_cases), is.finite(n_controls)) %>%
  distinct(sampleset_id, n_cases, n_controls) %>%
  mutate(sample_case_fraction = n_cases / (n_cases + n_controls))
n_pss_rows <- nrow(pss_count_rows)
n_pss_ids <- n_distinct(pss_count_rows$sampleset_id)
n_pss_rows_with_counts <- sum(is.finite(pss_count_rows$n_cases) & is.finite(pss_count_rows$n_controls))
n_pss_conflicting <- pss_case_variants %>%
  count(sampleset_id, name = "n_case_control_pairs") %>%
  summarise(n = sum(n_case_control_pairs > 1)) %>%
  pull(n)
n_prevalence_fields <- sum(grepl("prevalence|liability", pss_field_names, ignore.case = TRUE))
fig_pss_ids <- unique(cell_identifiers$sampleset_id)
n_fig_pss <- length(fig_pss_ids)
fig_pss_conflict_ids <- pss_case_variants %>%
  filter(sampleset_id %in% fig_pss_ids) %>%
  count(sampleset_id, name = "n_case_control_pairs") %>%
  filter(n_case_control_pairs > 1) %>%
  arrange(sampleset_id) %>%
  pull(sampleset_id)
n_fig_pss_conflicting <- length(fig_pss_conflict_ids)

class_n <- function(scope_name, class_name) {
  hit <- class_counts$n[class_counts$scope == scope_name & class_counts$r2_class == class_name]
  if (length(hit) == 0) 0L else hit
}

# Same AUC, three population prevalences: the mapping is not unique.
k_grid <- c(0.01, 0.10, 0.50)
r2_at_auc60 <- auc_to_liability_r2(rep(0.60, length(k_grid)), k_grid)
auc60_detail <- paste(
  sprintf("K = %.2f -> R2 = %.3f", k_grid, r2_at_auc60),
  collapse = "; "
)

liability_proof <- tibble(
  clause = c(
    "Performance metrics with a finite AUC",
    "AUROC metrics that mention the liability scale",
    "AUROC metrics with nested or incremental liability R2",
    "AUROC metrics that name liability R2 without a further definition",
    "AUROC metrics reporting Nagelkerke pseudo-R2",
    "AUROC metrics reporting another pseudo-R2",
    "AUROC metrics reporting an unspecified R2",
    "AUROC metrics reporting no R2",
    "AUROC metrics whose other metric is not an R2",
    "Figure 1C performance metrics",
    "Figure 1C metrics that mention the liability scale",
    "Figure 1C metrics reporting an unspecified R2",
    "Figure 1C metrics reporting Nagelkerke pseudo-R2",
    "Figure 1C metrics reporting no R2",
    "Performance metrics in the release that mention liability anywhere",
    "Performance metrics that state an author-assumed prevalence",
    "Evaluation sample-set rows in the release file",
    "Distinct evaluation sample-set identifiers",
    "Sample-set rows with case and control counts",
    "Distinct sample-set identifiers with more than one case/control pair",
    "Figure 1C sample sets with more than one case/control pair",
    "Sample-set columns that record population prevalence",
    "Figure 1C cells with a paired OR or HR",
    "Liability-scale R2 implied by AUC 0.60 at three prevalences"
  ),
  n = c(
    n_auroc,
    n_auroc_liability,
    n_auroc_liability_nested,
    n_auroc_liability_other,
    class_n("auroc_metrics", "nagelkerke_pseudo_r2"),
    class_n("auroc_metrics", "other_pseudo_r2"),
    class_n("auroc_metrics", "unspecified_r2"),
    class_n("auroc_metrics", "no_r2_reported"),
    class_n("auroc_metrics", "other_metric_not_r2"),
    n_fig_ppm,
    n_fig_liability,
    class_n("figure_1c_metrics", "unspecified_r2"),
    class_n("figure_1c_metrics", "nagelkerke_pseudo_r2"),
    class_n("figure_1c_metrics", "no_r2_reported"),
    n_ppm_liability,
    n_assumed_prevalence,
    n_pss_rows,
    n_pss_ids,
    n_pss_rows_with_counts,
    n_pss_conflicting,
    n_fig_pss_conflicting,
    n_prevalence_fields,
    or_summary$n_figure1c_cells_with_or_estimate,
    NA_real_
  ),
  denominator = c(
    n_ppm,
    n_auroc,
    n_auroc,
    n_auroc,
    n_auroc,
    n_auroc,
    n_auroc,
    n_auroc,
    n_auroc,
    n_fig_ppm,
    n_fig_ppm,
    n_fig_ppm,
    n_fig_ppm,
    n_fig_ppm,
    n_ppm,
    n_ppm,
    n_pss_rows,
    n_pss_rows,
    n_pss_rows,
    n_pss_ids,
    n_fig_pss,
    length(pss_field_names),
    or_summary$n_figure1c_cells,
    NA_real_
  ),
  detail = c(
    "Structured AUROC. This is the denominator of the liability count.",
    "Free-text match for 'liability' in Other Metric(s) or the performance note. Subclasses are the next two rows.",
    "The string also says nested or incremental, so the value is not a marginal liability-scale R2.",
    "Named liability R2 with no nested model and no stated prevalence.",
    "Nagelkerke pseudo-R2 is not a liability-scale R2.",
    "Covariate-adjusted or other pseudo-R2, not labelled liability-scale.",
    "The string says R2 or variance explained and does not name the scale.",
    "No R2 in the other-metric field or the performance note.",
    "Another metric is present (for example AUPRC) and it is not an R2.",
    "Distinct PPM identifiers in 13_cell_identifiers_long.csv.",
    "None of the Figure 1C metrics names the liability scale.",
    "Bare R2 string; scale not stated.",
    "Info field identifies Nagelkerke's method.",
    "Empty R2 fields.",
    "Includes records whose AUC sits only in free text, so they sit outside the finite-AUROC count.",
    "The prevalence is an author assumption written into the metric string, not a Catalog field. See 17_liability_r2_examples.csv.",
    paste0(basename(fp_pss), ". One sample set can occupy more than one row."),
    "Unique PGS Sample Set (PSS) identifiers in that file.",
    "These rows support an evaluation case fraction P = cases / (cases + controls). P is not the population prevalence K.",
    paste0(
      "The same PSS identifier is stored with different numeric case and control counts, so P is not a single Catalog value. In Figure 1C: ",
      paste(fig_pss_conflict_ids, collapse = ", "),
      "."
    ),
    paste0("Figure 1C uses 13_cell_identifiers_long.csv. ", n_fig_pss_conflicting, " of its sample sets have more than one case/control pair in the release file."),
    "No sample-set column matches prevalence or liability.",
    "An odds ratio or hazard ratio still needs K, and per-SD scaling is not a structured field. See 12_or_per_sd_summary.csv.",
    paste0(
      "Wray et al. 2010 population-sample mapping. ",
      auc60_detail,
      ". Lee et al. 2012 changes these values again when the sample is ascertained. See 17_liability_r2_prevalence_sensitivity.csv."
    )
  )
)
write_tbl(liability_proof, "17_liability_r2_not_derivable")
write_tbl(class_counts, "17_liability_r2_class_counts")

pick_example <- function(df, role) {
  df %>%
    arrange(performance_id) %>%
    slice(1) %>%
    transmute(
      role = role,
      performance_id, pgs_id, reported_trait, has_auc, auc, r2_class,
      other_metric_txt, other_info_txt
    )
}
liability_examples <- bind_rows(
  pick_example(auroc_cls %>% filter(r2_class == "liability_scale_other"), "auroc_liability_r2_named"),
  pick_example(auroc_cls %>% filter(r2_class == "liability_nested_or_incremental"), "auroc_liability_r2_nested"),
  pick_example(r2_classified %>% filter(!has_auc, mentions_liability, !states_assumed_prevalence), "liability_scale_only_in_free_text"),
  pick_example(r2_classified %>% filter(states_assumed_prevalence), "author_assumed_prevalence"),
  pick_example(fig_cls %>% filter(r2_class == "unspecified_r2"), "figure_1c_unspecified_r2"),
  pick_example(fig_cls %>% filter(r2_class == "nagelkerke_pseudo_r2"), "figure_1c_nagelkerke"),
  pick_example(fig_cls %>% filter(r2_class == "no_r2_reported"), "figure_1c_no_r2")
)
write_tbl(liability_examples, "17_liability_r2_examples")

# Figure 1C ascertainment. AUC is taken from the performance metric.
# Case fractions are every distinct case/control pair the release file
# stores for that sample set, including disagreements.
ascertainment_examples <- cell_identifiers %>%
  filter(
    (pgs_id == "PGS000004" & trait_label == "Breast cancer" & target_ancestry == "African") |
      (pgs_id == "PGS002250" & grepl("ovarian", trait_label, ignore.case = TRUE))
  ) %>%
  distinct(pgs_id, trait_label, side, performance_id, sampleset_id, auc) %>%
  filter(is.finite(auc)) %>%
  inner_join(pss_case_variants, by = "sampleset_id", relationship = "many-to-many") %>%
  group_by(sampleset_id) %>%
  mutate(n_case_control_pairs = n_distinct(paste(n_cases, n_controls))) %>%
  ungroup()

shared_k <- c(0.01, 0.05, 0.10, 0.50)
liability_sensitivity <- bind_rows(
  tibble(
    setting = "fixed_auc",
    pgs_id = NA_character_,
    trait_label = NA_character_,
    side = NA_character_,
    performance_id = NA_character_,
    sampleset_id = NA_character_,
    auc = 0.60,
    n_cases = NA_real_,
    n_controls = NA_real_,
    prevalence_used = k_grid,
    prevalence_role = "external_population_prevalence",
    n_case_control_pairs_in_catalog = NA_integer_
  ),
  ascertainment_examples %>%
    transmute(
      setting = "catalog_case_fraction_used_as_k",
      pgs_id, trait_label, side, performance_id, sampleset_id, auc,
      n_cases, n_controls,
      prevalence_used = sample_case_fraction,
      prevalence_role = ifelse(
        n_case_control_pairs > 1,
        "catalog_lists_more_than_one_case_fraction",
        "single_catalog_case_fraction_used_as_k"
      ),
      n_case_control_pairs_in_catalog = n_case_control_pairs
    ),
  ascertainment_examples %>%
    distinct(pgs_id, trait_label, side, performance_id, sampleset_id, auc, n_case_control_pairs) %>%
    crossing(prevalence_used = shared_k) %>%
    transmute(
      setting = "shared_external_k",
      pgs_id, trait_label, side, performance_id, sampleset_id, auc,
      n_cases = NA_real_,
      n_controls = NA_real_,
      prevalence_used,
      prevalence_role = "one_external_k_applied_to_the_auc",
      n_case_control_pairs_in_catalog = n_case_control_pairs
    )
) %>%
  mutate(
    r2_liability = auc_to_liability_r2(auc, prevalence_used),
    formula = "Wray et al. 2010 population-sample liability-threshold mapping. Not reported as a result; it shows that K must be supplied."
  )
write_tbl(liability_sensitivity, "17_liability_r2_prevalence_sensitivity")

stopifnot(
  sum(class_counts$n[class_counts$scope == "auroc_metrics"]) == n_auroc,
  sum(class_counts$n[class_counts$scope == "figure_1c_metrics"]) == n_fig_ppm,
  n_prevalence_fields == 0
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
    "Figure 1C rows whose cohort is All of Us"
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
    sum(cell_identifiers$catalog_all_of_us, na.rm = TRUE)
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
    "flag only; not a count lock"
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
message("\n==== Liability-scale R2 ====")
print(as.data.frame(liability_proof %>% select(clause, n, denominator)), row.names = FALSE)
message(glue("\n>> Wrote review-sensitivity pack to {normalizePath(root_out, mustWork = FALSE)}"))
message(">> Done.")
