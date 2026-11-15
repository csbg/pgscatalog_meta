#!/usr/bin/env Rscript

# Forest slides for type 2 diabetes, European-only training.
# Top: one Δ AUC per score (pipeline 1, stage 2).
# Bottom: the AUC models from analysis/10_delta_auc_models.Rmd:
# pipeline 1, pipeline 3, and pipeline 5.
#
# Run from the repository root:
#   Rscript tests/additional_forests.R
#
# Words in titles and captions: evaluation, score, sample set, comparison,
# European, non-European, average of all scores. In a model formula,
# (1|cohort) means the sample set. Do not write cell, code, pool,
# cohorts left, key, shape, bucket, or k.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(janitor)
  library(lme4)
  library(metafor)
  library(patchwork)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

if (!file.exists("R/load.R")) {
  stop("Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

stamp <- pipeline_stamp
base_dir <- file.path("results", paste0("review_sensitivity_", stamp))
table_dir <- file.path(base_dir, "tables")
root_dir <- file.path(base_dir, "additional_forests")
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

read_need <- function(path) {
  if (!file.exists(path)) stop("Missing ", path, ".", call. = FALSE)
  read_csv(path, show_col_types = FALSE)
}

# Same construction as the load-data chunk in analysis/10_delta_auc_models.Rmd.
fp_eval <- file.path("results", paste0("pgs_auc_ci_audit_", stamp), "eval_df_final_auc_ci.csv")
fp_s1 <- file.path(
  catalog_results_dir(stamp), "stage1",
  paste0("stage1_pooled_cells_I2filtered_", stamp, ".csv")
)
fp_s2 <- file.path(
  catalog_results_dir(stamp), "stage2",
  paste0("deltaAUC_forest_faceted_data_", stamp, ".csv")
)
eval_raw <- read_need(fp_eval)
stage1_saved <- read_need(fp_s1)
stage2_saved <- read_need(fp_s2)

pss <- read_csv(
  catalog_bulk_file("evaluation_sample_sets", stamp),
  show_col_types = FALSE, col_types = cols(.default = col_character())
) %>%
  clean_names() %>%
  transmute(
    sampleset_id = as.character(pgs_sample_set_pss),
    n_individuals = suppressWarnings(as.numeric(number_of_individuals)),
    n_cases = suppressWarnings(as.numeric(number_of_cases)),
    n_controls = suppressWarnings(as.numeric(number_of_controls)),
    cohort_txt = as.character(cohort_s)
  ) %>%
  distinct(sampleset_id, .keep_all = TRUE)

ppm_cov <- read_csv(
  catalog_bulk_file("performance_metrics", stamp),
  show_col_types = FALSE, col_types = cols(.default = col_character())
) %>%
  clean_names() %>%
  transmute(
    performance_id = as.character(pgs_performance_metric_ppm_id),
    covariates = covariates_included_in_the_model
  ) %>%
  distinct(performance_id, .keep_all = TRUE)

eval_proc <- eval_raw %>%
  mutate(
    pgs_id = as.character(pgs_id),
    sampleset_id = as.character(sampleset_id),
    performance_id = as.character(performance_id),
    train_bucket = as.character(train_bucket),
    trait_label = harmonize_trait_label(
      NA_character_,
      coalesce(na_if(trimws(as.character(reported_trait)), ""), NA_character_)
    ),
    ancestry_display = as.character(ancestry_eval)
  ) %>%
  mutate(across(c(auc, estimate_ci_lower, estimate_ci_upper), as.numeric)) %>%
  filter(
    .data$ancestry_display %in% evaluation_display_levels,
    is.finite(.data$auc), is.finite(.data$estimate_ci_lower), is.finite(.data$estimate_ci_upper)
  ) %>%
  left_join(
    pss %>% select(sampleset_id, n = n_individuals, cases = n_cases, controls = n_controls),
    by = "sampleset_id"
  ) %>%
  left_join(ppm_cov, by = "performance_id") %>%
  mutate(covariate_class = classify_covariates(covariates))

one_eval <- keep_one_evaluation(eval_proc, class_col = "covariate_class", return_dropped = FALSE)
t2d_eval <- one_eval %>%
  filter(.data$trait_label == "Type 2 diabetes", .data$train_bucket == "Train: European-only")

t2d_stage1 <- stage1_saved %>%
  filter(.data$trait_label == "Type 2 diabetes") %>%
  mutate(
    pgs_id = as.character(pgs_id),
    ancestry_display = as.character(ancestry_display),
    covariate_class = as.character(covariate_class),
    auc_pooled = plogis(eta)
  ) %>%
  inner_join(t2d_eval %>% distinct(pgs_id, trained_bucket = train_bucket), by = "pgs_id")
stage2_pairs <- pair_stage1_by_class(t2d_stage1, se_col = "se") %>%
  pick_one_class_per_score()

make_master <- function(target) {
  keep <- stage2_pairs %>%
    filter(.data$target_ancestry == target) %>%
    select(pgs_id, covariate_class)
  t2d_eval %>%
    filter(.data$ancestry_display %in% c("European", target)) %>%
    inner_join(keep, by = c("pgs_id", "covariate_class")) %>%
    left_join(pss %>% select(sampleset_id, cohort_txt, n_individuals), by = "sampleset_id") %>%
    mutate(tmp = pmap(list(auc, estimate_ci_lower, estimate_ci_upper), auc_ci_to_logit)) %>%
    unnest_wider(tmp) %>%
    transmute(
      pgs_id, ancestry = factor(ancestry_display, levels = c("European", target)),
      sampleset_id, study = normalize_study(cohort_txt), n_individuals, covariate_class,
      auc, ci_lower = estimate_ci_lower, ci_upper = estimate_ci_upper,
      logit_auc = eta, se_logit = se, se_auc = se * auc * (1 - auc)
    ) %>%
    mutate(PGS = factor(pgs_id), cohort = factor(sampleset_id)) %>%
    arrange(ancestry, pgs_id, sampleset_id)
}

# Names match the model column in the report.
auc_model_names <- c(
  "1. IVW stage 1 → stage 2",
  "3.1 lm(AUC ~ ancestry), weighted",
  "3.2 lm(AUC ~ ancestry + PGS)",
  "3.3 lm(AUC ~ ancestry + PGS), weighted",
  "3.4 lmer(AUC ~ ancestry + PGS + (1|cohort)), weighted",
  "3.5 lmer(AUC ~ ancestry + (1|PGS) + (1|cohort)), weighted",
  "5.1 rma.mv(AUC ~ ancestry)",
  "5.2 rma.mv(AUC ~ ancestry + PGS)",
  "5.3 rma.mv(AUC ~ ancestry + PGS, random = 1|cohort)",
  "5.4 rma.mv(AUC ~ ancestry, random = 1|PGS + 1|cohort)"
)

anc_from_fit <- function(fit) {
  if (inherits(fit, "rma")) {
    b <- coef(fit)
    V <- vcov(fit)
  } else if (inherits(fit, "merMod")) {
    b <- fixef(fit)
    V <- as.matrix(vcov(fit))
  } else {
    b <- coef(fit)
    V <- vcov(fit)
  }
  if (length(b) < 2 || !grepl("^ancestry", names(b)[2])) {
    stop(
      "Ancestry coefficient was not the second term: ",
      paste(names(b), collapse = ", "), ".",
      call. = FALSE
    )
  }
  est <- unname(b[[2]])
  se <- sqrt(unname(V[2, 2]))
  c(estimate = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se)
}

mer_note <- function(fit) {
  if (!inherits(fit, "merMod")) return("")
  out <- character()
  if (isSingular(fit)) out <- c(out, "singular fit")
  if (length(fit@optinfo$conv$lme4$messages)) out <- c(out, "convergence warning")
  paste(out, collapse = "; ")
}

row_from_est <- function(model, kind, est, se, note = "") {
  if (!is.finite(est) || !is.finite(se) || se < 0) {
    stop("Non-finite estimate for ", model, ".", call. = FALSE)
  }
  tibble(
    model = model, kind = kind, note = note,
    estimate = est, se = se,
    lo = est - 1.96 * se, hi = est + 1.96 * se,
    hollow = kind != "avg"
  )
}

# Per-score Δ AUC, the stage-2 rows of pipeline 1, plus the ten AUC models.
prepare_comparison <- function(target) {
  m <- make_master(target)
  s1 <- m %>%
    group_by(.data$pgs_id, .data$ancestry) %>%
    summarise(pooled = list(ivw_pool_logit(logit_auc, se_logit)), .groups = "drop") %>%
    unnest_wider(pooled) %>%
    mutate(auc_s1 = plogis(.data$eta), se_auc_s1 = .data$se * .data$auc_s1 * (1 - .data$auc_s1))
  wide <- s1 %>%
    mutate(ancestry = as.character(.data$ancestry)) %>%
    select(pgs_id, ancestry, auc_s1, se_auc_s1) %>%
    pivot_wider(names_from = ancestry, values_from = c(auc_s1, se_auc_s1))
  tgt_auc <- paste0("auc_s1_", target)
  tgt_se <- paste0("se_auc_s1_", target)
  sets <- m %>%
    filter(as.character(.data$ancestry) == target) %>%
    group_by(.data$pgs_id) %>%
    summarise(
      tgt_sampleset_ids = paste(sort(unique(.data$sampleset_id)), collapse = ";"),
      .groups = "drop"
    )
  scores <- wide %>%
    transmute(
      pgs_id,
      delta = .data[[tgt_auc]] - .data$auc_s1_European,
      delta_se = sqrt(.data[[tgt_se]]^2 + .data$se_auc_s1_European^2)
    ) %>%
    inner_join(sets, by = "pgs_id") %>%
    arrange(.data$tgt_sampleset_ids, .data$pgs_id)
  if (any(!is.finite(scores$delta)) || any(!is.finite(scores$delta_se))) {
    stop("Non-finite per-score Δ AUC for ", target, ".", call. = FALSE)
  }

  saved <- stage2_saved %>%
    filter(
      .data$trait_label == "Type 2 diabetes",
      .data$trained_bucket == "Train: European-only",
      .data$target_ancestry == target
    )
  if (nrow(saved) != 1L) {
    stop("Expected one saved stage-2 row for ", target, ".", call. = FALSE)
  }
  w <- 1 / scores$delta_se^2
  avg <- sum(w * scores$delta) / sum(w)
  avg_se <- sqrt(1 / sum(w))
  if (abs(avg - saved$delta_hat) > 1e-8 || abs(avg_se - saved$se_hat) > 1e-8) {
    stop("Pipeline 1 does not match the saved stage-2 table for ", target, ".", call. = FALSE)
  }

  fit_31 <- lm(auc ~ ancestry, data = m, weights = 1 / se_auc^2)
  fit_32 <- lm(auc ~ ancestry + PGS, data = m)
  fit_33 <- lm(auc ~ ancestry + PGS, data = m, weights = 1 / se_auc^2)
  fit_34 <- lmer(auc ~ ancestry + PGS + (1 | cohort), data = m, weights = 1 / se_auc^2)
  fit_35 <- lmer(auc ~ ancestry + (1 | PGS) + (1 | cohort), data = m, weights = 1 / se_auc^2)
  fit_51 <- rma.mv(auc, se_auc^2, mods = ~ ancestry, data = m)
  fit_52 <- rma.mv(auc, se_auc^2, mods = ~ ancestry + PGS, data = m)
  fit_53 <- rma.mv(auc, se_auc^2, mods = ~ ancestry + PGS, random = ~ 1 | cohort, data = m)
  fit_54 <- rma.mv(auc, se_auc^2, mods = ~ ancestry, random = list(~ 1 | PGS, ~ 1 | cohort), data = m)

  fitted <- list(fit_31, fit_32, fit_33, fit_34, fit_35, fit_51, fit_52, fit_53, fit_54)
  kinds <- c("lm", "lm", "lm", "lmer", "lmer", "rma", "rma", "rma", "rma")
  notes <- c("", "", "", mer_note(fit_34), mer_note(fit_35), "", "", "", "")
  model_rows <- purrr::pmap_dfr(list(auc_model_names[-1], kinds, fitted, notes), function(model, kind, fit, note) {
    a <- anc_from_fit(fit)
    row_from_est(model, kind, a[["estimate"]], a[["se"]], note)
  })
  models <- bind_rows(
    row_from_est(auc_model_names[1], "avg", avg, avg_se),
    model_rows
  )
  if (!identical(models$model, auc_model_names)) {
    stop("AUC model rows are out of order for ", target, ".", call. = FALSE)
  }

  meta <- m %>%
    filter(as.character(.data$ancestry) == target) %>%
    distinct(.data$sampleset_id, .data$study, .data$n_individuals)
  if (anyDuplicated(meta$sampleset_id) || any(!is.finite(meta$n_individuals))) {
    stop("Sample-set labels are incomplete for ", target, ".", call. = FALSE)
  }

  list(scores = scores, models = models, meta = meta)
}

dep_scores <- read_need(file.path(table_dir, "09_dependence_scores.csv")) %>%
  filter(
    .data$trait_label == "Type 2 diabetes",
    .data$trained_bucket == "Train: European-only",
    .data$target_ancestry %in% c("South Asian", "Hispanic or Latin American")
  )
if (nrow(dep_scores) != 24L) {
  stop("Expected 24 score rows for these two comparisons; got ", nrow(dep_scores), ".", call. = FALSE)
}

x_clip <- 0.15
x_lim <- c(-x_clip, x_clip)
x_breaks <- c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15)
slide_width <- 13.333
slide_height <- 7.5

scale_x_slide <- function() {
  scale_x_continuous(
    limits = x_lim,
    expand = c(0.01, 0.01),
    breaks = x_breaks,
    labels = function(x) formatC(x, format = "f", digits = 2)
  )
}

theme_lesson <- function(base = 13) {
  theme_minimal(base_size = base, base_family = "DejaVu Sans") +
    theme(
      plot.title = element_text(face = "bold", size = base + 2, colour = "#1a1a1a"),
      plot.subtitle = element_text(size = base - 1.5, colour = "#333333", lineheight = 1.15),
      plot.caption = element_text(size = base - 3.5, colour = "#444444", hjust = 0, lineheight = 1.15),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title = element_text(size = base - 1, colour = "#222222"),
      axis.text = element_text(colour = "#222222"),
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, size = base - 1, colour = "#1a1a1a"),
      strip.background = element_rect(fill = "#F3F3F3", colour = NA),
      plot.margin = margin(8, 14, 6, 10)
    )
}

wrap_cap <- function(text, width = 140) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  paste(vapply(lines, function(line) {
    paste(strwrap(line, width = width), collapse = "\n")
  }, character(1)), collapse = "\n")
}

fmt_n <- function(x) format(as.integer(round(x)), big.mark = ",")
fmt3 <- function(x) formatC(x, format = "f", digits = 3)
fmt4 <- function(x) formatC(x, format = "f", digits = 4)

save_slide <- function(plot, filename) {
  png_path <- file.path(root_dir, paste0(filename, ".png"))
  pdf_path <- file.path(root_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = slide_width, height = slide_height, dpi = 170, bg = "white", limitsize = FALSE)
  ggsave(
    pdf_path, plot,
    width = slide_width, height = slide_height, bg = "white",
    device = grDevices::cairo_pdf, limitsize = FALSE
  )
  message("Wrote ", png_path)
}

label_for_tgt_ids <- function(ids_txt, meta) {
  ids <- strsplit(ids_txt, ";", fixed = TRUE)[[1]]
  hit <- meta %>%
    filter(.data$sampleset_id %in% ids) %>%
    mutate(sampleset_id = factor(.data$sampleset_id, levels = ids)) %>%
    arrange(.data$sampleset_id)
  if (nrow(hit) != length(ids)) {
    stop("Missing study names for sample set(s): ", ids_txt, call. = FALSE)
  }
  paste(sprintf("%s, %s people", hit$study, fmt_n(hit$n_individuals)), collapse = " + ")
}

clip_intervals <- function(df) {
  df %>%
    mutate(
      lo_draw = pmax(.data$lo, x_lim[1]),
      hi_draw = pmin(.data$hi, x_lim[2]),
      estimate_draw = pmin(pmax(.data$estimate, x_lim[1]), x_lim[2]),
      arrow_left = .data$lo < x_lim[1],
      arrow_right = .data$hi > x_lim[2]
    )
}

summary_cols <- c(
  "avg" = "#222222",
  "lm" = "#7A7A7A",
  "lmer" = "#01665E",
  "rma" = "#3C6E9A"
)

# One comparison, two panels, shared x-axis.
forest_for_comparison <- function(ancestry_name, prepared, dep) {
  score_rows <- prepared$scores
  joined <- score_rows %>%
    inner_join(
      dep %>% select(pgs_id, delta_saved = delta, se_saved = delta_se, ids_saved = tgt_sampleset_ids),
      by = "pgs_id"
    )
  if (nrow(joined) != nrow(score_rows) || nrow(joined) != nrow(dep)) {
    stop("Per-score rows do not match 09_dependence_scores.csv for ", ancestry_name, ".", call. = FALSE)
  }
  if (
    max(abs(joined$delta - joined$delta_saved)) > 1e-8 ||
      max(abs(joined$delta_se - joined$se_saved)) > 1e-8 ||
      any(joined$tgt_sampleset_ids != joined$ids_saved)
  ) {
    stop("Per-score Δ AUC does not match 09_dependence_scores.csv for ", ancestry_name, ".", call. = FALSE)
  }

  meta <- prepared$meta
  score_plot <- score_rows %>%
    mutate(
      sample_set_label = vapply(.data$tgt_sampleset_ids, label_for_tgt_ids, character(1), meta = meta),
      lo = .data$delta - 1.96 * .data$delta_se,
      hi = .data$delta + 1.96 * .data$delta_se,
      estimate = .data$delta,
      axis = .data$pgs_id
    )
  n_sets <- n_distinct(score_plot$tgt_sampleset_ids)
  set_levels <- score_plot %>%
    distinct(.data$tgt_sampleset_ids, .data$sample_set_label) %>%
    arrange(.data$tgt_sampleset_ids) %>%
    pull(.data$sample_set_label)
  score_plot <- score_plot %>%
    mutate(
      sample_set_label = factor(.data$sample_set_label, levels = set_levels),
      axis = factor(.data$axis, levels = rev(.data$axis))
    )

  summary_plot <- prepared$models %>%
    mutate(axis = factor(.data$model, levels = rev(auc_model_names)))
  score_draw <- clip_intervals(score_plot)
  summary_draw <- clip_intervals(summary_plot)

  set_cols <- c("#332288", "#E69F00", "#009E73", "#CC79A7")[seq_len(n_sets)]
  names(set_cols) <- set_levels

  x_lab <- sprintf(
    "\u0394AUC (%s \u2212 European) with 95%% interval",
    ancestry_name
  )

  p_scores <- ggplot(score_draw, aes(y = .data$axis)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_errorbar(
      aes(xmin = .data$lo_draw, xmax = .data$hi_draw, colour = .data$sample_set_label),
      orientation = "y", width = 0.22, linewidth = 0.55
    ) +
    geom_point(aes(x = .data$estimate_draw, colour = .data$sample_set_label), size = 2.2) +
    geom_segment(
      data = score_draw %>% filter(.data$arrow_left),
      aes(
        x = x_lim[1] + 0.012, xend = x_lim[1],
        y = .data$axis, yend = .data$axis,
        colour = .data$sample_set_label
      ),
      arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed"),
      linewidth = 0.45, show.legend = FALSE
    ) +
    geom_segment(
      data = score_draw %>% filter(.data$arrow_right),
      aes(
        x = x_lim[2] - 0.012, xend = x_lim[2],
        y = .data$axis, yend = .data$axis,
        colour = .data$sample_set_label
      ),
      arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed"),
      linewidth = 0.45, show.legend = FALSE
    ) +
    scale_colour_manual(values = set_cols, drop = FALSE) +
    scale_x_slide() +
    coord_cartesian(xlim = x_lim, clip = "off") +
    labs(
      title = NULL,
      subtitle = sprintf(
        "%d scores, coloured by target sample set%s",
        nrow(score_plot),
        if (n_sets == 1L) "." else sprintf(" (%d in this comparison).", n_sets)
      ),
      x = NULL, y = NULL
    ) +
    guides(colour = guide_legend(nrow = if (n_sets > 1L) 2L else 1L, byrow = TRUE)) +
    theme_lesson(11) +
    theme(
      legend.text = element_text(size = 8),
      axis.text.y = element_text(size = 8.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(4, 14, 2, 10)
    )

  p_sum <- ggplot(summary_draw, aes(y = .data$axis)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_errorbar(
      aes(xmin = .data$lo_draw, xmax = .data$hi_draw, colour = .data$kind),
      orientation = "y", width = 0.25, linewidth = 0.6, show.legend = FALSE
    ) +
    geom_point(
      data = summary_draw %>% filter(!.data$hollow),
      aes(x = .data$estimate_draw, colour = .data$kind),
      size = 2.4, show.legend = FALSE
    ) +
    geom_point(
      data = summary_draw %>% filter(.data$hollow),
      aes(x = .data$estimate_draw, colour = .data$kind),
      size = 2.5, shape = 21, fill = "white", stroke = 1.05, show.legend = FALSE
    ) +
    geom_segment(
      data = summary_draw %>% filter(.data$arrow_left),
      aes(
        x = x_lim[1] + 0.012, xend = x_lim[1],
        y = .data$axis, yend = .data$axis,
        colour = .data$kind
      ),
      arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed"),
      linewidth = 0.5, show.legend = FALSE
    ) +
    geom_segment(
      data = summary_draw %>% filter(.data$arrow_right),
      aes(
        x = x_lim[2] - 0.012, xend = x_lim[2],
        y = .data$axis, yend = .data$axis,
        colour = .data$kind
      ),
      arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed"),
      linewidth = 0.5, show.legend = FALSE
    ) +
    scale_colour_manual(values = summary_cols) +
    scale_x_slide() +
    coord_cartesian(xlim = x_lim, clip = "off") +
    labs(
      title = NULL,
      subtitle = "AUC models from the report. Hollow points are ancestry coefficients. The filled point is the average of all scores.",
      x = x_lab, y = NULL
    ) +
    theme_lesson(11) +
    theme(
      axis.text.y = element_text(size = 7.4),
      plot.margin = margin(2, 14, 4, 10)
    )

  past <- summary_plot %>% filter(.data$lo < x_lim[1] | .data$hi > x_lim[2])
  past_line <- if (nrow(past) == 0L) {
    character()
  } else {
    paste0(
      "Intervals past the edge: ",
      paste(sprintf("%s (%s to %s)", past$model, fmt4(past$lo), fmt4(past$hi)), collapse = "; "),
      "."
    )
  }
  noted <- summary_plot %>% filter(nzchar(.data$note))
  note_line <- if (nrow(noted) == 0L) {
    character()
  } else {
    paste(sprintf("%s: %s.", noted$model, noted$note), collapse = " ")
  }
  caption_txt <- wrap_cap(paste(c(
    sprintf(
      "Catalog stamp %s. Shared x-axis clipped at \u00b1%s; arrows mark intervals that run past the edge.",
      stamp, fmt3(x_clip)
    ),
    "In (1|cohort), cohort means the sample set. Grey: lm. Green: lmer. Blue: rma.mv. The filled point is pipeline 1.",
    past_line,
    note_line
  ), collapse = "\n"))

  n_scores <- nrow(score_plot)
  (
    p_scores / p_sum +
      plot_layout(heights = c(max(n_scores, 8L), length(auc_model_names) + 1.4))
  ) +
    plot_annotation(
      title = sprintf(
        "Type 2 diabetes, %s participants, European-only training",
        ancestry_name
      ),
      caption = caption_txt,
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, family = "DejaVu Sans", colour = "#1a1a1a"),
        plot.caption = element_text(size = 8, family = "DejaVu Sans", colour = "#444444", hjust = 0, lineheight = 1.12)
      )
    )
}

targets <- c("South Asian", "Hispanic or Latin American")
prepared <- setNames(lapply(targets, prepare_comparison), targets)

p_sa <- forest_for_comparison(
  "South Asian",
  prepared[["South Asian"]],
  dep_scores %>% filter(.data$target_ancestry == "South Asian")
)
p_hi <- forest_for_comparison(
  "Hispanic or Latin American",
  prepared[["Hispanic or Latin American"]],
  dep_scores %>% filter(.data$target_ancestry == "Hispanic or Latin American")
)

save_slide(p_sa, "01_south_asian")
save_slide(p_hi, "02_hispanic")

message("Additional forests written to ", normalizePath(root_dir))
