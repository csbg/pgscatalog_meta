#!/usr/bin/env Rscript

# Forest slides: per-score ΔAUC beside lm, both lmer fits, the average of all
# scores, and the covariance model. Plotting only; no model refits.
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
  library(patchwork)
  library(readr)
  library(stringr)
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
diabetes_dir <- file.path(base_dir, "diabetes_cohort_models")
root_dir <- file.path(base_dir, "additional_forests")
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

read_need <- function(path) {
  if (!file.exists(path)) {
    stop("Missing ", path, ".", call. = FALSE)
  }
  read_csv(path, show_col_types = FALSE)
}

scores <- read_need(file.path(table_dir, "09_dependence_scores.csv"))
cells <- read_need(file.path(table_dir, "09_dependence_cells.csv"))
coefs <- read_need(file.path(diabetes_dir, "diabetes_model_coefficients.csv"))
eval_rows <- read_need(file.path(diabetes_dir, "diabetes_evaluation_rows.csv"))

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

wrap_cap <- function(text, width = 118) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  paste(vapply(lines, function(line) {
    paste(strwrap(line, width = width), collapse = "\n")
  }, character(1)), collapse = "\n")
}

fmt_n <- function(x) format(as.integer(round(x)), big.mark = ",")
fmt3 <- function(x) formatC(x, format = "f", digits = 3)

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

subset_comp <- function(df, trait, bucket, ancestry_name) {
  df %>%
    filter(
      .data$trait_label == .env$trait,
      .data$trained_bucket == .env$bucket,
      .data$target_ancestry == .env$ancestry_name
    )
}

# Study name and number of people for each sample set in these comparisons.
sample_set_meta <- eval_rows %>%
  filter(
    .data$trait_label == "Type 2 diabetes",
    .data$train_bucket == "Train: European-only",
    .data$ancestry_display %in% c("South Asian", "Hispanic or Latin American")
  ) %>%
  distinct(
    .data$ancestry_display, .data$sampleset_id, .data$study_name, .data$n_individuals
  )

label_for_tgt_ids <- function(ids_txt, ancestry_name) {
  ids <- strsplit(ids_txt, ";", fixed = TRUE)[[1]]
  meta <- sample_set_meta %>%
    filter(
      .data$ancestry_display == .env$ancestry_name,
      .data$sampleset_id %in% .env$ids
    )
  if (nrow(meta) != length(ids)) {
    stop("Missing study names for sample set(s): ", ids_txt, call. = FALSE)
  }
  meta <- meta %>%
    mutate(sampleset_id = factor(.data$sampleset_id, levels = ids)) %>%
    arrange(.data$sampleset_id)
  paste(
    sprintf("%s, %s people", meta$study_name, fmt_n(meta$n_individuals)),
    collapse = " + "
  )
}

model_labels <- c(
  lm = "lm(AUC ~ ancestry + PGS, weights = 1/SE\u00b2)",
  lmer_fixed = "lmer(AUC ~ ancestry + PGS + (1|cohort))",
  lmer_random = "lmer(AUC ~ ancestry + (1|PGS) + (1|cohort))",
  avg = "Average of all scores",
  cov = "Covariance model (\u03c1 = 0.6)"
)

pick_coef <- function(ancestry_name, model_pattern, cohort_ok = NULL) {
  rows <- coefs %>%
    filter(
      .data$train_bucket == "Train: European-only",
      .data$ancestry == .env$ancestry_name,
      grepl(model_pattern, .data$model)
    )
  if (!is.null(cohort_ok)) {
    rows <- rows %>% filter(.data$cohort_def == .env$cohort_ok)
  }
  if (nrow(rows) != 1L) {
    stop(
      "Expected one coefficient row for ", ancestry_name, " / ", model_pattern,
      "; got ", nrow(rows), ".",
      call. = FALSE
    )
  }
  rows
}

clip_intervals <- function(df) {
  df %>%
    mutate(
      lo_draw = pmax(.data$lo, x_lim[1]),
      hi_draw = pmin(.data$hi, x_lim[2]),
      estimate_draw = pmin(pmax(.data$estimate, x_lim[1]), x_lim[2]),
      arrow_left = .data$lo < x_lim[1],
      arrow_right = .data$hi > x_lim[2],
      point_outside = .data$estimate < x_lim[1] | .data$estimate > x_lim[2]
    )
}

# One comparison → one two-panel slide.
forest_for_comparison <- function(trait, bucket, ancestry_name) {
  score_rows <- subset_comp(scores, trait, bucket, ancestry_name) %>%
    arrange(.data$tgt_sampleset_ids, .data$pgs_id)
  cell <- subset_comp(cells, trait, bucket, ancestry_name)
  if (nrow(cell) != 1L) {
    stop("Expected one comparison row for ", ancestry_name, ".", call. = FALSE)
  }
  if (!isTRUE(cell$goes_through_covariance_model)) {
    stop(ancestry_name, " does not go through the covariance model.", call. = FALSE)
  }
  if (isTRUE(cell$robust_reported)) {
    stop(
      ancestry_name, " unexpectedly reports a robust interval; this slide leaves robust out.",
      call. = FALSE
    )
  }

  lm_row <- pick_coef(ancestry_name, "^lm\\(AUC ~ ancestry \\+ PGS, weights")
  lmer_fixed <- pick_coef(
    ancestry_name,
    "^lmer\\(AUC ~ ancestry \\+ PGS \\+ \\(1\\|cohort\\)",
    "sample-set id"
  )
  lmer_random <- pick_coef(
    ancestry_name,
    "^lmer\\(AUC ~ ancestry \\+ \\(1\\|PGS\\) \\+ \\(1\\|cohort\\)",
    "sample-set id"
  )

  score_plot <- score_rows %>%
    mutate(
      sample_set_label = vapply(
        .data$tgt_sampleset_ids,
        label_for_tgt_ids,
        character(1),
        ancestry_name = ancestry_name
      ),
      lo = .data$delta - 1.96 * .data$delta_se,
      hi = .data$delta + 1.96 * .data$delta_se,
      estimate = .data$delta,
      axis = .data$pgs_id,
      panel = "Scores",
      kind = "score",
      hollow = FALSE
    )

  summary_plot <- tibble(
    axis = unname(model_labels[c("lm", "lmer_fixed", "lmer_random", "avg", "cov")]),
    estimate = c(
      lm_row$estimate,
      lmer_fixed$estimate,
      lmer_random$estimate,
      cell$delta_fixed,
      cell$delta_cov
    ),
    lo = c(lm_row$lo, lmer_fixed$lo, lmer_random$lo, cell$lo_fixed, cell$lo_cov),
    hi = c(lm_row$hi, lmer_fixed$hi, lmer_random$hi, cell$hi_fixed, cell$hi_cov),
    sample_set_label = NA_character_,
    panel = "Summaries",
    kind = c("lm", "lmer", "lmer", "avg", "cov"),
    hollow = c(FALSE, TRUE, TRUE, FALSE, FALSE)
  )

  # Drawn covariance values must match the saved comparison table.
  cov_drawn <- summary_plot %>% filter(.data$kind == "cov")
  if (
    abs(cov_drawn$estimate - cell$delta_cov) > 1e-8 ||
      abs(cov_drawn$lo - cell$lo_cov) > 1e-8 ||
      abs(cov_drawn$hi - cell$hi_cov) > 1e-8
  ) {
    stop(
      "Drawn covariance values do not match delta_cov, lo_cov, hi_cov in 09_dependence_cells.csv for ",
      ancestry_name, ".",
      call. = FALSE
    )
  }

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
  summary_plot <- summary_plot %>%
    mutate(axis = factor(.data$axis, levels = rev(.data$axis)))

  score_draw <- clip_intervals(score_plot)
  summary_draw <- clip_intervals(summary_plot)

  set_cols <- c("#332288", "#E69F00", "#009E73", "#CC79A7")[seq_len(n_sets)]
  names(set_cols) <- set_levels

  summary_cols <- c(
    "lm" = "#7A7A7A",
    "lmer" = "#01665E",
    "avg" = "#222222",
    "cov" = "#8C2D04"
  )

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
    geom_point(
      aes(x = .data$estimate_draw, colour = .data$sample_set_label),
      size = 2.2
    ) +
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
      x = NULL,
      y = NULL
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
      orientation = "y", width = 0.28, linewidth = 0.65, show.legend = FALSE
    ) +
    geom_point(
      data = summary_draw %>% filter(!.data$hollow),
      aes(x = .data$estimate_draw, colour = .data$kind),
      size = 2.6, show.legend = FALSE
    ) +
    geom_point(
      data = summary_draw %>% filter(.data$hollow),
      aes(x = .data$estimate_draw, colour = .data$kind),
      size = 2.8, shape = 21, fill = "white", stroke = 1.1, show.legend = FALSE
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
      subtitle = "Summaries. Hollow points are joint-model coefficients, not averages of the scores above.",
      x = x_lab,
      y = NULL
    ) +
    theme_lesson(11) +
    theme(
      axis.text.y = element_text(size = 8.2),
      plot.margin = margin(2, 14, 4, 10)
    )

  title_txt <- sprintf(
    "Type 2 diabetes, %s participants, European-only training",
    ancestry_name
  )
  caption_txt <- wrap_cap(paste(
    sprintf(
      "Catalog stamp %s. Shared x-axis clipped at \u00b1%s; arrows mark intervals that run past the edge.",
      stamp, fmt3(x_clip)
    ),
    "In (1|cohort), cohort means the sample set. The two lmer rows are ancestry coefficients from Nikolaus's joint models.",
    sprintf(
      "Covariance model: %s (%s to %s) at \u03c1 = 0.6. Neither comparison reports a robust interval: each is one sharing group, so degrees of freedom are too few.",
      fmt3(cell$delta_cov), fmt3(cell$lo_cov), fmt3(cell$hi_cov)
    ),
    sep = "\n"
  ))

  n_scores <- nrow(score_plot)
  (
    p_scores / p_sum +
      plot_layout(heights = c(max(n_scores, 8L), 5.2))
  ) +
    plot_annotation(
      title = title_txt,
      caption = caption_txt,
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, family = "DejaVu Sans", colour = "#1a1a1a"),
        plot.caption = element_text(size = 8.5, family = "DejaVu Sans", colour = "#444444", hjust = 0, lineheight = 1.15)
      )
    )
}

# -----------------------------------------------------------------
# Both comparisons share the same drawing function and the same x-axis.
# -----------------------------------------------------------------

two <- scores %>%
  filter(
    .data$trait_label == "Type 2 diabetes",
    .data$trained_bucket == "Train: European-only",
    .data$target_ancestry %in% c("South Asian", "Hispanic or Latin American")
  )
if (nrow(two) != 24L) {
  stop("Expected 24 score rows for these two comparisons; got ", nrow(two), ".", call. = FALSE)
}

p_sa <- forest_for_comparison(
  "Type 2 diabetes", "Train: European-only", "South Asian"
)
p_hi <- forest_for_comparison(
  "Type 2 diabetes", "Train: European-only", "Hispanic or Latin American"
)

save_slide(p_sa, "01_south_asian")
save_slide(p_hi, "02_hispanic")

message("Additional forests written to ", normalizePath(root_dir))
