#!/usr/bin/env Rscript

# ============================================================
# Script: 11_auc_level_scatter.R
# Project: pgscatalog_meta
# Purpose: Show the Figure 1C European and non-European AUCs on one figure.
#          Panel A plots the two AUCs. Panel B plots the difference against
#          the European AUC. Large points are the inverse-variance weighted
#          average of all scores in the comparison, with the Figure 1C weights.
# Inputs:  results/pgs_catalog_<STAMP>/stage2/deltaAUC_forest_faceted_pairs_<STAMP>.csv
#          results/review_sensitivity_<STAMP>/tables/05_figure1c_cells.csv
# Outputs: results/review_sensitivity_<STAMP>/plots/18_auc_level_scatter.pdf/.png
#          results/review_sensitivity_<STAMP>/tables/18_auc_level_scores.csv
#          results/review_sensitivity_<STAMP>/tables/18_auc_level_comparisons.csv
# Run after: 08_generate_publication_figures.R and 09_review_sensitivity.R
# ============================================================
# Axis titles and the estimate columns name the same three quantities:
#   european_auc       European AUC
#   non_european_auc   Non-European AUC
#   difference         Difference (non-European − European)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(glue)
  library(stringr)
  library(scales)
  library(patchwork)
})

if (!file.exists("R/load.R")) {
  stop("R/load.R not found. Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

stamp <- pipeline_stamp
root_out   <- file.path("results", paste0("review_sensitivity_", stamp))
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

pairs_path <- file.path(
  catalog_results_dir(stamp), "stage2",
  paste0("deltaAUC_forest_faceted_pairs_", stamp, ".csv")
)
cells_path <- file.path(dir_tables, "05_figure1c_cells.csv")
if (!file.exists(pairs_path)) stop("Pairs file not found: ", pairs_path, call. = FALSE)
if (!file.exists(cells_path)) stop("Figure 1C file not found: ", cells_path, call. = FALSE)

# ---------------------------------------------------------------------------
# Scores and the inverse-variance weighted average of all scores in each comparison
# ---------------------------------------------------------------------------

pairs <- read_csv(pairs_path, show_col_types = FALSE) %>%
  mutate(
    trait_label = as.character(trait_label),
    trained_bucket = as.character(trained_bucket),
    target_ancestry = as.character(target_ancestry),
    pgs_id = as.character(pgs_id),
    is_pooled = as.logical(is_pooled)
  )

if (anyDuplicated(pairs[c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")])) {
  stop("A score contributed more than one difference to a comparison.", call. = FALSE)
}

# eur_se and tgt_se in the pairs file are on the logit scale. Figure 1C
# converts them with se_auc = se_logit * AUC * (1 - AUC) before weighting.
scored <- pairs %>%
  mutate(
    eur_se_auc = auc_scale_se(eur_auc, eur_se, "delta_method"),
    tgt_se_auc = auc_scale_se(tgt_auc, tgt_se, "delta_method"),
    difference = tgt_auc - eur_auc,
    difference_se = sqrt(eur_se_auc^2 + tgt_se_auc^2)
  ) %>%
  filter(is.finite(difference), is.finite(difference_se), difference_se > 0)

pooled <- pool_stage2_cells(pairs, se_scale = "delta_method") %>%
  mutate(
    trait_label = as.character(trait_label),
    trained_bucket = as.character(trained_bucket),
    target_ancestry = as.character(target_ancestry)
  )

averages <- scored %>%
  group_by(trait_label, trained_bucket, target_ancestry) %>%
  summarise(
    n_scores = n(),
    european_auc = sum(eur_auc / difference_se^2) / sum(1 / difference_se^2),
    non_european_auc = sum(tgt_auc / difference_se^2) / sum(1 / difference_se^2),
    .groups = "drop"
  ) %>%
  mutate(difference = non_european_auc - european_auc) %>%
  inner_join(pooled, by = c("trait_label", "trained_bucket", "target_ancestry"))

fig1c <- read_csv(cells_path, show_col_types = FALSE) %>%
  transmute(
    trait_label = as.character(trait_label),
    trained_bucket = as.character(trained_bucket),
    target_ancestry = as.character(target_ancestry),
    n_pairs_fig = n_pairs,
    delta_hat_fig = delta_hat,
    lo_fig = lo,
    hi_fig = hi,
    ci_excludes_zero_fig = ci_excludes_zero
  )

comparisons <- averages %>%
  inner_join(fig1c, by = c("trait_label", "trained_bucket", "target_ancestry"))

max_abs <- function(x) max(abs(x), na.rm = TRUE)
tol <- 1e-8
gap_vs_pool <- max_abs(comparisons$difference - comparisons$delta_hat)
gap_vs_fig <- max_abs(comparisons$difference - comparisons$delta_hat_fig)
gap_lo <- max_abs(comparisons$lo - comparisons$lo_fig)
gap_hi <- max_abs(comparisons$hi - comparisons$hi_fig)

if (nrow(scored) != 127L) {
  stop(glue("Expected 127 scores; found {nrow(scored)}."), call. = FALSE)
}
if (nrow(comparisons) != 32L) {
  stop(glue("Expected 32 comparisons; found {nrow(comparisons)}."), call. = FALSE)
}
if (any(comparisons$n_scores != comparisons$n_pairs) || any(comparisons$n_scores != comparisons$n_pairs_fig)) {
  stop("Score counts disagree with Figure 1C.", call. = FALSE)
}
if (gap_vs_pool > tol || gap_vs_fig > tol || gap_lo > tol || gap_hi > tol) {
  stop(glue(
    "Comparison averages disagree with Figure 1C ",
    "(difference {gap_vs_fig}, interval {max(gap_lo, gap_hi)})."
  ), call. = FALSE)
}
if (sum((comparisons$lo > 0) | (comparisons$hi < 0)) != 21L) {
  stop("Expected 21 corrected intervals to exclude zero.", call. = FALSE)
}
if (!all(comparisons$ci_excludes_zero_fig == ((comparisons$lo > 0) | (comparisons$hi < 0)))) {
  stop("Interval flags disagree with Figure 1C.", call. = FALSE)
}

message(glue(
  ">> 127 scores, 32 comparisons. ",
  "Largest gap between the inverse-variance weighted average AUCs and the Figure 1C difference: {signif(gap_vs_fig, 3)}."
))

trait_display <- c(
  "Incident chronic obstructive pulmonary disease before age 50 years" = "Incident COPD before age 50",
  "Coronary heart disease (incident and prevalent)" = "Coronary heart disease (incident + prevalent)",
  "Incident coronary heart disease" = "Incident coronary heart disease",
  "Stroke excluding subarachnoid hemorrhage" = "Stroke (excl. subarachnoid hemorrhage)",
  "Epithelial non-mucinous ovarian cancer" = "Ovarian cancer (epithelial, non-mucinous)",
  "Apparent Treatment-Resistant Hypertension" = "Treatment-resistant hypertension"
)

short_trait <- function(x) {
  x <- as.character(x)
  shown <- unname(trait_display[x])
  ifelse(is.na(shown), x, shown)
}

# These four comparisons are the ones at or below a difference of -0.10.
# Anchors sit off the points. hjust 1 ends the text at the anchor; vjust 1
# puts the top of the text on the anchor, so a label below its point does
# not run back along the leader.
label_pos <- tibble::tribble(
  ~trait_label, ~label, ~x_a, ~y_a, ~h_a, ~v_a, ~x_b, ~y_b, ~h_b, ~v_b,
  "Colorectal cancer", "Colorectal\ncancer", 0.600, 0.488, 1, 1, 0.545, -0.185, 1, 0.5,
  "Breast cancer", "Breast cancer", 0.730, 0.488, 0, 1, 0.790, -0.185, 0, 0.5,
  "5-year incident prostate cancer", "5-year incident\nprostate cancer", 0.640, 0.760, 1, 0, 0.660, -0.185, 0.5, 0.5,
  "Incident chronic obstructive pulmonary disease before age 50 years", "Incident COPD\nbefore age 50", 0.640, 0.860, 1, 0, 0.815, -0.132, 0, 0.5
)

scores <- scored %>%
  transmute(
    trait_label,
    trained_bucket,
    target_ancestry,
    pgs_id,
    european_auc = eur_auc,
    non_european_auc = tgt_auc,
    difference
  )

comparisons <- comparisons %>%
  transmute(
    trait_label,
    trait_short = short_trait(trait_label),
    trained_bucket,
    target_ancestry,
    n_scores,
    european_auc,
    non_european_auc,
    difference,
    difference_lo = lo,
    difference_hi = hi,
    interval_excludes_zero = (lo > 0) | (hi < 0)
  ) %>%
  arrange(trait_label, target_ancestry, trained_bucket)

labels <- comparisons %>%
  filter(difference <= -0.10) %>%
  inner_join(label_pos, by = "trait_label")

if (nrow(labels) != 4L || any(labels$target_ancestry != "African")) {
  stop("Expected four African comparisons at or below a difference of -0.10.", call. = FALSE)
}
expected_labels <- c(
  "Breast cancer",
  "5-year incident prostate cancer",
  "Colorectal cancer",
  "Incident chronic obstructive pulmonary disease before age 50 years"
)
if (!setequal(labels$trait_label, expected_labels)) {
  stop("The comparisons at or below -0.10 are not the expected four.", call. = FALSE)
}
if (any(labels$trait_short[labels$trait_label == "Incident chronic obstructive pulmonary disease before age 50 years"] != "Incident COPD before age 50")) {
  stop("The COPD label does not use the short trait name.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

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

anc_levels <- intersect(names(pal_ancestry), unique(scores$target_ancestry))
scores <- scores %>% mutate(target_ancestry = factor(target_ancestry, levels = anc_levels))
comparisons <- comparisons %>% mutate(target_ancestry = factor(target_ancestry, levels = anc_levels))
labels <- labels %>% mutate(target_ancestry = factor(target_ancestry, levels = anc_levels))

# A little past the scores, so the labels below 0.50 stay inside the panel.
auc_limits <- c(0.42, 0.94)
auc_breaks <- seq(0.5, 0.9, by = 0.1)
legend_name <- "Non-European ancestry group"

theme_auc <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(4, 8, 4, 4)
    )
}

scale_auc_x <- function() {
  scale_x_continuous(
    limits = auc_limits,
    breaks = auc_breaks,
    labels = label_number(accuracy = 0.1),
    expand = expansion(0, 0)
  )
}

draw_scores <- function(y_var) {
  geom_point(
    data = scores,
    aes(x = european_auc, y = .data[[y_var]], colour = target_ancestry),
    size = 1.7, alpha = 0.4, show.legend = TRUE
  )
}

draw_averages <- function(y_var) {
  geom_point(
    data = comparisons,
    aes(x = european_auc, y = .data[[y_var]], fill = target_ancestry),
    shape = 21, colour = "black", size = 3.1, stroke = 0.45, show.legend = FALSE
  )
}

draw_labels <- function(x_col, y_col, h_col, v_col, y_point) {
  list(
    geom_segment(
      data = labels,
      aes(
        x = european_auc, y = .data[[y_point]],
        xend = .data[[x_col]], yend = .data[[y_col]]
      ),
      colour = "grey60", linewidth = 0.25, inherit.aes = FALSE
    ),
    geom_label(
      data = labels,
      aes(
        x = .data[[x_col]], y = .data[[y_col]],
        label = label, hjust = .data[[h_col]], vjust = .data[[v_col]]
      ),
      size = 2.5, lineheight = 0.9, colour = "grey10", fill = "white",
      linewidth = 0, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines"),
      inherit.aes = FALSE
    )
  )
}

p_a <- ggplot() +
  geom_hline(yintercept = 0.5, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0.5, colour = "grey70", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey25", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = -0.05, colour = "grey45", linewidth = 0.4, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -0.10, colour = "grey25", linewidth = 0.45, linetype = "dashed") +
  draw_scores("non_european_auc") +
  draw_averages("non_european_auc") +
  draw_labels("x_a", "y_a", "h_a", "v_a", "non_european_auc") +
  annotate(
    "text",
    x = 0.912,
    y = c(0.908, 0.858, 0.808),
    label = c("same AUC", "0.05 lower", "0.10 lower"),
    hjust = 1, vjust = 1, size = 2.5, colour = "grey25"
  ) +
  scale_auc_x() +
  scale_y_continuous(
    limits = auc_limits,
    breaks = auc_breaks,
    labels = label_number(accuracy = 0.1),
    expand = expansion(0, 0)
  ) +
  scale_colour_manual(values = pal_ancestry, breaks = anc_levels, name = legend_name) +
  scale_fill_manual(values = pal_ancestry, breaks = anc_levels, guide = "none") +
  coord_fixed(clip = "on") +
  labs(
    title = "Non-European AUC against European AUC",
    x = "European AUC",
    y = "Non-European AUC"
  ) +
  theme_auc() +
  guides(colour = guide_legend(nrow = 2, override.aes = list(alpha = 1, size = 3)))

p_b <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey25", linewidth = 0.5) +
  geom_hline(yintercept = -0.10, colour = "grey25", linewidth = 0.45, linetype = "dashed") +
  draw_scores("difference") +
  geom_errorbar(
    data = comparisons,
    aes(
      x = european_auc, ymin = difference_lo, ymax = difference_hi,
      colour = target_ancestry
    ),
    width = 0, linewidth = 0.4, alpha = 0.9, show.legend = FALSE
  ) +
  draw_averages("difference") +
  draw_labels("x_b", "y_b", "h_b", "v_b", "difference") +
  scale_auc_x() +
  scale_y_continuous(
    limits = c(-0.20, 0.20),
    breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
    labels = label_number(accuracy = 0.1, style_negative = "minus"),
    expand = expansion(0, 0)
  ) +
  scale_colour_manual(values = pal_ancestry, breaks = anc_levels, name = legend_name) +
  scale_fill_manual(values = pal_ancestry, breaks = anc_levels, guide = "none") +
  coord_cartesian(clip = "on") +
  labs(
    title = "Difference against European AUC",
    x = "European AUC",
    y = "Difference (non-European \u2212 European)"
  ) +
  theme_auc() +
  guides(colour = guide_legend(nrow = 2, override.aes = list(alpha = 1, size = 3)))

caption_txt <- paste(
  "127 scores and 32 comparisons. The difference is the non-European AUC minus the European AUC.",
  "Large points are the inverse-variance weighted average of all scores in the comparison.",
  "The weights are the Figure 1C weights, so the gap between the European and non-European inverse-variance weighted averages is the Figure 1C difference.",
  "Vertical intervals are the corrected 95% intervals for that inverse-variance weighted average. 21 of the 32 intervals exclude zero.",
  "On the left, the solid diagonal is equal AUC, the dashed lines are 0.05 lower and 0.10 lower, and the grey lines are AUC 0.50.",
  "On the right, the solid horizontal line is a difference of zero and the dashed line is \u22120.10. Labels mark the comparisons at or below a difference of \u22120.10.",
  sep = "\n"
)

p <- p_a + p_b +
  plot_layout(guides = "collect", widths = c(1, 1)) +
  plot_annotation(
    title = "European and non-European AUC",
    subtitle = "Small points are scores. Large points are the inverse-variance weighted average of all scores in the comparison.",
    caption = caption_txt,
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.05),
      plot.tag = element_text(face = "bold", size = 12)
    )
  ) &
  theme(legend.position = "bottom")

scores_out <- scores %>%
  arrange(trait_label, target_ancestry, trained_bucket, pgs_id) %>%
  mutate(target_ancestry = as.character(target_ancestry))

write_csv(scores_out, file.path(dir_tables, "18_auc_level_scores.csv"))
write_csv(comparisons %>% mutate(target_ancestry = as.character(target_ancestry)),
          file.path(dir_tables, "18_auc_level_comparisons.csv"))
save_plot("18_auc_level_scatter", p, width = 12, height = 8.4)

message(">> Wrote 18_auc_level_scores.csv, 18_auc_level_comparisons.csv, and 18_auc_level_scatter.")
