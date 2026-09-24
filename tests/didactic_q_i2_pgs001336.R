#!/usr/bin/env Rscript

# ============================================================
# Didactic walkthrough of Cochran's Q and I²
# Two European breast-cancer cells, each with k = 2:
#   PGS001336  Q below df, so I² = 0
#   PGS000004  Q far above df, so I² is about 99.5
#
# Run from the repository root:
#   Rscript tests/didactic_q_i2_pgs001336.R
#
# Plots land in
#   results/pgs_catalog_<stamp>/q_i2_pgs001336/<PGS id>/
#
# The arithmetic uses only qlogis(), plogis(), and sums.
# It is the same calculation as ivw_pool_logit(), written out
# so each piece can be plotted.
#
# Q is the weighted scatter of the study estimates around their
# common value. It leaves the pooled AUC where the weights put it.
#
# If the studies share one true value, that scatter is a chi-square
# draw with df = k - 1. Its average is df: some disagreement is
# what sampling error looks like.
#
# I² = max(0, (Q - df) / Q) * 100
# The subtraction removes that chance benchmark. The division turns
# the leftover into a percent of Q. A negative leftover is reported
# as 0, because a percent of scatter cannot be negative.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tibble)
})

if (!file.exists("R/load.R")) {
  stop("Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

stamp <- pipeline_stamp
root_dir <- file.path(catalog_results_dir(stamp), "q_i2_pgs001336")
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_w <- function(x) format(round(x), big.mark = ",")
fmt_n <- function(x) format(as.integer(round(x)), big.mark = ",")

theme_lesson <- function(base = 13) {
  theme_minimal(base_size = base, base_family = "DejaVu Sans") +
    theme(
      plot.title = element_text(face = "bold", size = base + 2, colour = "#1a1a1a"),
      plot.subtitle = element_text(size = base - 1.5, colour = "#333333", lineheight = 1.15),
      plot.caption = element_text(
        size = base - 3.5, colour = "#444444", hjust = 0, lineheight = 1.15
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title = element_text(size = base - 1, colour = "#222222"),
      axis.text = element_text(colour = "#222222"),
      legend.position = "none",
      plot.margin = margin(12, 18, 12, 12)
    )
}

ann_theme <- function() {
  theme(
    plot.title = element_text(face = "bold", size = 15, family = "DejaVu Sans", colour = "#1a1a1a"),
    plot.subtitle = element_text(size = 11, family = "DejaVu Sans", colour = "#333333", lineheight = 1.15),
    plot.caption = element_text(size = 9.5, family = "DejaVu Sans", colour = "#444444", hjust = 0, lineheight = 1.15),
    plot.margin = margin(8, 12, 8, 12)
  )
}

bracket_path <- function(xmin, xmax, y, height = 0.08) {
  tibble(
    x = c(xmin, xmin, xmax, xmax),
    y = c(y - height, y, y, y - height)
  )
}

# Named by performance id. arrange(auc) later sets the drawing order.
examples <- list(
  list(
    pgs_id = "PGS001336",
    trait = "Breast cancer",
    ancestry = "European",
    style = "close",
    shorts = c(PPM009168 = "Non-white British", PPM009170 = "White British")
  ),
  list(
    pgs_id = "PGS000004",
    trait = "Breast cancer",
    ancestry = "European",
    style = "split",
    shorts = c(PPM005168 = "eMERGE", PPM038141 = "UCLA")
  )
)

eval_path <- file.path(
  "results", paste0("pgs_auc_ci_audit_", stamp), "eval_df_final_auc_ci.csv"
)
sample_path <- catalog_bulk_file("evaluation_sample_sets", stamp)
eval_df <- read_csv(eval_path, show_col_types = FALSE)
samples <- read_csv(sample_path, show_col_types = FALSE) %>%
  transmute(
    sampleset_id = `PGS Sample Set (PSS)`,
    n_individuals = `Number of Individuals`,
    n_cases = `Number of Cases`,
    cohort_note = `Additional Ancestry Description`,
    cohort_name = `Cohort(s)`
  ) %>%
  distinct(sampleset_id, .keep_all = TRUE)

published_all <- read_csv(
  file.path(catalog_results_dir(stamp), "stage1", paste0("stage1_pooled_cells_", stamp, ".csv")),
  show_col_types = FALSE
)

render_lesson <- function(ex) {
  out_dir <- file.path(root_dir, ex$pgs_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  save_lesson <- function(plot, filename, width, height) {
    png_path <- file.path(out_dir, paste0(filename, ".png"))
    pdf_path <- file.path(out_dir, paste0(filename, ".pdf"))
    ggsave(png_path, plot, width = width, height = height, dpi = 170, bg = "white")
    ggsave(
      pdf_path, plot, width = width, height = height, bg = "white",
      device = grDevices::cairo_pdf
    )
    message("Wrote ", png_path)
  }

  studies <- eval_df %>%
    filter(performance_id %in% names(ex$shorts)) %>%
    mutate(short = unname(ex$shorts[performance_id])) %>%
    left_join(samples, by = "sampleset_id") %>%
    arrange(auc)

  if (nrow(studies) != 2L) {
    stop("Expected two rows for ", ex$pgs_id, ", found ", nrow(studies), call. = FALSE)
  }

  studies <- studies %>%
    mutate(
      axis = sprintf(
        "%s\n%s people, %s cases",
        short, fmt_n(n_individuals), fmt_n(n_cases)
      )
    )

  # Logit AUC. SE = (logit(upper) - logit(lower)) / (2 * 1.96)
  auc <- studies$auc
  lo <- studies$estimate_ci_lower
  hi <- studies$estimate_ci_upper
  eta <- qlogis(auc)
  eta_lo <- qlogis(lo)
  eta_hi <- qlogis(hi)
  se <- (eta_hi - eta_lo) / (2 * 1.96)

  w <- 1 / se^2
  sw <- sum(w)
  eta_hat <- sum(w * eta) / sw
  se_hat <- sqrt(1 / sw)
  distance <- eta - eta_hat
  pull <- w * distance
  q_piece <- w * distance^2
  Q <- sum(q_piece)
  k <- length(eta)
  df <- k - 1
  excess_raw <- (Q - df) / Q
  I2 <- max(0, excess_raw) * 100
  auc_hat <- plogis(eta_hat)
  auc_lo <- plogis(eta_hat - 1.96 * se_hat)
  auc_hi <- plogis(eta_hat + 1.96 * se_hat)
  p_het <- pchisq(Q, df = df, lower.tail = FALSE)
  p_txt <- if (p_het < 0.001) "< 0.001" else fmt(p_het, 2)
  close <- identical(ex$style, "close")

  pipeline <- ivw_pool_logit(eta, se)
  stopifnot(abs(pipeline$Q - Q) < 1e-8)
  stopifnot(abs(pipeline$I2 - I2) < 1e-8)
  stopifnot(abs(pipeline$eta - eta_hat) < 1e-8)

  published <- published_all %>%
    filter(
      pgs_id == ex$pgs_id,
      ancestry_display == ex$ancestry,
      trait_label == ex$trait
    )
  stopifnot(nrow(published) == 1L)
  stopifnot(abs(published$Q - Q) < 1e-6)
  stopifnot(abs(published$I2 - I2) < 1e-6)

  worked <- studies %>%
    transmute(
      pgs_id, performance_id, sampleset_id,
      cohort = short, cohort_note, cohort_name,
      n_individuals, n_cases, auc,
      ci_lower = estimate_ci_lower, ci_upper = estimate_ci_upper,
      eta, eta_ci_lower = eta_lo, eta_ci_upper = eta_hi, se,
      weight = w, weight_share = w / sw,
      distance_from_pooled = distance,
      weight_times_distance = pull,
      Q_contribution = q_piece
    )
  summary_row <- tibble(
    pgs_id = ex$pgs_id,
    trait = ex$trait,
    ancestry = ex$ancestry,
    k = k,
    eta_hat = eta_hat,
    se_hat = se_hat,
    auc_pooled = auc_hat,
    auc_lo = auc_lo,
    auc_hi = auc_hi,
    Q = Q,
    df = df,
    expected_Q_if_only_chance = df,
    Q_minus_df = Q - df,
    raw_fraction = excess_raw,
    I2 = I2,
    heterogeneity_p = p_het
  )
  write_csv(worked, file.path(out_dir, paste0("worked_studies_", ex$pgs_id, ".csv")))
  write_csv(summary_row, file.path(out_dir, paste0("worked_summary_", ex$pgs_id, ".csv")))

  col_study <- setNames(c("#0072B2", "#E69F00"), studies$short)
  col_pooled <- "#222222"
  col_chance <- "#6B6B6B"
  col_excess <- "#D55E00"
  heavier <- studies$short[which.max(w)]
  heavier_share <- max(w) / sw

  forest_rows <- tibble(
    short = c(studies$short, "Pooled"),
    axis = c(studies$axis, "Pooled"),
    kind = c("study", "study", "pooled"),
    estimate = c(auc, auc_hat),
    lo = c(lo, auc_lo),
    hi = c(hi, auc_hi),
    colour = c(unname(col_study[studies$short]), col_pooled)
  )
  forest_rows$axis <- factor(forest_rows$axis, levels = rev(forest_rows$axis))

  # -----------------------------------------------------------------
  # 1. The two AUCs
  # -----------------------------------------------------------------
  if (close) {
    p1_title <- "1.  The two European evaluations of PGS001336"
    p1_sub <- "Breast cancer. Stage 1 compares evaluations that share a score and an ancestry.\nEast Asian and South Asian each have one evaluation of this score, so they have no Q."
    p1_cap <- "Points are the reported AUCs. The intervals are real, but at this scale they are thinner than the points.\nThe dotted line is AUC = 0.5, a score with no discrimination. Both studies, and the pooled diamond, sit near 0.81."
    zoom_xlim <- c(0.788, 0.862)
    zoom_sub <- "Each study interval covers the other study's point estimate. That is the pattern Q will call agreement."
  } else {
    p1_title <- "1.  The two European evaluations of PGS000004"
    p1_sub <- "Breast cancer. Stage 1 pools evaluations that share a score and an ancestry.\nThese two intervals sit apart: 0.59 to 0.61, and 0.69 to 0.71."
    p1_cap <- "Points are the reported AUCs. Whiskers are the reported 95% intervals. The diamond is the inverse-variance pooled AUC.\nThe dotted line is AUC = 0.5. The gap between the studies is large next to either interval, which is what a large Q looks like."
    zoom_xlim <- c(min(lo) - 0.025, max(hi) + 0.085)
    zoom_sub <- "Neither interval reaches the other study's AUC. Q will call this disagreement."
  }

  p1 <- ggplot(forest_rows, aes(x = estimate, y = axis, colour = colour)) +
    geom_vline(xintercept = 0.5, linetype = "dotted", colour = "#bbbbbb", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.12, linewidth = 0.9) +
    geom_point(data = function(d) dplyr::filter(d, kind == "study"), size = 3.2) +
    geom_point(data = function(d) dplyr::filter(d, kind == "pooled"), shape = 18, size = 5) +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0.45, 0.90), breaks = seq(0.5, 0.9, by = 0.1), name = "AUC") +
    scale_y_discrete(name = NULL) +
    labs(title = p1_title, subtitle = p1_sub, caption = p1_cap) +
    theme_lesson() +
    theme(panel.grid.major.x = element_line(colour = "#eeeeee"))

  p1_zoom <- ggplot(forest_rows, aes(x = estimate, y = axis, colour = colour)) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.12, linewidth = 0.9) +
    geom_point(data = function(d) dplyr::filter(d, kind == "study"), size = 3.2) +
    geom_point(data = function(d) dplyr::filter(d, kind == "pooled"), shape = 18, size = 5) +
    geom_text(
      aes(x = hi, label = sprintf("%s to %s", fmt(lo, 3), fmt(hi, 3))),
      hjust = -0.08, size = 3.3, colour = "#333333", family = "DejaVu Sans"
    ) +
    scale_colour_identity() +
    scale_x_continuous(limits = zoom_xlim, name = "AUC, zoomed") +
    scale_y_discrete(name = NULL) +
    labs(title = "The same three intervals, zoomed", subtitle = zoom_sub) +
    theme_lesson()

  save_lesson(p1 / p1_zoom + plot_layout(heights = c(1.05, 1)), "01_two_evaluations", 10.2, 7.4)

  # -----------------------------------------------------------------
  # 2. Logit scale
  # -----------------------------------------------------------------
  logit_rows <- tibble(
    short = studies$short,
    axis = factor(studies$axis, levels = rev(studies$axis)),
    eta = eta, eta_lo = eta_lo, eta_hi = eta_hi,
    colour = unname(col_study[studies$short])
  )
  if (close) {
    logit_xlim <- c(1.32, 1.62)
    pool_x <- 1.505
    pool_hjust <- 0
  } else {
    logit_xlim <- c(min(eta_lo) - 0.06, max(eta_hi) + 0.18)
    pool_x <- eta_hat + 0.03
    pool_hjust <- 0
  }

  p2_forest <- ggplot(logit_rows, aes(x = eta, y = axis, colour = colour)) +
    geom_vline(xintercept = eta_hat, linetype = "dashed", colour = col_pooled, linewidth = 0.5) +
    geom_errorbar(aes(xmin = eta_lo, xmax = eta_hi), orientation = "y", width = 0.15, linewidth = 0.9) +
    geom_point(size = 3.4) +
    annotate(
      "text", x = pool_x, y = 1.5, hjust = pool_hjust,
      label = sprintf("pooled logit\n%s", fmt(eta_hat, 3)),
      family = "DejaVu Sans", size = 3.3, colour = col_pooled, lineheight = 0.95
    ) +
    scale_colour_identity() +
    scale_x_continuous(name = "logit(AUC)", limits = logit_xlim) +
    scale_y_discrete(name = NULL) +
    labs(title = "Both studies on the logit scale") +
    theme_lesson(12) +
    coord_cartesian(clip = "off")

  heavier_line <- if (close) {
    c(
      "The larger study has the narrower interval,",
      "so its SE is smaller. That is the only reason",
      "it will receive more weight."
    )
  } else {
    c(
      sprintf("%s has the smaller SE,", heavier),
      "so it receives more of the weight.",
      "The logit, not the head count, sets that weight."
    )
  }
  arith_lines <- c(
    "Basic functions, one study at a time",
    "",
    studies$short[1],
    sprintf("  AUC = %s", fmt(auc[1], 5)),
    sprintf("  eta = log(AUC / (1 - AUC)) = %s", fmt(eta[1], 3)),
    sprintf("  logit interval = %s to %s", fmt(eta_lo[1], 3), fmt(eta_hi[1], 3)),
    sprintf("  SE = (upper - lower) / (2 x 1.96) = %s", fmt(se[1], 4)),
    "",
    studies$short[2],
    sprintf("  AUC = %s", fmt(auc[2], 5)),
    sprintf("  eta = log(AUC / (1 - AUC)) = %s", fmt(eta[2], 3)),
    sprintf("  logit interval = %s to %s", fmt(eta_lo[2], 3), fmt(eta_hi[2], 3)),
    sprintf("  SE = (upper - lower) / (2 x 1.96) = %s", fmt(se[2], 4)),
    "",
    heavier_line
  )
  p2_text <- ggplot() +
    annotate(
      "text", x = 0, y = rev(seq_along(arith_lines)), label = arith_lines,
      hjust = 0, vjust = 0.5, family = "DejaVu Sans Mono", size = 3.15, colour = "#1a1a1a"
    ) +
    xlim(-0.02, 1) +
    ylim(0.3, length(arith_lines) + 0.7) +
    theme_void(base_family = "DejaVu Sans") +
    theme(plot.margin = margin(28, 8, 8, 8))

  p2 <- p2_forest + p2_text +
    plot_layout(widths = c(1.05, 1.15)) +
    plot_annotation(
      title = "2.  Q is computed on the logit of AUC, not on AUC itself",
      subtitle = "A difference of 0.01 in AUC is not the same precision everywhere between 0 and 1.\nThe logit is the scale on which the weights, the mean, and Q are calculated.",
      caption = "The dashed line is the weighted mean of the two logits. Mapping it back with plogis() returns the pooled AUC in plot 1.",
      theme = ann_theme()
    )
  save_lesson(p2, "02_logit_scale", 12.2, 6.2)

  # -----------------------------------------------------------------
  # 3. Weights and the balance point
  # -----------------------------------------------------------------
  weight_df <- tibble(
    short = factor(studies$short, levels = rev(studies$short)),
    weight = w,
    share = w / sw,
    colour = unname(col_study[studies$short])
  )
  p3_w <- ggplot(weight_df, aes(x = weight, y = short, fill = colour)) +
    geom_col(width = 0.62, colour = NA) +
    geom_text(
      aes(label = sprintf("%s   (%s%%)", fmt_w(weight), fmt(100 * share, 0))),
      hjust = -0.06, family = "DejaVu Sans", size = 3.6
    ) +
    scale_fill_identity() +
    scale_x_continuous(
      name = "Weight  =  1 / SE²",
      limits = c(0, max(w) * 1.55),
      expand = c(0, 0)
    ) +
    scale_y_discrete(name = NULL) +
    labs(title = "Weight of each study") +
    theme_lesson(12) +
    theme(panel.grid.major.x = element_line(colour = "#eeeeee"))

  if (close) {
    product_x <- c(eta_hat - 0.010, eta[2] + 0.006)
    product_hjust <- c(1, 0)
    bal_xlim <- c(1.400, 1.530)
    p3_sub <- "Weight = 1 / SE². The white British evaluation is larger and has the tighter interval, so it holds about 73% of the weight.\nThe common value moves toward that study until weight times distance is the same on both sides."
  } else {
    product_x <- c(eta_hat - 0.08, eta[2] + 0.03)
    product_hjust <- c(1, 0)
    bal_xlim <- c(min(eta) - 0.22, max(eta) + 0.32)
    p3_sub <- sprintf(
      "Weight = 1 / SE². %s has the smaller SE, so it holds about %s%% of the weight.\nThe common value moves until weight times distance is the same on both sides.",
      heavier, fmt(100 * heavier_share, 0)
    )
  }
  balance_df <- tibble(
    short = studies$short,
    y = c(3, 1),
    eta = eta,
    colour = unname(col_study[studies$short]),
    product = sprintf(
      "%s\n%s  ×  %s  =  %s",
      studies$short, fmt_w(w), fmt(abs(distance), 4), fmt(abs(pull), 2)
    ),
    product_y = if (close) c(3.55, 1) else c(3.72, 1),
    product_x = product_x,
    product_hjust = product_hjust
  )
  p3_bal <- ggplot(balance_df) +
    geom_vline(xintercept = eta_hat, linetype = "dashed", colour = col_pooled, linewidth = 0.6) +
    geom_segment(
      aes(x = eta, xend = eta_hat, y = y, yend = y, colour = colour),
      linewidth = 1.05,
      arrow = arrow(length = unit(2.4, "mm"), type = "closed")
    ) +
    geom_point(aes(x = eta, y = y, colour = colour), size = 3.6) +
    annotate("point", x = eta_hat, y = 2, shape = 18, size = 5, colour = col_pooled) +
    annotate(
      "text", x = eta_hat + if (close) 0.003 else 0.02, y = 2.22, hjust = 0,
      label = "pooled", family = "DejaVu Sans", size = 3.3, colour = col_pooled
    ) +
    geom_text(
      aes(x = product_x, y = product_y, label = product, colour = colour, hjust = product_hjust),
      family = "DejaVu Sans", size = 3.3, lineheight = 1.05
    ) +
    scale_colour_identity() +
    scale_x_continuous(name = "logit(AUC)", limits = bal_xlim) +
    scale_y_continuous(limits = c(0.35, 4.05), name = NULL, breaks = NULL) +
    labs(
      title = "The pooled logit is where the pulls cancel",
      subtitle = sprintf(
        "Both products equal %s. The heavier study sits closer to the dashed line.",
        fmt(abs(pull[1]), 2)
      )
    ) +
    theme_lesson(12) +
    theme(panel.grid.major = element_blank())

  p3 <- p3_w / p3_bal +
    plot_layout(heights = c(0.72, 1)) +
    plot_annotation(
      title = "3.  Inverse-variance weighting is a balance, not a vote",
      subtitle = p3_sub,
      caption = sprintf(
        "Pooled logit = (w1 × eta1 + w2 × eta2) / (w1 + w2) = %s.    SE of that mean = 1 / sqrt(sum of weights) = %s.    AUC = plogis(%s) = %s.",
        fmt(eta_hat, 3), fmt(se_hat, 4), fmt(eta_hat, 3), fmt(auc_hat, 3)
      ),
      theme = ann_theme()
    )
  save_lesson(p3, "03_weights_and_balance", 11.2, 7.6)

  # -----------------------------------------------------------------
  # 4. Pieces of Q
  # -----------------------------------------------------------------
  piece_df <- tibble(
    short = factor(studies$short, levels = studies$short),
    q_piece = q_piece,
    colour = unname(col_study[studies$short]),
    formula = sprintf(
      "%s  ×  (%s)²   =   %s",
      fmt_w(w), fmt(distance, 4), fmt(q_piece, 3)
    )
  )
  q_xlim <- if (close) max(q_piece) * 2.15 else max(q_piece) * 3.3
  sum_x <- if (close) max(q_piece) * 0.04 else max(q_piece) * 1.65
  p4_cap <- if (close) {
    sprintf(
      "Q = %s. It describes scatter around the pooled logit. It does not replace that logit, and the pooled AUC stays %s.\nThe less precise study contributes more (%s versus %s) because it sits further from the mean, and squaring magnifies that gap.",
      fmt(Q, 3), fmt(auc_hat, 3), fmt(q_piece[1], 3), fmt(q_piece[2], 3)
    )
  } else {
    sprintf(
      "Q = %s. It describes scatter around the pooled logit. The pooled AUC stays %s.\n%s contributes %s and %s contributes %s. Squaring gives the study that sits further from the mean the larger piece.",
      fmt(Q, 3), fmt(auc_hat, 3),
      studies$short[1], fmt(q_piece[1], 3),
      studies$short[2], fmt(q_piece[2], 3)
    )
  }
  p4 <- ggplot(piece_df, aes(x = q_piece, y = short, fill = colour)) +
    geom_col(width = 0.58, colour = NA) +
    geom_text(aes(label = formula), hjust = -0.04, family = "DejaVu Sans Mono", size = 3.5) +
    annotate(
      "text", x = sum_x, y = 1.5, hjust = 0,
      family = "DejaVu Sans", size = 4.1, fontface = "bold", colour = "#1a1a1a",
      label = sprintf("Q  =  %s  +  %s  =  %s", fmt(q_piece[1], 3), fmt(q_piece[2], 3), fmt(Q, 3))
    ) +
    scale_fill_identity() +
    scale_x_continuous(
      name = "Contribution to Cochran's Q     w × (logit − pooled logit)²",
      limits = c(0, q_xlim), expand = c(0, 0)
    ) +
    scale_y_discrete(name = NULL, expand = expansion(add = 0.55)) +
    labs(
      title = "4.  Cochran's Q is the sum of the weighted squared distances",
      subtitle = "Square each distance from the balance point in plot 3, so the sign drops out.\nMultiply by that study's weight, then add the two products. That sum is Q.",
      caption = p4_cap
    ) +
    theme_lesson() +
    theme(panel.grid.major.x = element_line(colour = "#eeeeee"))
  save_lesson(p4, "04_building_Q", 11.2, 5.6)

  # -----------------------------------------------------------------
  # 5. Q against the chi-square benchmark
  # -----------------------------------------------------------------
  y_top <- 1.15
  x_panel <- 0.11
  x_right <- if (close) 6.2 else max(8, Q * 1.08)
  xs <- seq(x_panel, x_right, length.out = 800)
  chisq_df <- tibble(x = xs, density = dchisq(xs, df = df))
  if (close) {
    shade_df <- tibble(
      x = seq(Q, x_right, length.out = 400),
      density = dchisq(seq(Q, x_right, length.out = 400), df = df)
    )
    p5_label_x <- 2.15
    p5_label <- sprintf(
      "Blue line: observed Q = %s\nDashed line: expected Q if the studies\nshare one AUC,  df = k - 1 = 1\nShaded area: P(chi-square >= Q) = %s",
      fmt(Q, 3), p_txt
    )
    p5_sub <- "Two studies give one contrast, so df = 1. Chance alone produces a Q whose average is 1, not 0.\nThe median of this curve is 0.45. The observed Q, 0.44, sits on that median."
    p5_cap <- "The heterogeneity p-value is the shaded area: how often chance alone would produce a Q at least this large. Here, about half the time.\nThe density rises without bound as Q approaches 0. The panel is cut at 1.15 so the rest of the curve stays visible."
    x_breaks <- 0:6
  } else {
    shade_df <- tibble(x = numeric(0), density = numeric(0))
    p5_label_x <- 25
    p5_label <- sprintf(
      "Dashed line: expected Q if the studies\nshare one AUC,  df = k - 1 = 1\nBlue line: observed Q = %s\nP(chi-square >= Q) %s",
      fmt(Q, 1), p_txt
    )
    p5_sub <- "Two studies give one contrast, so df = 1. Chance alone produces a Q whose average is 1, not 0.\nThe observed Q is 187, far to the right of that benchmark."
    p5_cap <- "The chi-square curve is the spike at the left. Almost all of its area sits between 0 and 5.\nA Q of 187 is not a chance fluctuation of two studies that share one AUC. The p-value is below 0.001."
    x_breaks <- c(0, 50, 100, 150, 187)
  }

  p5 <- ggplot(chisq_df, aes(x = x, y = density)) +
    { if (nrow(shade_df)) geom_ribbon(
      data = shade_df, aes(x = x, ymin = 0, ymax = density),
      fill = "#0072B2", alpha = 0.15, inherit.aes = FALSE
    ) } +
    geom_line(linewidth = 0.8, colour = "#333333") +
    geom_vline(xintercept = Q, colour = "#0072B2", linewidth = 0.9) +
    geom_vline(xintercept = df, colour = col_chance, linewidth = 0.8, linetype = "dashed") +
    annotate(
      "label", x = p5_label_x, y = 1.05, hjust = 0, vjust = 1,
      fill = "white", linewidth = 0.25, colour = "#1a1a1a",
      family = "DejaVu Sans", size = 3.5, lineheight = 1.05,
      label = p5_label
    ) +
    scale_x_continuous(
      name = "Cochran's Q, and the chi-square distribution it follows when the studies share one AUC",
      limits = c(0, x_right), breaks = x_breaks, expand = c(0, 0)
    ) +
    scale_y_continuous(name = "Chi-square density (df = 1)", limits = c(0, y_top), expand = c(0, 0)) +
    labs(
      title = "5.  Under pure sampling error, Q averages df, not zero",
      subtitle = p5_sub,
      caption = p5_cap
    ) +
    theme_lesson() +
    theme(panel.grid.major = element_line(colour = "#f2f2f2")) +
    coord_cartesian(xlim = c(0, x_right), ylim = c(0, y_top), clip = "on")
  save_lesson(p5, "05_Q_versus_chance", 11.2, 6.4)

  # -----------------------------------------------------------------
  # 6. What I² keeps
  # -----------------------------------------------------------------
  if (close) {
    p6_real <- ggplot() +
      geom_path(data = bracket_path(0, Q, 2.15), aes(x = x, y = y), linewidth = 1.15, colour = "#0072B2", lineend = "round") +
      geom_path(data = bracket_path(0, df, 1.05), aes(x = x, y = y), linewidth = 1.15, colour = col_chance, lineend = "round") +
      annotate("text", x = Q / 2, y = 2.42, label = sprintf("Observed Q = %s", fmt(Q, 3)),
               family = "DejaVu Sans", size = 3.7, colour = "#0072B2", fontface = "bold") +
      annotate("text", x = df / 2, y = 0.72, label = "Chance benchmark, df = 1",
               family = "DejaVu Sans", size = 3.7, colour = "#333333", fontface = "bold") +
      annotate("segment", x = Q, xend = df, y = 1.6, yend = 1.6, colour = col_excess,
               arrow = arrow(length = unit(2.2, "mm"), ends = "both", type = "closed")) +
      annotate("text", x = (Q + df) / 2, y = 1.82,
               label = sprintf("shortfall\nQ - df = %s", fmt(Q - df, 3)),
               family = "DejaVu Sans", size = 3.4, colour = col_excess, lineheight = 0.95) +
      annotate("text", x = 0, y = 0.22, hjust = 0,
               label = sprintf("(Q − df) / Q = %s    →    I² = max(0, %s) × 100 = 0", fmt(excess_raw, 2), fmt(excess_raw, 2)),
               family = "DejaVu Sans Mono", size = 3.3, colour = "#1a1a1a") +
      scale_x_continuous(limits = c(-0.02, 1.35), expand = c(0, 0), name = NULL, breaks = NULL) +
      scale_y_continuous(limits = c(0, 2.7), name = NULL, breaks = NULL) +
      labs(title = "PGS001336, the real Q") +
      theme_lesson(12) +
      theme(panel.grid = element_blank(), axis.text.x = element_blank())

    p6_demo <- ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0.85, ymax = 1.55, fill = "#d0d0d0", colour = NA) +
      annotate("rect", xmin = 1, xmax = 4, ymin = 0.85, ymax = 1.55, fill = "#f4c7a8", colour = NA) +
      annotate("text", x = 0.5, y = 1.2, label = "df = 1\nchance", family = "DejaVu Sans", size = 3.3, lineheight = 0.95) +
      annotate("text", x = 2.5, y = 1.2, label = "Q − df = 3\nexcess", family = "DejaVu Sans", size = 3.3, lineheight = 0.95) +
      annotate("text", x = 0, y = 1.9, hjust = 0, label = "Suppose Q = 4 and df = 1.  I² = 3 / 4 × 100 = 75",
               family = "DejaVu Sans", size = 3.5, fontface = "bold", colour = "#1a1a1a") +
      annotate("text", x = 0, y = 0.45, hjust = 0,
               label = "I² is the excess slice as a percent of the whole bar.\nThis bar is not PGS001336. It only shows the same formula.",
               family = "DejaVu Sans", size = 3.2, colour = "#333333", lineheight = 1.05) +
      scale_x_continuous(limits = c(0, 4.3), expand = c(0.02, 0), name = NULL, breaks = NULL) +
      scale_y_continuous(limits = c(0.1, 2.3), name = NULL, breaks = NULL) +
      labs(title = "Same correction, if Q had exceeded df") +
      theme_lesson(12) +
      theme(panel.grid = element_blank(), axis.text.x = element_blank())

    p6 <- p6_real / p6_demo +
      plot_layout(heights = c(1.15, 0.95)) +
      plot_annotation(
        title = "6.  I² corrects Q by removing the scatter chance is expected to produce",
        subtitle = "I² = max(0, (Q - df) / Q) × 100.\nThe numerator drops df, the average Q when sampling error is the only difference.\nWhat remains is a percent of Q. Here that remainder is negative, so I² is 0.",
        caption = "I² does not move the pooled AUC. The pooled AUC stays 0.809.\nI² = 0 means none of Q is extra heterogeneity: the studies agree at least as tightly as their intervals require.\nQ is the scatter. I² is that scatter after the chi-square benchmark has been removed.",
        theme = ann_theme()
      )
  } else {
    p6_real <- ggplot() +
      annotate("rect", xmin = 0, xmax = df, ymin = 0.9, ymax = 1.6, fill = "#bdbdbd", colour = NA) +
      annotate("rect", xmin = df, xmax = Q, ymin = 0.9, ymax = 1.6, fill = "#f4c7a8", colour = NA) +
      annotate(
        "text", x = Q * 0.52, y = 1.25,
        label = sprintf("excess = Q − df = %s\nI² = %s", fmt(Q - df, 1), fmt(I2, 1)),
        family = "DejaVu Sans", size = 3.6, lineheight = 0.95
      ) +
      annotate(
        "text", x = 0, y = 2.05, hjust = 0,
        label = sprintf("Whole bar = observed Q = %s.   Grey edge at the left = df = 1.", fmt(Q, 1)),
        family = "DejaVu Sans", size = 3.5, fontface = "bold", colour = "#1a1a1a"
      ) +
      scale_x_continuous(
        limits = c(0, Q * 1.04), breaks = c(0, 1, 50, 100, 150, round(Q)),
        name = "Q", expand = c(0.01, 0)
      ) +
      scale_y_continuous(limits = c(0.45, 2.35), name = NULL, breaks = NULL) +
      labs(title = "PGS000004, the real Q") +
      theme_lesson(12) +
      theme(panel.grid = element_blank(), axis.text.y = element_blank())

    p6_zoom <- ggplot() +
      annotate("rect", xmin = 0, xmax = df, ymin = 0.85, ymax = 1.55, fill = "#bdbdbd", colour = NA) +
      annotate("rect", xmin = df, xmax = 8, ymin = 0.85, ymax = 1.55, fill = "#f4c7a8", colour = NA) +
      annotate("text", x = 0.5, y = 1.2, label = "df = 1\nchance", family = "DejaVu Sans", size = 3.2, lineheight = 0.95) +
      annotate("text", x = 4.5, y = 1.2, label = "excess, continuing\nout to Q = 187", family = "DejaVu Sans", size = 3.3, lineheight = 0.95) +
      annotate(
        "text", x = 0, y = 0.35, hjust = 0,
        label = sprintf("(Q − df) / Q = %s    →    I² = %s", fmt(excess_raw, 3), fmt(I2, 1)),
        family = "DejaVu Sans Mono", size = 3.3, colour = "#1a1a1a"
      ) +
      scale_x_continuous(limits = c(0, 8.4), breaks = 0:8, name = "Q, first 8 units only", expand = c(0.01, 0)) +
      scale_y_continuous(limits = c(0.05, 1.9), name = NULL, breaks = NULL) +
      labs(title = "Left edge of the same bar, zoomed") +
      theme_lesson(12) +
      theme(panel.grid = element_blank(), axis.text.y = element_blank())

    p6 <- p6_real / p6_zoom +
      plot_layout(heights = c(1, 0.85)) +
      plot_annotation(
        title = "6.  I² corrects Q by removing the scatter chance is expected to produce",
        subtitle = "I² = max(0, (Q - df) / Q) × 100.\nThe numerator drops df, the average Q when sampling error is the only difference.\nWhat remains is a percent of Q. Here that percent is 99.5.",
        caption = sprintf(
          "I² does not move the pooled AUC. The pooled AUC stays %s.\nI² = %s means almost all of Q is heterogeneity beyond sampling error, so this European cell is removed by the I² filter before Figure 1C.\nThe eMERGE logit, 0.405, is the European number in the Stage-2 scale check. Q here is the within-ancestry scatter, computed after the second study is included.",
          fmt(auc_hat, 3), fmt(I2, 1)
        ),
        theme = ann_theme()
      )
  }
  save_lesson(p6, "06_what_I2_corrects", 11.2, 8.2)

  message(sprintf(
    "%s %s: Q = %.3f, df = %d, (Q - df) / Q = %.3f, I2 = %.2f, pooled AUC = %.3f",
    ex$pgs_id, ex$ancestry, Q, df, excess_raw, I2, auc_hat
  ))
  invisible(out_dir)
}

for (ex in examples) render_lesson(ex)

stale <- list.files(root_dir, pattern = "\\.(png|pdf|csv)$", full.names = TRUE)
if (length(stale)) {
  file.remove(stale)
  message("Removed ", length(stale), " files left in ", root_dir, " before the per-score folders.")
}
message("Plots and worked tables: ", normalizePath(root_dir))
