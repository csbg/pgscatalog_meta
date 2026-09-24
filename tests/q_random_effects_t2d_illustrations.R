#!/usr/bin/env Rscript

# ============================================================
# Didactic walkthrough of the type 2 diabetes random-effects ΔAUC
# DerSimonian-Laird is fit only after keeping one PGS per target
# sample set, and only when at least two ΔAUC values remain.
#
# Worked cells:
#   European-only training, South Asian target
#     12 scores, 3 target sample-set groups.
#     The pool of 12 includes 0. The pool of 3, and the random-effects
#     pool of those 3, are both negative and exclude 0. I² is about 51.
#   European-only training, Hispanic or Latin American target
#     10 scores, 1 target sample set. Random effects are not fit.
#   Multi-ancestry training, Hispanic or Latin American target
#     5 scores, 4 target sample sets. I² is about 97.
#     The random-effects interval is wide and still excludes 0.
#
# Run from the repository root, after analysis/09_review_sensitivity.R:
#   Rscript tests/q_random_effects_t2d_illustrations.R
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
table_dir <- file.path("results", paste0("review_sensitivity_", stamp), "tables")
root_dir <- file.path("results", paste0("review_sensitivity_", stamp), "q_random_effects_t2d_illustrations")
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

pair_path <- file.path(table_dir, "10_covariates_pairs.csv")
dep_path <- file.path(table_dir, "09_dependence_cells.csv")
id_path <- file.path(table_dir, "13_cell_identifiers_long.csv")
stopifnot(file.exists(pair_path), file.exists(dep_path), file.exists(id_path))

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
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
      legend.position = "bottom",
      legend.title = element_blank(),
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

save_lesson <- function(plot, filename, width, height) {
  png_path <- file.path(root_dir, paste0(filename, ".png"))
  pdf_path <- file.path(root_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 170, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white", device = grDevices::cairo_pdf)
  message("Wrote ", png_path)
}

pool_delta <- function(delta, se, model = c("fixed", "random_dl")) {
  model <- match.arg(model)
  w <- 1 / se^2
  sw <- sum(w)
  dh <- sum(w * delta) / sw
  Q <- sum(w * (delta - dh)^2)
  k <- length(delta)
  df <- k - 1
  tau2 <- 0
  if (model == "random_dl" && k >= 2) {
    C <- sw - sum(w^2) / sw
    tau2 <- if (is.finite(C) && C > 0) max(0, (Q - df) / C) else 0
  }
  ws <- 1 / (se^2 + tau2)
  dh2 <- sum(ws * delta) / sum(ws)
  se_hat <- sqrt(1 / sum(ws))
  I2 <- if (k >= 2 && Q > 0) max(0, (Q - df) / Q) * 100 else NA_real_
  tibble(
    delta_hat = dh2, se_hat = se_hat,
    lo = dh2 - 1.96 * se_hat, hi = dh2 + 1.96 * se_hat,
    tau2 = tau2, Q = Q, df = df, I2 = I2, k = k,
    C = if (k >= 2) sw - sum(w^2) / sw else NA_real_
  )
}

group_colours <- c("#332288", "#E41A1C", "#01665E", "#6A3D9A")
col_ivw12 <- "#222222"
col_ivw_indep <- "#6A3D9A"
col_re <- "#01665E"

pairs <- read_csv(pair_path, show_col_types = FALSE) %>%
  filter(trait_label == "Type 2 diabetes")
dep <- read_csv(dep_path, show_col_types = FALSE) %>%
  filter(trait_label == "Type 2 diabetes")
ids <- read_csv(id_path, show_col_types = FALSE) %>%
  filter(trait_label == "Type 2 diabetes", side == "Target")

tgt_meta <- ids %>%
  group_by(trained_bucket, target_ancestry, pgs_id, sampleset_id) %>%
  summarise(
    n_individuals = max(n_individuals),
    cohort_txt = paste(unique(trimws(cohort_txt)), collapse = "+"),
    .groups = "drop"
  ) %>%
  group_by(trained_bucket, target_ancestry, pgs_id) %>%
  summarise(
    tgt_key = paste(sort(unique(sampleset_id)), collapse = ";"),
    group_label = paste(
      sprintf("%s, %s people", cohort_txt, fmt_n(n_individuals)),
      collapse = " + "
    ),
    .groups = "drop"
  )

cell_pairs <- function(bucket, ancestry) {
  pairs %>%
    filter(trained_bucket == bucket, target_ancestry == ancestry) %>%
    left_join(
      tgt_meta %>% filter(trained_bucket == bucket, target_ancestry == ancestry),
      by = c("trained_bucket", "target_ancestry", "pgs_id")
    )
}

keep_one <- function(d) {
  d %>%
    group_by(tgt_key) %>%
    arrange(pgs_id, .by_group = TRUE) %>%
    mutate(kept = row_number() == 1L, n_in_group = n()) %>%
    ungroup()
}

expect_cell <- function(bucket, ancestry, got_all, got_one, got_re) {
  published <- dep %>% filter(trained_bucket == bucket, target_ancestry == ancestry)
  stopifnot(nrow(published) == 1L)
  stopifnot(abs(got_all$delta_hat - published$delta_fixed) < 1e-8)
  stopifnot(got_one$k == published$n_pairs_one_per_pss)
  stopifnot(abs(got_one$delta_hat - published$delta_one_per_pss) < 1e-8)
  if (published$re_applicable) {
    stopifnot(abs(got_re$delta_hat - published$delta_random) < 1e-8)
    stopifnot(abs(got_re$I2 - published$I2_stage2) < 1e-6)
    stopifnot(abs(got_re$tau2 - published$tau2_random) < 1e-12)
  } else {
    stopifnot(got_one$k == 1L)
    stopifnot(is.na(published$delta_random))
  }
  published
}

# -----------------------------------------------------------------
# European-only, South Asian: 12 scores, 3 independent ΔAUC
# -----------------------------------------------------------------
sa <- keep_one(cell_pairs("Train: European-only", "South Asian"))
stopifnot(nrow(sa) == 12L, sum(sa$kept) == 3L)
sa_all <- pool_delta(sa$delta, sa$delta_se, "fixed")
sa_one <- pool_delta(sa$delta[sa$kept], sa$delta_se[sa$kept], "fixed")
sa_re <- pool_delta(sa$delta[sa$kept], sa$delta_se[sa$kept], "random_dl")
sa_pub <- expect_cell("Train: European-only", "South Asian", sa_all, sa_one, sa_re)
stopifnot(isTRUE(sa_pub$sig_one_per_pss), isTRUE(sa_pub$sig_random), !isTRUE(sa_pub$sig_fixed))

sa <- keep_one(cell_pairs("Train: European-only", "South Asian")) %>%
  group_by(tgt_key) %>%
  mutate(group_label = if_else(
    n() > 1,
    sprintf("%s  (%d scores on this sample set)", group_label[1], n()),
    group_label[1]
  )) %>%
  ungroup()

write_csv(
  sa %>% select(pgs_id, delta, delta_se, tgt_key, group_label, kept, n_in_group),
  file.path(root_dir, "worked_south_asian_scores.csv")
)
write_csv(
  bind_rows(
    sa_all %>% mutate(estimate = "IVW, all 12 scores"),
    sa_one %>% mutate(estimate = "IVW, 3 independent scores"),
    sa_re %>% mutate(estimate = "Random effects, 3 independent scores")
  ),
  file.path(root_dir, "worked_south_asian_pools.csv")
)

sa_forest <- sa %>%
  mutate(
    axis = factor(pgs_id, levels = rev(pgs_id[order(group_label, pgs_id)])),
    role = if_else(kept, "Kept for this sample set", "Extra score on that sample set")
  )

p1 <- ggplot(sa_forest, aes(x = delta, y = axis, colour = group_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(xmin = delta - 1.96 * delta_se, xmax = delta + 1.96 * delta_se),
    orientation = "y", width = 0.2, linewidth = 0.6
  ) +
  geom_point(aes(shape = role), size = 2.8) +
  geom_vline(xintercept = sa_all$delta_hat, linetype = "solid", colour = col_ivw12, linewidth = 0.4) +
  scale_shape_manual(values = c(
    "Kept for this sample set" = 17,
    "Extra score on that sample set" = 16
  ), guide = "none") +
  labs(
    title = "1.  Twelve South Asian ΔAUC values, and three target sample sets",
    subtitle = "Type 2 diabetes, European-only training. Triangles are the score kept for each sample set (lowest PGS ID).\nCircles are further scores on a sample set that already has a triangle. The black line is the inverse-variance pool of all twelve.",
    x = "ΔAUC (South Asian − European) with 95% CI",
    y = NULL,
    caption = "Ten scores share one UK Biobank South Asian sample set. Pooling all twelve counts that sample ten times.\nThe pool of twelve is about 0 and its interval includes 0. Plots 2 and 3 use the three triangles only."
  ) +
  guides(colour = guide_legend(nrow = 3)) +
  theme_lesson(11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 8)
  )
save_lesson(p1, "01_twelve_scores_three_sample_sets", 11.6, 8.6)

# -----------------------------------------------------------------
# 2. Fixed effect of 12, fixed effect of 3, random effects of 3
# -----------------------------------------------------------------
three <- sa %>%
  filter(kept) %>%
  arrange(delta)

summary_rows <- bind_rows(
  three %>% transmute(
    axis = pgs_id,
    estimate = delta, lo = delta - 1.96 * delta_se, hi = delta + 1.96 * delta_se,
    kind = "score"
  ),
  tibble(
    axis = c("IVW of all 12", "IVW of these 3", "Random effects of these 3"),
    estimate = c(sa_all$delta_hat, sa_one$delta_hat, sa_re$delta_hat),
    lo = c(sa_all$lo, sa_one$lo, sa_re$lo),
    hi = c(sa_all$hi, sa_one$hi, sa_re$hi),
    kind = c("ivw12", "ivw3", "re")
  )
) %>%
  mutate(axis = factor(axis, levels = rev(axis)))

p2 <- ggplot(summary_rows, aes(x = estimate, y = axis, colour = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.18, linewidth = 0.75) +
  geom_point(aes(shape = kind), size = 3.2) +
  scale_colour_manual(values = c(
    score = "#332288", ivw12 = col_ivw12, ivw3 = col_ivw_indep, re = col_re
  ), guide = "none") +
  scale_shape_manual(values = c(score = 17, ivw12 = 18, ivw3 = 18, re = 15), guide = "none") +
  labs(
    title = "2.  The twelve-score pool and the three-score pools answer different questions",
    subtitle = "Purple diamond: inverse-variance fixed effect of the three independent scores.\nGreen square: DerSimonian-Laird random effects of those same three. Both exclude 0. The black diamond, all twelve, does not.",
    x = "ΔAUC (South Asian − European) with 95% CI",
    y = NULL,
    caption = sprintf(
      "IVW of 12: %s (%s to %s).  IVW of 3: %s (%s to %s).  Random effects: %s (%s to %s).\nPGS004106, evaluated in Genes & Health and in UK Biobank, carries most of the fixed-effect weight among the three.\nRandom effects gives the other two scores a little more say. The interval widens and still excludes 0.",
      fmt(sa_all$delta_hat, 3), fmt(sa_all$lo, 3), fmt(sa_all$hi, 3),
      fmt(sa_one$delta_hat, 3), fmt(sa_one$lo, 3), fmt(sa_one$hi, 3),
      fmt(sa_re$delta_hat, 3), fmt(sa_re$lo, 3), fmt(sa_re$hi, 3)
    )
  ) +
  theme_lesson(12)
save_lesson(p2, "02_three_summaries", 11.2, 6.4)

# -----------------------------------------------------------------
# 3. The τ² arithmetic
# -----------------------------------------------------------------
kept3 <- three
w <- 1 / kept3$delta_se^2
ws <- 1 / (kept3$delta_se^2 + sa_re$tau2)
weight_df <- bind_rows(
  tibble(pgs_id = kept3$pgs_id, share = w / sum(w), model = "Fixed effect\nweight = 1 / SE²"),
  tibble(pgs_id = kept3$pgs_id, share = ws / sum(ws), model = "Random effects\nweight = 1 / (SE² + τ²)")
)

p3_w <- ggplot(weight_df, aes(x = model, y = share, fill = pgs_id)) +
  geom_col(width = 0.62, colour = "white") +
  geom_text(
    data = function(d) dplyr::filter(d, share >= 0.08),
    aes(label = sprintf("%s\n%s%%", pgs_id, fmt(100 * share, 0))),
    position = position_stack(vjust = 0.5),
    family = "DejaVu Sans", size = 3.1, colour = "white", lineheight = 0.9
  ) +
  scale_y_continuous(name = "Share of the weight", labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0)) +
  scale_x_discrete(name = NULL) +
  scale_fill_manual(values = setNames(c("#332288", "#E41A1C", "#01665E"), kept3$pgs_id), guide = "none") +
  labs(title = "Who holds the weight") +
  theme_lesson(12)

arith <- c(
  "Three independent ΔAUC values",
  "",
  sprintf("Q  =  Σ  w (Δ − Δ_IVW)²   =   %s", fmt(sa_re$Q, 2)),
  "df  =  k − 1  =  2",
  sprintf("τ²  =  max(0, (Q − df) / C)   =   %s", fmt(sa_re$tau2, 5)),
  sprintf("I²  =  max(0, (Q − df) / Q) × 100   =   %s", fmt(sa_re$I2, 1)),
  "",
  "C is the usual DerSimonian-Laird denominator,",
  "sum(w) − sum(w²) / sum(w).",
  "",
  "τ² is added to every SE² before the weights",
  "are rebuilt. The pooled ΔAUC moves only",
  "because those weights move.",
  "",
  sprintf("Pooled random-effects ΔAUC = %s", fmt(sa_re$delta_hat, 3)),
  sprintf("95%% interval %s to %s", fmt(sa_re$lo, 3), fmt(sa_re$hi, 3))
)
p3_txt <- ggplot() +
  annotate(
    "text", x = 0, y = rev(seq_along(arith)), label = arith,
    hjust = 0, family = "DejaVu Sans Mono", size = 3.15, colour = "#1a1a1a"
  ) +
  xlim(-0.02, 1.05) +
  ylim(0.2, length(arith) + 0.8) +
  theme_void() +
  theme(plot.margin = margin(18, 8, 8, 8))

p3 <- p3_w + p3_txt +
  plot_layout(widths = c(0.85, 1.15)) +
  plot_annotation(
    title = "3.  Random effects on these three scores adds a small τ² and keeps the deficit",
    subtitle = "I² near 51 means about half of Q is scatter beyond the sampling error expected for three estimates.\nτ² is small next to the variance of the most precise score, so the weights shift and the conclusion stays a deficit.",
    caption = "Q and I² are computed from the fixed-effect weights, then τ² rebuilds the weights. The same arithmetic, on the same three ΔAUC values, is what the review figure draws in green for this cell.",
    theme = ann_theme()
  )
save_lesson(p3, "03_tau_squared", 12.2, 6.6)

# -----------------------------------------------------------------
# 4. One sample set: random effects are not fit
# -----------------------------------------------------------------
hisp_eur <- keep_one(cell_pairs("Train: European-only", "Hispanic or Latin American"))
hisp_eur_all <- pool_delta(hisp_eur$delta, hisp_eur$delta_se, "fixed")
hisp_eur_one <- pool_delta(hisp_eur$delta[hisp_eur$kept], hisp_eur$delta_se[hisp_eur$kept], "fixed")
hisp_eur_pub <- expect_cell(
  "Train: European-only", "Hispanic or Latin American",
  hisp_eur_all, hisp_eur_one, hisp_eur_one
)
stopifnot(n_distinct(hisp_eur$tgt_key) == 1L, sum(hisp_eur$kept) == 1L)

kept_id <- hisp_eur$pgs_id[hisp_eur$kept]
hisp_plot <- hisp_eur %>%
  mutate(
    axis = factor(pgs_id, levels = rev(sort(pgs_id))),
    role = if_else(kept, "Kept", "Same sample set")
  )

p4 <- ggplot(hisp_plot, aes(x = delta, y = axis, colour = role)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(xmin = delta - 1.96 * delta_se, xmax = delta + 1.96 * delta_se),
    orientation = "y", width = 0.2, linewidth = 0.6
  ) +
  geom_point(size = 2.6) +
  geom_vline(xintercept = hisp_eur_all$delta_hat, colour = col_ivw12, linewidth = 0.4) +
  scale_colour_manual(values = c(Kept = "#E41A1C", "Same sample set" = "#332288")) +
  labs(
    title = "4.  Ten Hispanic or Latin American scores, one target sample set: no random effect",
    subtitle = sprintf(
      "European-only training. Every score uses %s.\n%s is the kept score. The black line is the inverse-variance pool of all ten, which is the only pooled estimate drawn.",
      unique(hisp_eur$group_label), kept_id
    ),
    x = "ΔAUC (Hispanic or Latin American − European) with 95% CI",
    y = NULL,
    caption = "After one score per target sample set, a single ΔAUC remains. τ² is not estimated from one number.\nThe review figure draws the fixed effect of all ten and the single retained score. It does not draw a random-effects point."
  ) +
  theme_lesson(11)
save_lesson(p4, "04_one_sample_set", 11.0, 7.4)

# -----------------------------------------------------------------
# 5. Multi-ancestry Hispanic: four independent scores, large I²
# -----------------------------------------------------------------
hisp_multi <- keep_one(cell_pairs("Train: Multi incl. EUR", "Hispanic or Latin American"))
stopifnot(nrow(hisp_multi) == 5L, sum(hisp_multi$kept) == 4L)
hm_all <- pool_delta(hisp_multi$delta, hisp_multi$delta_se, "fixed")
hm_one <- pool_delta(hisp_multi$delta[hisp_multi$kept], hisp_multi$delta_se[hisp_multi$kept], "fixed")
hm_re <- pool_delta(hisp_multi$delta[hisp_multi$kept], hisp_multi$delta_se[hisp_multi$kept], "random_dl")
hm_pub <- expect_cell("Train: Multi incl. EUR", "Hispanic or Latin American", hm_all, hm_one, hm_re)
stopifnot(isTRUE(hm_pub$sig_fixed), isTRUE(hm_pub$sig_one_per_pss), isTRUE(hm_pub$sig_random))

hm_rows <- bind_rows(
  hisp_multi %>% transmute(
    axis = pgs_id,
    estimate = delta,
    lo = delta - 1.96 * delta_se,
    hi = delta + 1.96 * delta_se,
    kind = if_else(kept, "Independent score", "Second score on a sample set already kept")
  ),
  tibble(
    axis = c("IVW of all 5", "IVW of the 4", "Random effects of the 4"),
    estimate = c(hm_all$delta_hat, hm_one$delta_hat, hm_re$delta_hat),
    lo = c(hm_all$lo, hm_one$lo, hm_re$lo),
    hi = c(hm_all$hi, hm_one$hi, hm_re$hi),
    kind = c("IVW of all 5", "IVW of the 4", "Random effects of the 4")
  )
) %>%
  mutate(axis = factor(axis, levels = rev(axis)), kind = factor(kind, levels = unique(kind)))

write_csv(
  hisp_multi %>% select(pgs_id, delta, delta_se, tgt_key, group_label, kept),
  file.path(root_dir, "worked_hispanic_multi_scores.csv")
)

p5 <- ggplot(hm_rows, aes(x = estimate, y = axis, colour = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.18, linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = c(
    "Independent score" = "#332288",
    "Second score on a sample set already kept" = "#B0B0B0",
    "IVW of all 5" = col_ivw12,
    "IVW of the 4" = col_ivw_indep,
    "Random effects of the 4" = col_re
  )) +
  labs(
    title = "5.  Four independent Hispanic or Latin American ΔAUC values, and a wide random-effects interval",
    subtitle = sprintf(
      "Multi-ancestry training. Two scores share one All of Us sample set; the higher PGS ID is the grey point and leaves the random-effects pool.\nI² = %s. The fixed-effect interval of all five is narrow. The random-effects interval runs from %s to %s and still excludes 0.",
      fmt(hm_re$I2, 0), fmt(hm_re$lo, 3), fmt(hm_re$hi, 3)
    ),
    x = "ΔAUC (Hispanic or Latin American − European) with 95% CI",
    y = NULL,
    caption = "One of the four independent scores is negative and the other three are positive, which is why Q is large.\nRandom effects is the interval that treats that disagreement as extra uncertainty. Here the extra uncertainty is not enough to include 0.\nThis cell, and the South Asian cell in plots 1–3, are the type 2 diabetes random-effects results with three or more independent ΔAUC values."
  ) +
  theme_lesson(11) +
  theme(legend.text = element_text(size = 8))
save_lesson(p5, "05_four_independent_scores", 11.4, 7.2)

message(sprintf(
  "South Asian: IVW12 = %.4f, IVW3 = %.4f, RE = %.4f, I2 = %.1f",
  sa_all$delta_hat, sa_one$delta_hat, sa_re$delta_hat, sa_re$I2
))
message(sprintf(
  "Hispanic multi: IVW5 = %.4f, IVW4 = %.4f, RE = %.4f, I2 = %.1f",
  hm_all$delta_hat, hm_one$delta_hat, hm_re$delta_hat, hm_re$I2
))
message("Random-effects illustrations: ", normalizePath(root_dir))
