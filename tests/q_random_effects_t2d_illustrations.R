#!/usr/bin/env Rscript

# ============================================================
# Two shapes of shared sample sets, for the cohort model Georg
# is being asked about.
#
# A comparison is one trait, one training bucket, one target ancestry.
# Several scores on the same people are the same dependence in both shapes.
#
#   8 comparisons leave more than one target cohort after one score
#   is kept per cohort. Nikolaus's
#     AUC ~ ancestry + PGS + (1|cohort)
#   has more than one target cohort to attach to. Seven are type 2
#   diabetes. The eighth is breast cancer in South Asian participants,
#   European-only training.
#
#   17 comparisons leave one target cohort. Ancestry and that cohort
#   are the same column, so the cohort model cannot separate them.
#   The average stays. Its standard error becomes the smallest
#   standard error among the scores.
#
#   5 comparisons have one score. Nothing is shared.
#
# DerSimonian-Laird is a different random effect and is not drawn.
#
# Run from the repository root, after analysis/09_review_sensitivity.R:
#   Rscript tests/q_random_effects_t2d_illustrations.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
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

old_names <- c(
  "01_twelve_scores_three_sample_sets",
  "02_three_summaries",
  "03_tau_squared",
  "04_one_sample_set",
  "05_four_independent_scores",
  "worked_hispanic_multi_scores",
  "worked_south_asian_pools"
)
for (nm in old_names) {
  unlink(file.path(root_dir, paste0(nm, ".png")))
  unlink(file.path(root_dir, paste0(nm, ".pdf")))
  unlink(file.path(root_dir, paste0(nm, ".csv")))
}

pair_path <- file.path(table_dir, "10_covariates_pairs.csv")
dep_path <- file.path(table_dir, "09_dependence_cells.csv")
id_path <- file.path(table_dir, "13_cell_identifiers_long.csv")
coef_path <- file.path(
  "results", paste0("review_sensitivity_", stamp),
  "diabetes_cohort_models", "diabetes_model_coefficients.csv"
)
stopifnot(file.exists(pair_path), file.exists(dep_path), file.exists(id_path), file.exists(coef_path))

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_n <- function(x) format(as.integer(round(x)), big.mark = ",")
excludes_zero <- function(lo, hi) lo > 0 | hi < 0

train_short <- function(x) {
  ifelse(x == "Train: European-only", "European-only training", "Multi-ancestry training")
}

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
      strip.text = element_text(face = "bold", hjust = 0, size = base - 1, colour = "#1a1a1a"),
      strip.background = element_rect(fill = "#F3F3F3", colour = NA),
      plot.margin = margin(12, 18, 12, 12)
    )
}

save_lesson <- function(plot, filename, width, height) {
  png_path <- file.path(root_dir, paste0(filename, ".png"))
  pdf_path <- file.path(root_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 170, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white", device = grDevices::cairo_pdf)
  message("Wrote ", png_path)
}

pool_delta <- function(delta, se) {
  w <- 1 / se^2
  sw <- sum(w)
  dh <- sum(w * delta) / sw
  tibble(
    delta_hat = dh,
    se_hat = sqrt(1 / sw),
    lo = dh - 1.96 * sqrt(1 / sw),
    hi = dh + 1.96 * sqrt(1 / sw),
    k = length(delta)
  )
}

col_share <- "#332288"
col_ivw <- "#222222"
col_min <- "#E41A1C"
col_kept <- "#6A3D9A"

# -----------------------------------------------------------------
# Classify the 30 Figure 1C comparisons
# -----------------------------------------------------------------
pairs <- read_csv(pair_path, show_col_types = FALSE)
dep <- read_csv(dep_path, show_col_types = FALSE)

se_by_cell <- pairs %>%
  group_by(trait_label, trained_bucket, target_ancestry) %>%
  summarise(
    se_min = min(delta_se),
    se_median = median(delta_se),
    .groups = "drop"
  )

cells <- dep %>%
  left_join(se_by_cell, by = c("trait_label", "trained_bucket", "target_ancestry")) %>%
  mutate(
    shape = case_when(
      n_pairs == 1 ~ "One score",
      single_pgs_sensitivity ~ "One target cohort left",
      re_applicable ~ "More than one target cohort",
      TRUE ~ "Unclassified"
    ),
    lo_min = delta_fixed - 1.96 * se_min,
    hi_min = delta_fixed + 1.96 * se_min,
    sig_min = excludes_zero(lo_min, hi_min),
    lo_med = delta_fixed - 1.96 * se_median,
    hi_med = delta_fixed + 1.96 * se_median,
    sig_med = excludes_zero(lo_med, hi_med),
    call_changes_min_se = shape == "One target cohort left" & sig_fixed != sig_min,
    call_changes_one_per = shape == "More than one target cohort" & sig_fixed != sig_one_per_pss,
    label = sprintf("%s — %s, %s", trait_label, target_ancestry, train_short(trained_bucket))
  )

stopifnot(!any(cells$shape == "Unclassified"))
stopifnot(nrow(cells) == 30L)
stopifnot(sum(cells$shape == "More than one target cohort") == 8L)
stopifnot(sum(cells$shape == "One target cohort left") == 17L)
stopifnot(sum(cells$shape == "One score") == 5L)
stopifnot(all(cells$sig_min[cells$shape == "One target cohort left"] ==
  cells$sig_med[cells$shape == "One target cohort left"]))

eight <- cells %>% filter(shape == "More than one target cohort")
seventeen <- cells %>% filter(shape == "One target cohort left")

stopifnot(sum(eight$trait_label == "Type 2 diabetes") == 7L)
stopifnot(any(eight$trait_label == "Breast cancer" & eight$target_ancestry == "South Asian"))
stopifnot(any(
  seventeen$trait_label == "Type 2 diabetes" &
    seventeen$target_ancestry == "Hispanic or Latin American" &
    seventeen$trained_bucket == "Train: European-only"
))
stopifnot(any(
  seventeen$trait_label == "Stroke excluding subarachnoid hemorrhage" &
    seventeen$tgt_n_sample_sets == 2L
))

flip17 <- seventeen %>% filter(call_changes_min_se) %>% arrange(trait_label, target_ancestry)
flip8 <- eight %>% filter(call_changes_one_per) %>% arrange(trait_label, target_ancestry)
stopifnot(nrow(flip17) == 3L, nrow(flip8) == 3L)

write_csv(
  cells %>%
    select(
      trait_label, trained_bucket, target_ancestry, shape,
      n_pairs, n_pairs_one_per_pss, tgt_n_sample_sets,
      delta_fixed, se_fixed, lo_fixed, hi_fixed, sig_fixed,
      se_min, lo_min, hi_min, sig_min,
      delta_one_per_pss, lo_one_per_pss, hi_one_per_pss, sig_one_per_pss,
      call_changes_min_se, call_changes_one_per
    ),
  file.path(root_dir, "comparison_shapes.csv")
)

# -----------------------------------------------------------------
# 1. The 30 comparisons, split by whether a cohort term has
#    more than one target cohort
# -----------------------------------------------------------------
shape_levels <- c(
  "More than one target cohort",
  "One target cohort left",
  "One score"
)
map_df <- cells %>%
  mutate(
    shape = factor(shape, levels = shape_levels),
    y = factor(label, levels = rev(label[order(shape, n_pairs, trait_label)])),
    cohorts_left = sprintf(
      "%d cohort%s left",
      n_pairs_one_per_pss,
      ifelse(n_pairs_one_per_pss == 1, "", "s")
    )
  )

p_map <- ggplot(map_df, aes(x = n_pairs, y = y)) +
  geom_segment(aes(x = 0, xend = n_pairs, yend = y), colour = "#D9D9D9", linewidth = 0.6) +
  geom_point(size = 2.8, colour = col_share) +
  geom_text(
    aes(label = cohorts_left),
    hjust = -0.15, size = 2.7, family = "DejaVu Sans", colour = "#333333"
  ) +
  scale_x_continuous(
    name = "Scores in the comparison",
    limits = c(0, 18),
    breaks = c(1, 2, 5, 10, 12)
  ) +
  facet_grid(shape ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "1.  Thirty comparisons. The cohort model can be fit in eight of them.",
    subtitle = "Each row is one trait, one training bucket, and one target ancestry.\nA cohort term needs more than one target cohort. That is the top panel. The middle panel is the same sharing with one target cohort left.",
    y = NULL,
    caption = "Catalog stamp 20260923, after one evaluation per score, sample set, and ancestry.\nStroke is in the middle panel: ten scores sit on two South Asian sample sets, and one score per set of cohorts leaves a single ancestry difference.\nThe bottom panel has one score. Nothing is shared, so neither check applies."
  ) +
  theme_lesson(11) +
  theme(panel.grid.major.x = element_line(colour = "#EEEEEE"))
save_lesson(p_map, "01_two_shapes", 12.4, 13.2)

# -----------------------------------------------------------------
# Worked cell: South Asian type 2 diabetes, one of the 8
# -----------------------------------------------------------------
ids <- read_csv(id_path, show_col_types = FALSE) %>%
  filter(side == "Target")

tgt_meta <- ids %>%
  group_by(trait_label, trained_bucket, target_ancestry, pgs_id, sampleset_id) %>%
  summarise(
    n_individuals = max(n_individuals),
    cohort_txt = paste(unique(trimws(cohort_txt)), collapse = "+"),
    .groups = "drop"
  ) %>%
  group_by(trait_label, trained_bucket, target_ancestry, pgs_id) %>%
  summarise(
    tgt_key = paste(sort(unique(sampleset_id)), collapse = ";"),
    group_label = paste(
      sprintf("%s, %s people", cohort_txt, fmt_n(n_individuals)),
      collapse = " + "
    ),
    .groups = "drop"
  )

cell_pairs <- function(trait, bucket, ancestry) {
  pairs %>%
    filter(trait_label == trait, trained_bucket == bucket, target_ancestry == ancestry) %>%
    left_join(
      tgt_meta %>% filter(trait_label == trait, trained_bucket == bucket, target_ancestry == ancestry),
      by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")
    )
}

keep_one <- function(d) {
  d %>%
    group_by(tgt_key) %>%
    arrange(pgs_id, .by_group = TRUE) %>%
    mutate(kept = row_number() == 1L, n_in_group = n()) %>%
    ungroup()
}

sa <- keep_one(cell_pairs("Type 2 diabetes", "Train: European-only", "South Asian"))
stopifnot(nrow(sa) == 12L, sum(sa$kept) == 3L, any(grepl("PSS011749", sa$tgt_key) & sa$n_in_group == 10L))
sa_all <- pool_delta(sa$delta, sa$delta_se)
sa_one <- pool_delta(sa$delta[sa$kept], sa$delta_se[sa$kept])
sa_pub <- cells %>%
  filter(trait_label == "Type 2 diabetes", trained_bucket == "Train: European-only", target_ancestry == "South Asian")
stopifnot(abs(sa_all$delta_hat - sa_pub$delta_fixed) < 1e-8)
stopifnot(abs(sa_one$delta_hat - sa_pub$delta_one_per_pss) < 1e-8)
stopifnot(isTRUE(sa_pub$sig_one_per_pss), !isTRUE(sa_pub$sig_fixed))

sa <- sa %>%
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

georg_sa <- read_csv(coef_path, show_col_types = FALSE) %>%
  filter(
    train_bucket == "Train: European-only",
    ancestry == "South Asian",
    cohort_def == "sample-set id",
    model == "lmer(AUC ~ ancestry + PGS + (1|cohort), weights = 1/SE^2)"
  )
stopifnot(nrow(georg_sa) == 1L)

sa_forest <- sa %>%
  mutate(
    axis = factor(pgs_id, levels = rev(pgs_id[order(group_label, pgs_id)])),
    role = if_else(kept, "Kept for this cohort", "Further score on that cohort")
  )

p_sa <- ggplot(sa_forest, aes(x = delta, y = axis, colour = group_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(xmin = delta - 1.96 * delta_se, xmax = delta + 1.96 * delta_se),
    orientation = "y", width = 0.2, linewidth = 0.6
  ) +
  geom_point(aes(shape = role), size = 2.8) +
  geom_vline(xintercept = sa_all$delta_hat, colour = col_ivw, linewidth = 0.45) +
  scale_shape_manual(values = c(
    "Kept for this cohort" = 17,
    "Further score on that cohort" = 16
  ), guide = "none") +
  labs(
    title = "2.  One of the eight. Ten scores share PSS011749, and other target cohorts remain.",
    subtitle = "Type 2 diabetes, South Asian target, European-only training.\nTriangles are the score kept for each target cohort. Circles are further scores on a cohort that already has a triangle. The black line is the inverse-variance pool of all twelve.",
    x = "ΔAUC (South Asian − European) with 95% CI",
    y = NULL,
    caption = sprintf(
      "The pool of twelve is %s (%s to %s) and includes 0. One score per target cohort leaves three values. Their pool is %s (%s to %s) and excludes 0.\nOn the page for Georg, lmer(AUC ~ ancestry + PGS + (1|cohort)) gives %s for South Asian in this training bucket (%s to %s).\nThat number is an ancestry coefficient from the joint model. It is not a point on this ΔAUC axis. The interval includes 0.",
      fmt(sa_all$delta_hat, 4), fmt(sa_all$lo, 4), fmt(sa_all$hi, 4),
      fmt(sa_one$delta_hat, 3), fmt(sa_one$lo, 3), fmt(sa_one$hi, 3),
      fmt(georg_sa$estimate, 3), fmt(georg_sa$lo, 3), fmt(georg_sa$hi, 3)
    )
  ) +
  guides(colour = guide_legend(nrow = 3)) +
  theme_lesson(11) +
  theme(
    legend.box = "vertical",
    legend.text = element_text(size = 8)
  )
save_lesson(p_sa, "02_south_asian_one_of_eight", 11.8, 8.8)

# -----------------------------------------------------------------
# Worked cell: Hispanic type 2 diabetes, one of the 17
# -----------------------------------------------------------------
hisp <- keep_one(cell_pairs(
  "Type 2 diabetes", "Train: European-only", "Hispanic or Latin American"
))
stopifnot(nrow(hisp) == 10L, n_distinct(hisp$tgt_key) == 1L, all(hisp$tgt_key == "PSS011736"))
hisp_all <- pool_delta(hisp$delta, hisp$delta_se)
hisp_pub <- cells %>%
  filter(
    trait_label == "Type 2 diabetes",
    trained_bucket == "Train: European-only",
    target_ancestry == "Hispanic or Latin American"
  )
stopifnot(abs(hisp_all$delta_hat - hisp_pub$delta_fixed) < 1e-8)
stopifnot(isTRUE(hisp_pub$sig_fixed), isTRUE(hisp_pub$sig_min), !isTRUE(hisp_pub$call_changes_min_se))

hisp_rows <- bind_rows(
  hisp %>% transmute(
    axis = pgs_id,
    estimate = delta,
    lo = delta - 1.96 * delta_se,
    hi = delta + 1.96 * delta_se,
    kind = "Score on PSS011736"
  ),
  tibble(
    axis = c("Inverse-variance of all 10", "Same average, smallest SE"),
    estimate = hisp_all$delta_hat,
    lo = c(hisp_all$lo, hisp_pub$lo_min),
    hi = c(hisp_all$hi, hisp_pub$hi_min),
    kind = c("Inverse-variance of all 10", "Same average, smallest SE")
  )
) %>%
  mutate(axis = factor(axis, levels = rev(axis)))

p_hisp <- ggplot(hisp_rows, aes(x = estimate, y = axis, colour = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.2, linewidth = 0.65) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c(
    "Score on PSS011736" = col_share,
    "Inverse-variance of all 10" = col_ivw,
    "Same average, smallest SE" = col_min
  )) +
  labs(
    title = "3.  One of the seventeen. Ten scores, and the only target cohort is PSS011736.",
    subtitle = "Type 2 diabetes, Hispanic or Latin American target, European-only training.\nThe rows that are Hispanic are the rows that are this sample set. A cohort term and the ancestry term are the same column.",
    x = "ΔAUC (Hispanic or Latin American − European) with 95% CI",
    y = NULL,
    caption = sprintf(
      "The black interval is the inverse-variance pool, %s (%s to %s). It counts the same people ten times, so it is too narrow.\nThe red interval keeps that average and uses the smallest standard error, %s. It runs from %s to %s and still excludes 0.\nThis comparison does not change its call. Three of the seventeen do. They are marked on the next plot.",
      fmt(hisp_all$delta_hat, 3), fmt(hisp_all$lo, 3), fmt(hisp_all$hi, 3),
      fmt(hisp_pub$se_min, 4), fmt(hisp_pub$lo_min, 3), fmt(hisp_pub$hi_min, 3)
    )
  ) +
  theme_lesson(11) +
  theme(legend.text = element_text(size = 8))
save_lesson(p_hisp, "03_hispanic_one_of_seventeen", 11.4, 7.6)

# -----------------------------------------------------------------
# 4. All 17: keep the average, replace the standard error
# -----------------------------------------------------------------
sev_long <- bind_rows(
  seventeen %>% transmute(
    label, delta_fixed, call_changes_min_se,
    estimate = delta_fixed, lo = lo_fixed, hi = hi_fixed,
    kind = "Inverse-variance interval"
  ),
  seventeen %>% transmute(
    label, delta_fixed, call_changes_min_se,
    estimate = delta_fixed, lo = lo_min, hi = hi_min,
    kind = "Same average, smallest SE"
  )
) %>%
  mutate(
    y_label = if_else(
      call_changes_min_se,
      paste0(label, "\ncall changes"),
      label
    ),
    y_label = factor(y_label, levels = rev(unique(y_label[order(delta_fixed, label)])))
  )

pd <- position_dodge(width = 0.62)
p_17 <- ggplot(sev_long, aes(x = estimate, y = y_label, colour = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(xmin = lo, xmax = hi),
    position = pd, orientation = "y", width = 0.28, linewidth = 0.65
  ) +
  geom_point(position = pd, size = 2.3) +
  scale_colour_manual(values = c(
    "Inverse-variance interval" = col_ivw,
    "Same average, smallest SE" = col_min
  )) +
  labs(
    title = "4.  The seventeen. The point stays. The interval returns to one score.",
    subtitle = "Black is the inverse-variance interval. Red is the same average with the smallest standard error among the scores.\nThree red intervals include 0 where the black interval does not.",
    x = "ΔAUC (target − European) with 95% CI",
    y = NULL,
    caption = paste(
      "The three whose call changes are coronary heart disease, incident and prevalent, in Hispanic or Latin American participants;",
      "incident coronary heart disease in South Asian participants; and seropositive rheumatoid arthritis in South Asian participants.",
      "Stroke has two South Asian sample sets. All ten scores use both, so one key remains and the cohort model still has a single ancestry difference.",
      "The median standard error gives the same call as the smallest in all seventeen.",
      sep = "\n"
    )
  ) +
  theme_lesson(10) +
  theme(legend.text = element_text(size = 8), axis.text.y = element_text(size = 8))
save_lesson(p_17, "04_seventeen_smallest_se", 12.6, 11.4)

# -----------------------------------------------------------------
# 5. All 8: one score per target cohort moves the estimate
# -----------------------------------------------------------------
eight_long <- bind_rows(
  eight %>% transmute(
    label, delta_fixed, call_changes_one_per,
    estimate = delta_fixed, lo = lo_fixed, hi = hi_fixed,
    kind = "All scores"
  ),
  eight %>% transmute(
    label, delta_fixed, call_changes_one_per,
    estimate = delta_one_per_pss, lo = lo_one_per_pss, hi = hi_one_per_pss,
    kind = "One score per target cohort"
  )
) %>%
  mutate(
    y_label = if_else(
      call_changes_one_per,
      paste0(label, "\ncall changes"),
      label
    ),
    y_label = factor(y_label, levels = rev(unique(y_label[order(delta_fixed, label)])))
  )

p_8 <- ggplot(eight_long, aes(x = estimate, y = y_label, colour = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(xmin = lo, xmax = hi),
    position = pd, orientation = "y", width = 0.28, linewidth = 0.7
  ) +
  geom_point(position = pd, size = 2.5) +
  scale_colour_manual(values = c(
    "All scores" = col_ivw,
    "One score per target cohort" = col_kept
  )) +
  labs(
    title = "5.  The eight. One score per target cohort is the check, because more than one cohort remains.",
    subtitle = "Black is every score. Purple is one score per target cohort, the lowest PGS identifier.\nThe cohort model Georg is looking at is the model for this panel. It is not the correction drawn in red on the previous plot.",
    x = "ΔAUC (target − European) with 95% CI",
    y = NULL,
    caption = paste(
      "The call changes in three: both European-only type 2 diabetes comparisons, and breast cancer in South Asian participants.",
      "South Asian, European-only, moves from about 0 to a deficit. East Asian, European-only, and breast cancer move from excluding 0 to including it.",
      "The other five, including Hispanic or Latin American type 2 diabetes with multi-ancestry training, keep the same call.",
      sep = "\n"
    )
  ) +
  theme_lesson(11) +
  theme(legend.text = element_text(size = 8), axis.text.y = element_text(size = 8))
save_lesson(p_8, "05_eight_one_score_per_cohort", 12.4, 8.2)

message(sprintf(
  "Shapes: %d with more than one target cohort, %d with one target cohort left, %d with one score.",
  nrow(eight), nrow(seventeen), sum(cells$shape == "One score")
))
message("Call changes in the 8: ", paste(flip8$label, collapse = " | "))
message("Call changes in the 17: ", paste(flip17$label, collapse = " | "))
message("Illustrations: ", normalizePath(root_dir))
