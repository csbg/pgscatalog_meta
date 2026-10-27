#!/usr/bin/env Rscript

# ============================================================
# Words used in the figure text. Column names in the tables are unchanged.
#
#   PGS            one score, one PGS identifier
#   evaluation     one AUC: one PGS, one sample set, one ancestry group
#   sample set     the people, a PSS identifier. In the model, (1|cohort)
#                  is this sample set.
#   pair           one PGS in Europeans and in one non-European ancestry
#                  group, and the difference in AUC
#   comparison     all the pairs for one trait, one training set, and one
#                  non-European ancestry group. Figure 1C has 30.
#   training set   European-only, or multi-ancestry
#
# Do not write "cell", "code", "cohorts left", "key", "bucket", or "k"
# in titles, subtitles, or captions.
#
# Eight comparisons use more than one target sample set, so the sample-set
# term can be fit. Seven are type 2 diabetes. The eighth is breast cancer
# in South Asian and European participants, European-only training.
#
# Seventeen comparisons have several PGS and one target sample set the
# model can tell apart. Ancestry group and that sample set are the same
# column. The average stays. Its standard error becomes the smallest
# standard error among the PGS. Stroke uses two sample sets, and every
# PGS uses both, so it counts with these seventeen.
#
# Five comparisons have one PGS. Nothing is shared.
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
      n_pairs == 1 ~ "One PGS",
      single_pgs_sensitivity ~ "One target sample set",
      re_applicable ~ "Several target sample sets",
      TRUE ~ "Unclassified"
    ),
    lo_min = delta_fixed - 1.96 * se_min,
    hi_min = delta_fixed + 1.96 * se_min,
    sig_min = excludes_zero(lo_min, hi_min),
    lo_med = delta_fixed - 1.96 * se_median,
    hi_med = delta_fixed + 1.96 * se_median,
    sig_med = excludes_zero(lo_med, hi_med),
    call_changes_min_se = shape == "One target sample set" & sig_fixed != sig_min,
    call_changes_one_per = shape == "Several target sample sets" & sig_fixed != sig_one_per_pss,
    label = sprintf("%s — %s, %s", trait_label, target_ancestry, train_short(trained_bucket))
  )

stopifnot(!any(cells$shape == "Unclassified"))
stopifnot(nrow(cells) == 30L)
stopifnot(sum(cells$shape == "Several target sample sets") == 8L)
stopifnot(sum(cells$shape == "One target sample set") == 17L)
stopifnot(sum(cells$shape == "One PGS") == 5L)
stopifnot(all(cells$sig_min[cells$shape == "One target sample set"] ==
  cells$sig_med[cells$shape == "One target sample set"]))

eight <- cells %>% filter(shape == "Several target sample sets")
seventeen <- cells %>% filter(shape == "One target sample set")

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
# 1. The 30 comparisons, split by how many target sample sets the PGS use
# -----------------------------------------------------------------
shape_levels <- c(
  "Several target sample sets",
  "One target sample set",
  "One PGS"
)
map_df <- cells %>%
  mutate(
    shape = factor(shape, levels = shape_levels),
    y = factor(label, levels = rev(label[order(shape, n_pairs, trait_label)]))
  )

p_map <- ggplot(map_df, aes(y = y)) +
  geom_point(aes(x = 1), size = 2.2, colour = col_share) +
  scale_x_continuous(limits = c(0.6, 1.4), breaks = NULL, name = NULL) +
  facet_grid(shape ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "1.  Which comparisons can use the sample-set term",
    subtitle = "Each row is one comparison: one trait, one training set, and Europeans compared with the ancestry group named here.\nThe top panel uses more than one target sample set. The middle panel does not. The bottom panel has one PGS.",
    y = NULL,
    caption = "Catalog stamp 20260923. In (1|cohort), cohort means the sample set.\nStroke is in the middle panel. It uses two South Asian sample sets, and every PGS uses both, so those sample sets are not separate from the ancestry group.\nThis list is not on the page for Georg."
  ) +
  theme_lesson(11) +
  theme(panel.grid.major.x = element_blank())
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
    role = if_else(kept, "Kept for this sample set", "Further PGS on that sample set")
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
    "Kept for this sample set" = 17,
    "Further PGS on that sample set" = 16
  ), guide = "none") +
  labs(
    title = "2.  One of the eight. Ten PGS share PSS011749, and two other target sample sets are in the comparison.",
    subtitle = "Type 2 diabetes, Europeans and South Asian participants, European-only training.\nTriangles are the PGS kept for each target sample set. Circles are further PGS on a sample set that already has a triangle. The black line is All PGS IVW, the inverse-variance pool of all twelve.",
    x = "ΔAUC (South Asian − European) with 95% CI",
    y = NULL,
    caption = sprintf(
      "The pool of twelve is %s (%s to %s) and includes 0. One PGS per target sample set leaves three values. Their pool is %s (%s to %s) and excludes 0.\nOn the page for Georg, lmer(AUC ~ ancestry + PGS + (1|cohort)) gives %s for South Asian participants in this training set (%s to %s). Here cohort means the sample set.\nThat number is an ancestry coefficient from the joint model. It is not a point on this axis. The interval includes 0.",
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
    kind = "PGS on PSS011736"
  ),
  tibble(
    axis = c("All PGS IVW", "Same average, smallest SE"),
    estimate = hisp_all$delta_hat,
    lo = c(hisp_all$lo, hisp_pub$lo_min),
    hi = c(hisp_all$hi, hisp_pub$hi_min),
    kind = c("All PGS IVW", "Same average, smallest SE")
  )
) %>%
  mutate(axis = factor(axis, levels = rev(axis)))

p_hisp <- ggplot(hisp_rows, aes(x = estimate, y = axis, colour = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.2, linewidth = 0.65) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c(
    "PGS on PSS011736" = col_share,
    "All PGS IVW" = col_ivw,
    "Same average, smallest SE" = col_min
  )) +
  labs(
    title = "3.  One of the seventeen. Ten PGS, and the only target sample set is PSS011736.",
    subtitle = "Type 2 diabetes, Europeans and Hispanic or Latin American participants, European-only training.\nThe rows that are Hispanic or Latin American are the rows that are this sample set. In the model, ancestry group and sample set are the same column.",
    x = "ΔAUC (Hispanic or Latin American − European) with 95% CI",
    y = NULL,
    caption = sprintf(
      "All PGS IVW, the black interval, is %s (%s to %s). It counts the same people ten times, so it is too narrow.\nThe red interval keeps that average and uses the smallest standard error, %s. It runs from %s to %s and still excludes 0.\nThis comparison does not change its call. Three of the seventeen do. They are marked on the next plot.",
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
    kind = "All PGS IVW"
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
    "All PGS IVW" = col_ivw,
    "Same average, smallest SE" = col_min
  )) +
  labs(
    title = "4.  The seventeen. The point stays. The interval returns to one PGS.",
    subtitle = "Black is All PGS IVW. Red is the same average with the smallest standard error among the PGS.\nThree red intervals include 0 where the black interval does not.",
    x = "ΔAUC (non-European − European) with 95% CI",
    y = NULL,
    caption = paste(
      "The three whose call changes are coronary heart disease, incident and prevalent, in Hispanic or Latin American participants;",
      "incident coronary heart disease in South Asian participants; and seropositive rheumatoid arthritis in South Asian participants.",
      "Stroke has two South Asian sample sets. All ten PGS use both, so the comparison has one ancestry difference.",
      "The median standard error gives the same call as the smallest in all seventeen.",
      sep = "\n"
    )
  ) +
  theme_lesson(10) +
  theme(legend.text = element_text(size = 8), axis.text.y = element_text(size = 8))
save_lesson(p_17, "04_seventeen_smallest_se", 12.6, 11.4)

# -----------------------------------------------------------------
# 5. All 8: one PGS per target sample set moves the estimate
# -----------------------------------------------------------------
eight_long <- bind_rows(
  eight %>% transmute(
    label, delta_fixed, call_changes_one_per,
    estimate = delta_fixed, lo = lo_fixed, hi = hi_fixed,
    kind = "All PGS IVW"
  ),
  eight %>% transmute(
    label, delta_fixed, call_changes_one_per,
    estimate = delta_one_per_pss, lo = lo_one_per_pss, hi = hi_one_per_pss,
    kind = "One PGS per target sample set"
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
    "All PGS IVW" = col_ivw,
    "One PGS per target sample set" = col_kept
  )) +
  labs(
    title = "5.  The eight. One PGS per target sample set is the check, because more than one sample set is in the comparison.",
    subtitle = "Black is All PGS IVW. Purple is one PGS per target sample set, the lowest PGS identifier.\nThe model Georg is looking at, with (1|cohort) equal to the sample set, is the model for this panel. It is not the red interval on the previous plot.",
    x = "ΔAUC (non-European − European) with 95% CI",
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
  "Comparisons: %d with several target sample sets, %d with one target sample set, %d with one PGS.",
  nrow(eight), nrow(seventeen), sum(cells$shape == "One PGS")
))
message("Call changes in the 8: ", paste(flip8$label, collapse = " | "))
message("Call changes in the 17: ", paste(flip17$label, collapse = " | "))
message("Illustrations: ", normalizePath(root_dir))
