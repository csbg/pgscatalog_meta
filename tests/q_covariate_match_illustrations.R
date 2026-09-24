#!/usr/bin/env Rscript

# ============================================================
# Didactic walkthrough of the covariate-matched ΔAUC
# The review figure draws a second point only when matching
# removes at least one pair and leaves at least one.
#
# Worked cells, all from the 20260923 review tables:
#   Coronary heart disease, multi-ancestry training
#     African, and Hispanic or Latin American
#     k = 4 pairs, 2 matched. The point moves. Inference does not.
#   Ovarian cancer, European-only, East Asian
#     k = 2, both sides "None (PRS only)". The matched pool is the same number.
#   Breast carcinoma, European-only, African
#     k = 4, every pair mixes two classes on a side. No matched pair.
#
# The two retained coronary heart disease scores still share one
# target sample set, so this script does not estimate a random effect.
#
# Run from the repository root, after analysis/09_review_sensitivity.R:
#   Rscript tests/q_covariate_match_illustrations.R
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
root_dir <- file.path("results", paste0("review_sensitivity_", stamp), "q_covariate_match_illustrations")
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

pair_path <- file.path(table_dir, "10_covariates_pairs.csv")
cell_path <- file.path(table_dir, "10_covariates_cells_matched_sensitivity.csv")
id_path <- file.path(table_dir, "13_cell_identifiers_long.csv")
stopifnot(file.exists(pair_path), file.exists(cell_path), file.exists(id_path))

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

pool_delta <- function(delta, se) {
  w <- 1 / se^2
  dh <- sum(w * delta) / sum(w)
  se_hat <- sqrt(1 / sum(w))
  tibble(
    delta_hat = dh,
    se_hat = se_hat,
    lo = dh - 1.96 * se_hat,
    hi = dh + 1.96 * se_hat,
    k = length(delta)
  )
}

col_kept <- "#332288"
col_drop <- "#B0B0B0"
col_pool_all <- "#222222"
col_pool_matched <- "#E41A1C"

pairs <- read_csv(pair_path, show_col_types = FALSE)
cells <- read_csv(cell_path, show_col_types = FALSE)
ids <- read_csv(id_path, show_col_types = FALSE)

chd_trait <- "Coronary heart disease (incident and prevalent)"
chd_bucket <- "Train: Multi incl. EUR"

chd <- pairs %>%
  filter(trait_label == chd_trait, trained_bucket == chd_bucket) %>%
  mutate(
    decision = if_else(covariate_matched, "Kept", "Dropped"),
    eur_short = if_else(
      covariate_matched,
      "Clinical or other",
      "Age/sex/PCs  |  Clinical"
    ),
    tgt_short = "Clinical or other"
  )

stopifnot(nrow(chd) == 8L)
stopifnot(sum(chd$covariate_matched) == 4L)

tgt_n <- ids %>%
  filter(trait_label == chd_trait, trained_bucket == chd_bucket, side == "Target") %>%
  distinct(target_ancestry, pgs_id, sampleset_id, n_individuals, cohort_txt)

chd <- chd %>%
  left_join(tgt_n, by = c("target_ancestry", "pgs_id"))

# One target sample set per ancestry, shared by every score in the cell.
shared_pss <- chd %>%
  distinct(target_ancestry, sampleset_id, n_individuals) %>%
  group_by(target_ancestry) %>%
  summarise(
    n_sets = n_distinct(sampleset_id),
    sampleset_id = paste(unique(sampleset_id), collapse = ","),
    n_individuals = unique(n_individuals),
    .groups = "drop"
  )
stopifnot(all(shared_pss$n_sets == 1L))

worked_pools <- chd %>%
  group_by(target_ancestry) %>%
  group_modify(\(d, key) {
    all_p <- pool_delta(d$delta, d$delta_se)
    kept <- d %>% filter(covariate_matched)
    mat_p <- pool_delta(kept$delta, kept$delta_se)
    bind_rows(
      all_p %>% mutate(pool = "All 4 pairs"),
      mat_p %>% mutate(pool = "2 matched pairs")
    )
  }) %>%
  ungroup()

published <- cells %>%
  filter(trait_label == chd_trait, trained_bucket == chd_bucket) %>%
  select(target_ancestry, n_pairs, n_pairs_covariate_matched, delta_all, delta_matched, sig_all, sig_matched)
stopifnot(all(published$n_pairs == 4L), all(published$n_pairs_covariate_matched == 2L))

check_pool <- function(ancestry, pool_name, published_delta) {
  got <- worked_pools %>% filter(target_ancestry == ancestry, pool == pool_name)
  stopifnot(nrow(got) == 1L)
  stopifnot(abs(got$delta_hat - published_delta) < 1e-8)
}
for (anc in unique(published$target_ancestry)) {
  row <- published %>% filter(target_ancestry == anc)
  check_pool(anc, "All 4 pairs", row$delta_all)
  check_pool(anc, "2 matched pairs", row$delta_matched)
}

write_csv(chd %>% select(
  target_ancestry, pgs_id, delta, delta_se, eur_classes, tgt_classes,
  covariate_matched, decision, sampleset_id, n_individuals
), file.path(root_dir, "worked_chd_pairs.csv"))
write_csv(worked_pools, file.path(root_dir, "worked_chd_pools.csv"))

# -----------------------------------------------------------------
# 1. Which pairs match
# -----------------------------------------------------------------
class_rows <- chd %>%
  filter(target_ancestry == "African") %>%
  arrange(desc(covariate_matched), pgs_id) %>%
  mutate(
    y = factor(pgs_id, levels = rev(pgs_id)),
    why = if_else(
      covariate_matched,
      "One class on the European side, and the target side uses that same class.",
      "The European side mixes CoLaus (age, sex) with the eMERGE network (clinical covariates)."
    )
  )

p1 <- ggplot(class_rows, aes(y = y)) +
  geom_text(aes(x = 0.02, label = pgs_id), hjust = 0, fontface = "bold", family = "DejaVu Sans", size = 3.6) +
  geom_text(aes(x = 0.22, label = eur_short, colour = decision), hjust = 0, family = "DejaVu Sans", size = 3.5) +
  geom_text(aes(x = 0.62, label = tgt_short), hjust = 0, family = "DejaVu Sans", size = 3.5, colour = "#222222") +
  geom_text(aes(x = 0.92, label = decision, colour = decision), hjust = 0, fontface = "bold", family = "DejaVu Sans", size = 3.5) +
  annotate("text", x = 0.02, y = 4.55, label = "Score", hjust = 0, family = "DejaVu Sans", size = 3.1, colour = "#666666") +
  annotate("text", x = 0.22, y = 4.55, label = "European class", hjust = 0, family = "DejaVu Sans", size = 3.1, colour = "#666666") +
  annotate("text", x = 0.62, y = 4.55, label = "Target class", hjust = 0, family = "DejaVu Sans", size = 3.1, colour = "#666666") +
  scale_colour_manual(values = c(Kept = col_kept, Dropped = "#9A3B3B"), guide = "none") +
  scale_x_continuous(limits = c(0, 1.15), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.1))) +
  labs(
    title = "1.  A pair stays in only when both sides have the same single covariate class",
    subtitle = "Coronary heart disease (incident and prevalent), multi-ancestry training, African target.\nThe Hispanic or Latin American cell uses the same four scores and the same keep-or-drop decisions.",
    caption = "Classes come from the Catalog field 'Covariates Included in the Model': None (PRS only), Age/sex/PCs/technical, or Clinical or other.\nPGS000013 and PGS000018 are evaluated in Europeans both in CoLaus and in the multi-site eMERGE network, so the European side is not one class.\nThe target side of every score is the clinical eMERGE covariate list. Matching requires one class, and the same class, on both sides."
  ) +
  theme_void(base_family = "DejaVu Sans") +
  theme(
    plot.title = element_text(face = "bold", size = 15, colour = "#1a1a1a"),
    plot.subtitle = element_text(size = 11, colour = "#333333", lineheight = 1.15),
    plot.caption = element_text(size = 9.5, colour = "#444444", hjust = 0, lineheight = 1.15),
    plot.margin = margin(12, 16, 12, 12)
  )
save_lesson(p1, "01_which_pairs_match", 11.4, 5.6)

# -----------------------------------------------------------------
# 2. The two pools
# -----------------------------------------------------------------
forest_pairs <- chd %>%
  transmute(
    target_ancestry,
    axis = pgs_id,
    kind = "pair",
    estimate = delta,
    lo = delta - 1.96 * delta_se,
    hi = delta + 1.96 * delta_se,
    decision
  )

forest_pools <- worked_pools %>%
  transmute(
    target_ancestry,
    axis = pool,
    kind = "pool",
    estimate = delta_hat,
    lo = lo,
    hi = hi,
    decision = if_else(pool == "All 4 pairs", "All", "Matched")
  )

forest <- bind_rows(forest_pairs, forest_pools) %>%
  mutate(
    target_ancestry = factor(
      target_ancestry,
      levels = c("African", "Hispanic or Latin American")
    ),
    axis = factor(axis, levels = rev(c(
      "PGS000011", "PGS000200", "PGS000013", "PGS000018",
      "All 4 pairs", "2 matched pairs"
    )))
  )

p2 <- ggplot(forest, aes(x = estimate, y = axis)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(xmin = lo, xmax = hi, colour = decision),
    orientation = "y", width = 0.18, linewidth = 0.7
  ) +
  geom_point(
    data = function(d) filter(d, kind == "pair"),
    aes(colour = decision), size = 2.8
  ) +
  geom_point(
    data = function(d) filter(d, kind == "pool"),
    aes(colour = decision, shape = decision), size = 4.2
  ) +
  scale_colour_manual(values = c(
    Kept = col_kept, Dropped = col_drop, All = col_pool_all, Matched = col_pool_matched
  )) +
  scale_shape_manual(values = c(Kept = 16, Dropped = 16, All = 18, Matched = 17)) +
  facet_wrap(~ target_ancestry, nrow = 1) +
  labs(
    title = "2.  Dropping the mixed-class scores moves the pooled ΔAUC",
    subtitle = "Indigo pairs are kept. Grey pairs leave. The black diamond is the inverse-variance pool of all four.\nThe red triangle is the pool of the two kept scores. That is the only second point the review figure draws for these cells.",
    x = "ΔAUC (target − European) with 95% CI",
    y = NULL,
    caption = "The four standard errors within a cell are almost the same, so the pool is close to a simple mean.\nAfrican: the two dropped scores sit near 0, and the pool moves from 0.007 to 0.014. Both intervals still include 0.\nHispanic or Latin American: the pool moves from 0.017 to 0.021. Both intervals still exclude 0. Matching changes the estimate and leaves the conclusion in place."
  ) +
  theme_lesson(12) +
  theme(strip.text = element_text(face = "bold"), legend.position = "none")
save_lesson(p2, "02_two_pools", 11.2, 6.4)

# -----------------------------------------------------------------
# 3. Cells that keep a single point
# -----------------------------------------------------------------
ov <- pairs %>%
  filter(
    trait_label == "Epithelial non-mucinous ovarian cancer",
    trained_bucket == "Train: European-only",
    target_ancestry == "East Asian"
  )
stopifnot(nrow(ov) == 2L, all(ov$covariate_matched), all(ov$eur_classes == "None (PRS only)"))
ov_all <- pool_delta(ov$delta, ov$delta_se)
ov_matched <- pool_delta(ov$delta[ov$covariate_matched], ov$delta_se[ov$covariate_matched])
stopifnot(abs(ov_all$delta_hat - ov_matched$delta_hat) < 1e-12)

bc <- pairs %>%
  filter(
    trait_label == "Breast Carcinoma",
    trained_bucket == "Train: European-only",
    target_ancestry == "African"
  )
stopifnot(nrow(bc) == 4L, !any(bc$covariate_matched), all(bc$either_side_mixed_classes))
bc_all <- pool_delta(bc$delta, bc$delta_se)
bc_published <- cells %>%
  filter(trait_label == "Breast Carcinoma", trained_bucket == "Train: European-only", target_ancestry == "African")
stopifnot(is.na(bc_published$delta_matched))
stopifnot(abs(bc_all$delta_hat - bc_published$delta_all) < 1e-8)

n_identical <- sum(cells$n_pairs_covariate_matched == cells$n_pairs)
n_reduced <- sum(cells$n_pairs_covariate_matched > 0 & cells$n_pairs_covariate_matched < cells$n_pairs)
n_empty <- sum(cells$n_pairs_covariate_matched == 0)
stopifnot(n_identical == 24L, n_reduced == 2L, n_empty == 2L)

one_point <- bind_rows(
  ov %>% transmute(
    panel = "Ovarian cancer, East Asian\nboth pairs already match",
    axis = pgs_id,
    estimate = delta, lo = delta - 1.96 * delta_se, hi = delta + 1.96 * delta_se,
    kind = "pair"
  ),
  tibble(
    panel = "Ovarian cancer, East Asian\nboth pairs already match",
    axis = "IVW of both pairs",
    estimate = ov_all$delta_hat, lo = ov_all$lo, hi = ov_all$hi,
    kind = "pool"
  ),
  bc %>% transmute(
    panel = "Breast carcinoma, African\nno pair has a single class",
    axis = pgs_id,
    estimate = delta, lo = delta - 1.96 * delta_se, hi = delta + 1.96 * delta_se,
    kind = "pair"
  ),
  tibble(
    panel = "Breast carcinoma, African\nno pair has a single class",
    axis = "IVW of all four",
    estimate = bc_all$delta_hat, lo = bc_all$lo, hi = bc_all$hi,
    kind = "pool"
  )
) %>%
  mutate(axis = factor(axis, levels = rev(unique(axis))))

p3 <- ggplot(one_point, aes(x = estimate, y = axis)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.18, linewidth = 0.7, colour = col_kept) +
  geom_point(data = function(d) filter(d, kind == "pair"), size = 2.6, colour = col_kept) +
  geom_point(data = function(d) filter(d, kind == "pool"), shape = 18, size = 4.2, colour = col_pool_all) +
  facet_wrap(~ panel, scales = "free_y", ncol = 1) +
  labs(
    title = "3.  Every other cell keeps a single point",
    subtitle = sprintf(
      "%d of 28 cells already match on every pair, so a second pool would repeat the first number.\n%d cells, both breast carcinoma, have no matched pair. The review figure keeps the full inverse-variance point and writes \"no matched pair\".",
      n_identical, n_empty
    ),
    x = "ΔAUC (target − European) with 95% CI",
    y = NULL,
    caption = "Ovarian cancer: both scores are fit with no covariates on either side. The matched pool equals the full pool.\nBreast carcinoma, African: each score mixes \"Age/sex/PCs/technical\" and \"None (PRS only)\" on both sides, so the matched set is empty.\nThe diamond is the only estimate drawn. The two coronary heart disease cells in plot 2 are the only cells with a second point."
  ) +
  theme_lesson(12) +
  theme(strip.text = element_text(face = "bold", hjust = 0))
save_lesson(p3, "03_cells_with_one_point", 10.4, 7.2)

# -----------------------------------------------------------------
# 4. Why these two retained scores are not a random-effects model
# -----------------------------------------------------------------
kept_af <- chd %>%
  filter(target_ancestry == "African", covariate_matched) %>%
  arrange(pgs_id)

pss_af <- shared_pss %>% filter(target_ancestry == "African")
stopifnot(n_distinct(kept_af$sampleset_id) == 1L)

share_df <- tibble(
  label = c(
    sprintf("%s\nΔAUC %s", kept_af$pgs_id[1], fmt(kept_af$delta[1], 3)),
    sprintf("One target sample\n%s people", fmt_n(pss_af$n_individuals)),
    sprintf("%s\nΔAUC %s", kept_af$pgs_id[2], fmt(kept_af$delta[2], 3))
  ),
  x = c(1, 2, 3),
  kind = c("score", "sample", "score")
)

p4 <- ggplot(share_df, aes(x = x, y = 1)) +
  annotate("segment", x = 1.25, xend = 1.75, y = 1, yend = 1, linewidth = 0.6, colour = "#888888") +
  annotate("segment", x = 2.25, xend = 2.75, y = 1, yend = 1, linewidth = 0.6, colour = "#888888") +
  geom_point(data = function(d) filter(d, kind == "score"), size = 16, colour = col_kept) +
  geom_point(data = function(d) filter(d, kind == "sample"), size = 16, colour = "#C4B8A5", shape = 15) +
  geom_text(aes(label = label, y = 1.28), family = "DejaVu Sans", size = 3.4, lineheight = 0.95, colour = "#1a1a1a") +
  scale_x_continuous(limits = c(0.35, 3.65)) +
  scale_y_continuous(limits = c(0.7, 1.55)) +
  labs(
    title = "4.  The two kept scores are still one target sample, so there is no random effect",
    subtitle = sprintf(
      "African target sample set %s, %s people.\nHispanic or Latin American is the same pattern: both kept scores use one shared target sample set.\nA DerSimonian-Laird model needs two or more independent ΔAUC values. Two scores on the same people are one ΔAUC with two versions.",
      pss_af$sampleset_id, fmt_n(pss_af$n_individuals)
    ),
    caption = "The review figure therefore stops at the red triangle in plot 2: the inverse-variance pool of the two matched pairs.\nIt does not add a random-effects estimate. The dependence lesson, for type 2 diabetes, is where several target sample sets remain and that model is fit."
  ) +
  theme_void(base_family = "DejaVu Sans") +
  theme(
    plot.title = element_text(face = "bold", size = 15, colour = "#1a1a1a"),
    plot.subtitle = element_text(size = 11, colour = "#333333", lineheight = 1.15),
    plot.caption = element_text(size = 9.5, colour = "#444444", hjust = 0, lineheight = 1.15),
    plot.margin = margin(12, 16, 12, 12)
  )
save_lesson(p4, "04_one_sample_set", 11.2, 4.8)

message("Covariate-match illustrations: ", normalizePath(root_dir))
