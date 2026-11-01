#!/usr/bin/env Rscript

# Didactic figures for the covariance model of shared samples.
# Run from the repository root after analysis/09_review_sensitivity.R.
#
# Words in titles and captions: evaluation, score, sample set, comparison,
# European, non-European. The average of all scores is the inverse-variance
# average. Do not write cell, code, pool, cohorts left, key, shape, or bucket.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tibble)
  library(tidyr)
})

if (!file.exists("R/load.R")) {
  stop("Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

stamp <- pipeline_stamp
table_dir <- file.path("results", paste0("review_sensitivity_", stamp), "tables")
root_dir <- file.path("results", paste0("review_sensitivity_", stamp), "q_dependence_illustrations")
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

read_need <- function(name) {
  path <- file.path(table_dir, name)
  if (!file.exists(path)) stop("Missing ", path, ". Run analysis/09_review_sensitivity.R first.", call. = FALSE)
  read_csv(path, show_col_types = FALSE)
}

scores_src <- read_need("07_source_sharing_scores.csv")
cells <- read_need("09_dependence_cells.csv")
rho_tbl <- read_need("09_dependence_rho.csv")
score_fit <- read_need("09_dependence_scores.csv")

theme_lesson <- function(base = 13) {
  theme_minimal(base_size = base, base_family = "DejaVu Sans") +
    theme(
      plot.title = element_text(face = "bold", size = base + 2, colour = "#1a1a1a"),
      plot.subtitle = element_text(size = base - 1.5, colour = "#333333", lineheight = 1.15),
      plot.caption = element_text(size = base - 3.5, colour = "#444444", hjust = 0, lineheight = 1.15),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
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

train_short <- function(x) {
  ifelse(grepl("European-only", x), "European-only training", "Multi-ancestry training")
}

study_palette <- function(studies) {
  base <- c(
    "UK Biobank" = "#332288",
    "All of Us" = "#E69F00",
    "Genes & Health" = "#009E73",
    "FinnGen" = "#D55E00",
    "EB" = "#56B4E9",
    "HUNT" = "#CC79A7",
    "MCCS" = "#882255",
    "QSKIN" = "#44AA99",
    "(blank)" = "#999999"
  )
  extra <- setdiff(unique(studies), names(base))
  extra_cols <- setNames(
    rep(c("#117733", "#88CCEE", "#661100", "#AA4499", "#999933"), length.out = length(extra)),
    extra
  )
  c(base, extra_cols)
}

pick <- function(trait, bucket, ancestry) {
  list(trait = trait, bucket = bucket, ancestry = ancestry)
}

sa <- pick("Type 2 diabetes", "Train: European-only", "South Asian")
ea <- pick("Type 2 diabetes", "Train: Multi incl. EUR", "East Asian")
stroke <- pick("Stroke excluding subarachnoid hemorrhage", "Train: Multi incl. EUR", "South Asian")

subset_comp <- function(df, spec) {
  df %>%
    filter(trait_label == spec$trait, trained_bucket == spec$bucket, target_ancestry == spec$ancestry)
}

unnest_sets <- function(spec) {
  raw <- subset_comp(scores_src, spec)
  n_id <- lengths(strsplit(raw$sampleset_ids, ";", fixed = TRUE))
  n_name <- lengths(strsplit(raw$display_name, ";", fixed = TRUE))
  if (!all(n_id == n_name)) stop("Sample-set ids and study names are not aligned.", call. = FALSE)
  raw %>%
    separate_rows(sampleset_ids, display_name, sep = ";") %>%
    rename(sampleset_id = sampleset_ids, study = display_name) %>%
    mutate(study = ifelse(startsWith(study, "unnamed:"), "unnamed", study))
}

who_shares <- function(spec, title, subtitle) {
  sets <- unnest_sets(spec)
  eur <- sets %>% filter(side == "European") %>% distinct(sampleset_id, study)
  tgt <- sets %>% filter(side == "Target") %>% distinct(sampleset_id, study)
  sc <- sets %>% distinct(pgs_id)
  place <- function(ids) {
    tibble(id = ids, y = rev(seq_along(ids)))
  }
  eur_y <- place(eur$sampleset_id) %>% rename(sampleset_id = id, y_set = y)
  tgt_y <- place(tgt$sampleset_id) %>% rename(sampleset_id = id, y_set = y)
  sc_y <- place(sc$pgs_id) %>% rename(pgs_id = id, y_score = y)
  links <- sets %>%
    distinct(pgs_id, side, sampleset_id, study) %>%
    left_join(sc_y, by = "pgs_id") %>%
    left_join(eur_y %>% rename(y_eur = y_set), by = "sampleset_id") %>%
    left_join(tgt_y %>% rename(y_tgt = y_set), by = "sampleset_id") %>%
    mutate(
      x = if_else(side == "European", 0, 2),
      xend = 1,
      y = if_else(side == "European", y_eur, y_tgt),
      yend = y_score
    )
  pal <- study_palette(links$study)
  nodes <- bind_rows(
    eur %>% left_join(eur_y, by = "sampleset_id") %>% transmute(x = 0, y = y_set, label = paste(sampleset_id, study, sep = "\n"), study),
    sc_y %>% transmute(x = 1, y = y_score, label = pgs_id, study = "Score"),
    tgt %>% left_join(tgt_y, by = "sampleset_id") %>% transmute(x = 2, y = y_set, label = paste(sampleset_id, study, sep = "\n"), study)
  )
  ggplot() +
    geom_segment(
      data = links,
      aes(x = x, xend = xend, y = y, yend = yend, colour = study),
      linewidth = 0.6, alpha = 0.85
    ) +
    geom_label(
      data = nodes %>% filter(study != "Score"),
      aes(x = x, y = y, label = label, fill = study),
      colour = "white", size = 2.6, linewidth = 0, label.padding = unit(0.15, "lines"),
      fontface = "bold"
    ) +
    geom_label(
      data = nodes %>% filter(study == "Score"),
      aes(x = x, y = y, label = label),
      fill = "#1a1a1a", colour = "white", size = 2.7, linewidth = 0,
      label.padding = unit(0.15, "lines")
    ) +
    scale_colour_manual(values = pal, breaks = setdiff(names(pal), "Score")) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(limits = c(-0.35, 2.35), breaks = c(0, 1, 2), labels = c("European sample", "Score", "Non-European sample")) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL, caption = "Colour is the source study. Two sample-set ids of one colour are the same people.") +
    theme_lesson(12) +
    theme(axis.text.y = element_blank(), legend.position = "bottom")
}

matrix_stages <- function(spec) {
  df <- subset_comp(score_fit, spec) %>% arrange(eur_cluster, pgs_id)
  cov <- difference_covariance(df$eur_se, df$tgt_se, df$eur_cluster, df$tgt_cluster, df$pgs_id, rho = 0.6)
  off <- function(M) {
    diag(M) <- 0
    M
  }
  stage1 <- diag(nrow(cov$V)) * diag(cov$V)
  stage2 <- stage1 + off(cov$V_eur)
  stage3 <- stage2 + off(cov$V_tgt)
  if (max(abs(stage3 - cov$V)) > 1e-10) stop("The three covariance stages do not add to V.", call. = FALSE)
  ids <- df$pgs_id
  melt <- function(M, stage) {
    dimnames(M) <- list(ids, ids)
    as.data.frame(M) %>%
      mutate(row = ids) %>%
      pivot_longer(-row, names_to = "col", values_to = "value") %>%
      mutate(stage = stage, row = factor(row, levels = rev(ids)), col = factor(col, levels = ids))
  }
  bind_rows(
    melt(stage1, "1.  Diagonal only"),
    melt(stage2, "2.  Plus the European block"),
    melt(stage3, "3.  Plus the non-European block")
  )
}

heatmap_plot <- function(long, title, subtitle) {
  ggplot(long, aes(x = col, y = row, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient(
      low = "#F7F7F7", high = "#332288", trans = "sqrt", name = NULL,
      guide = guide_colourbar(barwidth = unit(8, "cm"))
    ) +
    facet_wrap(~stage, nrow = 1) +
    labs(
      title = title, subtitle = subtitle, x = NULL, y = NULL,
      caption = "Each square is the covariance of two scores' ancestry differences. Colour is on a square-root scale. The European block is faint where those standard errors are small. The non-European block is the wide square. A score whose label is not UK Biobank stays out of that square. \u03c1 = 0.6."
    ) +
    theme_lesson(12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), axis.text.y = element_text(size = 7))
}

score_forest <- function(spec, title, subtitle, caption) {
  df <- subset_comp(score_fit, spec) %>% arrange(desc(delta))
  cell <- subset_comp(cells, spec)
  rows <- bind_rows(
    tibble(
      label = "Average of all scores",
      estimate = cell$delta_fixed, lo = cell$lo_fixed, hi = cell$hi_fixed, kind = "Average of all scores"
    ),
    tibble(
      label = "Covariance model",
      estimate = cell$delta_cov, lo = cell$lo_cov, hi = cell$hi_cov, kind = "Covariance model"
    ),
    if (isTRUE(cell$robust_reported)) {
      tibble(
        label = "Robust interval",
        estimate = cell$delta_robust, lo = cell$lo_robust, hi = cell$hi_robust, kind = "Robust interval"
      )
    },
    tibble(
      label = df$pgs_id,
      estimate = df$delta,
      lo = df$delta - 1.96 * df$delta_se,
      hi = df$delta + 1.96 * df$delta_se,
      kind = "Score"
    )
  ) %>%
    mutate(
      kind = factor(kind, levels = c("Average of all scores", "Covariance model", "Robust interval", "Score")),
      label = factor(label, levels = rev(label))
    )
  pal <- c(
    "Average of all scores" = "#222222",
    "Covariance model" = "#332288",
    "Robust interval" = "#01665E",
    "Score" = "#B0B0B0"
  )
  ggplot(rows, aes(x = estimate, y = label, colour = kind)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y", linewidth = 0.7) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = pal) +
    labs(title = title, subtitle = subtitle, x = "Non-European AUC \u2212 European AUC, with 95% interval", y = NULL, caption = caption) +
    theme_lesson(12)
}

rho_plot <- function(spec) {
  df <- subset_comp(rho_tbl, spec) %>%
    mutate(rho_lab = factor(sprintf("\u03c1 = %s", rho), levels = sprintf("\u03c1 = %s", rho)))
  ggplot(df, aes(x = delta_hat, y = rho_lab)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, orientation = "y", linewidth = 0.8, colour = "#332288") +
    geom_point(size = 3, colour = "#332288") +
    labs(
      title = "4.  The same comparison at four values of \u03c1",
      subtitle = "Type 2 diabetes, South Asian participants, European-only training.",
      x = "Non-European AUC \u2212 European AUC, with 95% interval", y = NULL,
      caption = "With one sharing group, the interval depends on the assumed \u03c1; we show 0.3\u20130.9. \u03c1 = 0 still estimates a between-score variance."
    ) +
    theme_lesson(12) +
    theme(legend.position = "none")
}

# -----------------------------------------------------------------
# 1. Who shares what
# -----------------------------------------------------------------
p_who <- who_shares(
  sa,
  "1.  Who shares a sample, type 2 diabetes, South Asian participants",
  "European-only training. Scores sit in the middle. A line joins each score to its European sample and its non-European sample."
)
save_lesson(p_who, "01_who_shares_south_asian", width = 12.4, height = 8.6)

# -----------------------------------------------------------------
# 2. Covariance matrix
# -----------------------------------------------------------------
heat <- matrix_stages(sa)
write_csv(heat, file.path(root_dir, "02_covariance_south_asian.csv"))
p_heat <- heatmap_plot(
  heat,
  "2.  Covariance of the ancestry differences",
  "Type 2 diabetes, South Asian participants, European-only training. Built from two calls to vcalc."
)
save_lesson(p_heat, "02_covariance_south_asian", width = 12.8, height = 7.2)

# -----------------------------------------------------------------
# 3. Forest
# -----------------------------------------------------------------
p_for <- score_forest(
  sa,
  "3.  What the covariance does to the interval",
  "Type 2 diabetes, South Asian participants, European-only training. This comparison is one sharing group, so the robust interval is not reported.",
  "Grey intervals are one score. Black is the average of all scores. Blue is the covariance model at \u03c1 = 0.6. The robust interval needs degrees of freedom of at least 4."
)
save_lesson(p_for, "03_forest_south_asian", width = 10.5, height = 7.4)

# -----------------------------------------------------------------
# 4. rho
# -----------------------------------------------------------------
p_rho <- rho_plot(sa)
save_lesson(p_rho, "04_rho_south_asian", width = 9.2, height = 4.6)
write_csv(subset_comp(rho_tbl, sa), file.path(root_dir, "04_rho_south_asian.csv"))

# -----------------------------------------------------------------
# 5. Two other patterns
# -----------------------------------------------------------------
p_ea <- score_forest(
  ea,
  "European side only",
  "Type 2 diabetes, East Asian participants, multi-ancestry training.",
  "Two scores. The between-score variance is unstable and is not interpreted."
)
p_st <- score_forest(
  stroke,
  "Both sides, same studies for every score",
  "Stroke, South Asian participants, multi-ancestry training.",
  "The covariance model keeps every score. One sharing group, so the interval depends on \u03c1."
)
p_pair <- (p_ea | p_st) +
  plot_annotation(
    title = "5.  Two other patterns of shared samples",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, family = "DejaVu Sans", colour = "#1a1a1a")
    )
  )
save_lesson(p_pair, "05_two_other_patterns", width = 13.4, height = 6.8)

# -----------------------------------------------------------------
# 6. Dumbbell of every shared comparison
# -----------------------------------------------------------------
shared <- cells %>%
  filter(goes_through_covariance_model) %>%
  mutate(
    label = sprintf("%s — %s, %s", trait_label, target_ancestry, train_short(trained_bucket)),
    conclusion = case_when(
      linked_changes_conclusion ~ "Linked-group bound changes the conclusion",
      sig_fixed != sig_cov ~ "Covariance model changes the conclusion",
      TRUE ~ "Conclusion unchanged"
    ),
    conclusion = factor(conclusion, levels = c(
      "Conclusion unchanged",
      "Covariance model changes the conclusion",
      "Linked-group bound changes the conclusion"
    ))
  ) %>%
  arrange(delta_fixed) %>%
  mutate(label = factor(label, levels = label))

write_csv(
  shared %>%
    select(trait_label, trained_bucket, target_ancestry, n_pairs, share_side, n_sharing_groups,
           delta_fixed, lo_fixed, hi_fixed, sig_fixed,
           delta_cov, lo_cov, hi_cov, sig_cov, width_ratio_shared, width_ratio_cov,
           delta_cov_linked, lo_cov_linked, hi_cov_linked, sig_cov_linked,
           linked_changes_conclusion, robust_df, robust_reported, conclusion),
  file.path(root_dir, "06_dumbbell_shared.csv")
)

dumb_long <- bind_rows(
  shared %>% transmute(label, conclusion, estimate = delta_fixed, lo = lo_fixed, hi = hi_fixed, which = "Average of all scores"),
  shared %>% transmute(label, conclusion, estimate = delta_cov, lo = lo_cov, hi = hi_cov, which = "Covariance model")
)
p_dumb <- ggplot(dumb_long, aes(y = label, colour = which)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.28, linewidth = 0.7, position = position_dodge(width = 0.55)) +
  geom_point(aes(x = estimate, shape = conclusion), size = 2.2, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c("Average of all scores" = "#222222", "Covariance model" = "#332288")) +
  scale_shape_manual(values = c(
    "Conclusion unchanged" = 16,
    "Covariance model changes the conclusion" = 17,
    "Linked-group bound changes the conclusion" = 18
  )) +
  labs(
    title = "6.  Every comparison whose scores share people",
    subtitle = "Black is the average of all scores. Blue is the covariance model at \u03c1 = 0.6. A triangle marks a comparison whose interval changes whether it excludes zero.",
    x = "Non-European AUC \u2212 European AUC, with 95% interval", y = NULL,
    caption = "The width ratio uses the delta-method inverse-variance interval. The linked-group bound puts \u03c1 on every pair that shares any study. Comparisons with no shared people are not drawn."
  ) +
  guides(
    colour = guide_legend(nrow = 1, order = 1),
    shape = guide_legend(nrow = 1, order = 2)
  ) +
  theme_lesson(11) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.box = "vertical",
    plot.margin = margin(12, 18, 28, 12)
  )
save_lesson(p_dumb, "06_dumbbell_shared", width = 12.2, height = 11.4)

message("Dependence illustrations done: ", normalizePath(root_dir))
