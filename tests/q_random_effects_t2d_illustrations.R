#!/usr/bin/env Rscript

# Walk through the metafor pipeline for an ancestry difference.
# One comparison is drawn at every step: type 2 diabetes, South Asian
# participants, European-only training. Other comparisons appear only
# where that example cannot show the step.
#
# Run from the repository root, after analysis/09_review_sensitivity.R:
#   Rscript tests/q_random_effects_t2d_illustrations.R
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
table_dir <- file.path("results", paste0("review_sensitivity_", stamp), "tables")
root_dir <- file.path("results", paste0("review_sensitivity_", stamp), "q_random_effects_t2d_illustrations")
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
old_files <- list.files(root_dir, pattern = "\\.(png|pdf|csv)$", full.names = TRUE)
if (length(old_files)) unlink(old_files)

read_need <- function(name) {
  path <- file.path(table_dir, name)
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run analysis/09_review_sensitivity.R first.", call. = FALSE)
  }
  read_csv(path, show_col_types = FALSE)
}

scores <- read_need("09_dependence_scores.csv")
cells <- read_need("09_dependence_cells.csv")
rho_tbl <- read_need("09_dependence_rho.csv")

rho_main <- 0.6
min_robust_df <- 4

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

wrap_cap <- function(text, width = 104) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  paste(vapply(lines, function(line) {
    paste(strwrap(line, width = width), collapse = "\n")
  }, character(1)), collapse = "\n")
}

annotation_theme <- theme(
  plot.title = element_text(face = "bold", size = 15, family = "DejaVu Sans", colour = "#1a1a1a"),
  plot.subtitle = element_text(size = 11, family = "DejaVu Sans", colour = "#333333", lineheight = 1.15),
  plot.caption = element_text(size = 9, family = "DejaVu Sans", colour = "#444444", hjust = 0, lineheight = 1.15)
)

save_lesson <- function(plot, filename, width, height) {
  png_path <- file.path(root_dir, paste0(filename, ".png"))
  pdf_path <- file.path(root_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 170, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, bg = "white", device = grDevices::cairo_pdf, limitsize = FALSE)
  message("Wrote ", png_path)
}

write_numbers <- function(df, filename) {
  path <- file.path(root_dir, paste0(filename, ".csv"))
  write_csv(df, path)
  invisible(path)
}

train_short <- function(x) {
  ifelse(grepl("European-only", x), "European-only training", "Multi-ancestry training")
}

comp_phrase <- function(trait, ancestry, bucket) {
  sprintf("%s, %s participants, %s", trait, ancestry, train_short(bucket))
}

subset_comp <- function(df, trait, bucket, ancestry) {
  df %>%
    filter(
      .data$trait_label == trait,
      .data$trained_bucket == bucket,
      .data$target_ancestry == ancestry
    )
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

split_tokens <- function(x) {
  lapply(strsplit(as.character(x), ";", fixed = TRUE), function(z) unique(z[nzchar(z)]))
}

# -----------------------------------------------------------------
# Worked comparison, and the two comparisons the example cannot show
# -----------------------------------------------------------------

sa <- subset_comp(scores, "Type 2 diabetes", "Train: European-only", "South Asian") %>%
  arrange(.data$pgs_id)
if (nrow(sa) != 13L) stop("The worked comparison does not have 13 scores.", call. = FALSE)

sa_cell <- subset_comp(cells, "Type 2 diabetes", "Train: European-only", "South Asian")
if (nrow(sa_cell) != 1L) stop("The worked comparison is missing from the comparison table.", call. = FALSE)

fit_one <- function(df, cluster_eur, cluster_tgt, sigma2 = NULL) {
  fit_difference_model(
    df$delta, df$eur_se, df$tgt_se, cluster_eur, cluster_tgt, df$pgs_id,
    rho = rho_main, sigma2 = sigma2
  )
}

fit_plain_rho0 <- fit_difference_model(
  sa$delta, sa$eur_se, sa$tgt_se, sa$eur_cluster, sa$tgt_cluster, sa$pgs_id,
  rho = 0, sigma2 = 0
)
fit_corr <- fit_one(sa, sa$eur_cluster, sa$tgt_cluster, sigma2 = 0)
fit_full <- fit_one(sa, sa$eur_cluster, sa$tgt_cluster)
fit_linked <- fit_one(sa, sa$eur_group, sa$tgt_group)

recover <- function(fit) abs(sum(fit$weight * sa$delta) - fit$estimate) < 1e-8
if (!recover(fit_plain_rho0) || !recover(fit_corr) || !recover(fit_full) || !recover(fit_linked)) {
  stop("Row-sum weights do not recover the estimate.", call. = FALSE)
}
if (abs(fit_full$estimate - sa_cell$delta_cov) > 1e-8 || abs(fit_full$sigma2 - sa_cell$sigma2_cov) > 1e-8) {
  stop("The worked comparison does not match delta_cov in 09_dependence_cells.csv.", call. = FALSE)
}
if (abs(fit_corr$estimate - sa_cell$delta_shared) > 1e-8) {
  stop("The correlation-only fit does not match the saved comparison.", call. = FALSE)
}
if (abs(fit_linked$estimate - sa_cell$delta_cov_linked) > 1e-8) {
  stop("The linked fit does not match the saved comparison.", call. = FALSE)
}
if (fit_corr$weight[["PGS004106"]] < 0.7 || fit_corr$weight[["PGS001295"]] >= 0) {
  stop("Correlation-only weights are not the row-sum shares.", call. = FALSE)
}

shared_cells <- cells %>% filter(.data$goes_through_covariance_model)
n_comp <- nrow(cells)
n_shared <- nrow(shared_cells)
n_one_group <- sum(shared_cells$n_sharing_groups == 1L)
n_unstable <- sum(shared_cells$sigma2_unstable)
n_robust <- sum(shared_cells$robust_reported %in% TRUE)

bc <- subset_comp(scores, "Breast cancer", "Train: European-only", "Hispanic or Latin American") %>%
  arrange(.data$pgs_id)
bc_cell <- subset_comp(cells, "Breast cancer", "Train: European-only", "Hispanic or Latin American")
if (nrow(bc) != 2L || nrow(bc_cell) != 1L || bc_cell$sigma2_cov > 1e-8) {
  stop("The two-score comparison is not stuck at a between-score variance of 0.", call. = FALSE)
}
fit_bc <- fit_one(bc, bc$eur_cluster, bc$tgt_cluster)

# -----------------------------------------------------------------
# 1. Pipeline map
# -----------------------------------------------------------------

fmt_est <- function(x) formatC(x, format = "f", digits = 3)
sigma_txt <- formatC(fit_full$sigma2, format = "e", digits = 2)
robust_row <- shared_cells %>% filter(is.finite(.data$robust_df))
if (nrow(robust_row) != 1L) stop("Expected one comparison with robust degrees of freedom.", call. = FALSE)
df_txt <- formatC(robust_row$robust_df, format = "f", digits = 1)

map_df <- tibble(
  y = 8:1,
  title = c(
    "Evaluations",
    "Score pair",
    "Labels",
    "Two calls to vcalc",
    "rma.mv",
    "Estimate",
    "robust",
    "Figure 09"
  ),
  body = c(
    sprintf("%d comparisons. This worked comparison has %d scores.", n_comp, nrow(sa)),
    "An ancestry difference, a European standard error, and a non-European standard error.",
    sprintf(
      "Sample sets become study names, then one label on each side.\nThis comparison is %d sharing group.",
      length(unique(sa$sharing_group))
    ),
    "A European covariance and a non-European covariance, at \u03c1 = 0.6. Their sum is V.",
    sprintf("Between-score variance by REML. Here \u03c3\u00b2 = %s.", sigma_txt),
    sprintf(
      "This comparison's estimate is %s.\n%d comparisons with shared people take this path. The other %d keep the average of all scores.",
      fmt_est(fit_full$estimate), n_shared, n_comp - n_shared
    ),
    sprintf(
      "Reported when the degrees of freedom are at least %d.\n%d of the %d have one sharing group, so the interval is the model interval.\n%s has %d sharing groups and %s degrees of freedom.\nBelow that threshold, the reported interval is the model interval.",
      min_robust_df, n_one_group, n_shared,
      comp_phrase(robust_row$trait_label, robust_row$target_ancestry, robust_row$trained_bucket),
      robust_row$n_sharing_groups, df_txt
    ),
    "Draws the average of all scores and this model at \u03c1 = 0.6."
  )
)
write_numbers(
  map_df %>% mutate(n_comparisons = n_comp, n_shared = n_shared, n_one_group = n_one_group,
                    n_unstable = n_unstable, n_robust_reported = n_robust),
  "01_pipeline_map"
)

p_map <- ggplot(map_df) +
  geom_rect(
    aes(xmin = 0, xmax = 1, ymin = .data$y - 0.46, ymax = .data$y + 0.46),
    fill = "#F4F7FB", colour = "#332288", linewidth = 0.4
  ) +
  geom_segment(
    data = map_df %>% filter(.data$y > 1),
    aes(x = 0.5, xend = 0.5, y = .data$y - 0.46, yend = .data$y - 0.54),
    arrow = grid::arrow(length = grid::unit(0.12, "cm")),
    colour = "#332288", linewidth = 0.4
  ) +
  geom_text(aes(x = 0.03, y = .data$y + 0.22, label = .data$title), hjust = 0, fontface = "bold", size = 4.2, family = "DejaVu Sans") +
  geom_text(aes(x = 0.03, y = .data$y - 0.14, label = .data$body), hjust = 0, vjust = 0.5, size = 3.1, lineheight = 1.02, family = "DejaVu Sans") +
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 8.6), expand = c(0, 0)) +
  labs(
    title = "1.  From evaluations to the estimate",
    subtitle = "Type 2 diabetes, South Asian participants, European-only training. Counts are every comparison.",
    caption = wrap_cap(sprintf(
      "%d of the %d comparisons with shared people have 2 or 3 scores, so the between-score variance is flagged unstable.",
      n_unstable, n_shared
    ))
  ) +
  theme_void(base_family = "DejaVu Sans") +
  theme(
    plot.title = element_text(face = "bold", size = 15, colour = "#1a1a1a", hjust = 0),
    plot.subtitle = element_text(size = 11, colour = "#333333", hjust = 0),
    plot.caption = element_text(size = 9, colour = "#444444", hjust = 0),
    plot.margin = margin(12, 18, 12, 12)
  )
save_lesson(p_map, "01_pipeline_map", width = 11.2, height = 12.4)

# -----------------------------------------------------------------
# 2. What metafor receives
# -----------------------------------------------------------------

score_y <- sa %>%
  mutate(w = fit_corr$weight[as.character(.data$pgs_id)]) %>%
  arrange(desc(.data$w)) %>%
  pull(.data$pgs_id)
sa_view <- sa %>%
  mutate(
    var_eur = .data$eur_se^2,
    var_tgt = .data$tgt_se^2,
    lo = .data$delta - 1.96 * .data$delta_se,
    hi = .data$delta + 1.96 * .data$delta_se,
    y_lab = factor(.data$pgs_id, levels = rev(score_y))
  )

write_numbers(
  sa_view %>%
    transmute(
      pgs_id, delta, delta_se, eur_se, tgt_se, var_eur, var_tgt,
      eur_cluster, tgt_cluster, sharing_group,
      weight_plain = fit_plain_rho0$weight[as.character(pgs_id)],
      weight_correlation = fit_corr$weight[as.character(pgs_id)],
      weight_full = fit_full$weight[as.character(pgs_id)],
      weight_linked = fit_linked$weight[as.character(pgs_id)]
    ),
  "02_what_metafor_receives"
)

var_long <- sa_view %>%
  select("y_lab", "European" = "var_eur", "Non-European" = "var_tgt") %>%
  pivot_longer(-"y_lab", names_to = "side", values_to = "variance") %>%
  mutate(side = factor(.data$side, levels = c("European", "Non-European")))

p_labs <- ggplot(sa_view, aes(y = .data$y_lab)) +
  geom_text(aes(x = 0, label = .data$pgs_id), hjust = 0, fontface = "bold", size = 2.8, family = "DejaVu Sans") +
  geom_text(aes(x = 1.05, label = .data$eur_cluster), hjust = 0, size = 2.5, family = "DejaVu Sans") +
  geom_text(aes(x = 3.35, label = .data$tgt_cluster), hjust = 0, size = 2.5, family = "DejaVu Sans") +
  scale_x_continuous(limits = c(-0.05, 5.7), expand = c(0, 0)) +
  labs(title = "Score and both labels", x = NULL, y = NULL) +
  theme_lesson(11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

p_delta <- ggplot(sa_view, aes(x = .data$delta, y = .data$y_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = .data$lo, xmax = .data$hi), orientation = "y", width = 0.2, linewidth = 0.6, colour = "#444444") +
  geom_point(size = 2.4, colour = "#222222") +
  labs(title = "Ancestry difference", x = "Non-European AUC \u2212 European AUC", y = NULL) +
  theme_lesson(11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_var <- ggplot(var_long, aes(x = .data$variance, y = .data$y_lab, fill = .data$side)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c("European" = "#88CCEE", "Non-European" = "#332288")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Variance, standard error squared", x = NULL, y = NULL) +
  theme_lesson(11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_input <- p_labs + p_delta + p_var +
  plot_layout(widths = c(1.55, 1.05, 0.95)) +
  plot_annotation(
    title = "2.  What the model receives",
    subtitle = "Type 2 diabetes, South Asian participants, European-only training. One sharing group.\nEach row is one score. The label panel gives the score, then the European label, then the non-European label.",
    caption = wrap_cap("The non-European variance is the wider bar for every score, so that side dominates the variance of the ancestry difference.", width = 140),
    theme = annotation_theme
  )
save_lesson(p_input, "02_what_metafor_receives", width = 14.4, height = 8.4)

# -----------------------------------------------------------------
# 3. From sample sets to correlation
# -----------------------------------------------------------------

long_sets <- bind_rows(lapply(seq_len(nrow(sa)), function(i) {
  r <- sa[i, ]
  eur_ids <- strsplit(r$eur_sampleset_ids, ";", fixed = TRUE)[[1]]
  eur_nm <- strsplit(r$eur_display, ";", fixed = TRUE)[[1]]
  tgt_ids <- strsplit(r$tgt_sampleset_ids, ";", fixed = TRUE)[[1]]
  tgt_nm <- strsplit(r$tgt_display, ";", fixed = TRUE)[[1]]
  if (length(eur_ids) != length(eur_nm) || length(tgt_ids) != length(tgt_nm)) {
    stop("Sample-set ids and study names are not aligned.", call. = FALSE)
  }
  bind_rows(
    tibble(pgs_id = r$pgs_id, side = "European", sampleset_id = eur_ids, study = eur_nm),
    tibble(pgs_id = r$pgs_id, side = "Target", sampleset_id = tgt_ids, study = tgt_nm)
  )
})) %>%
  mutate(study = ifelse(startsWith(.data$study, "unnamed:"), "unnamed", .data$study))

place <- function(ids) tibble(id = ids, y = rev(seq_along(ids)))
eur_nodes <- long_sets %>% filter(.data$side == "European") %>% distinct(.data$sampleset_id, .data$study)
tgt_nodes <- long_sets %>% filter(.data$side == "Target") %>% distinct(.data$sampleset_id, .data$study)
score_nodes <- long_sets %>% distinct(.data$pgs_id)
eur_y <- place(eur_nodes$sampleset_id) %>% rename(sampleset_id = "id", y_eur = "y")
tgt_y <- place(tgt_nodes$sampleset_id) %>% rename(sampleset_id = "id", y_tgt = "y")
sc_y <- place(score_nodes$pgs_id) %>% rename(pgs_id = "id", y_score = "y")
links <- long_sets %>%
  distinct(.data$pgs_id, .data$side, .data$sampleset_id, .data$study) %>%
  left_join(sc_y, by = "pgs_id") %>%
  left_join(eur_y, by = "sampleset_id") %>%
  left_join(tgt_y, by = "sampleset_id") %>%
  mutate(
    x = if_else(.data$side == "European", 0, 2),
    xend = 1,
    y = if_else(.data$side == "European", .data$y_eur, .data$y_tgt),
    yend = .data$y_score
  )
pal <- study_palette(links$study)
who_nodes <- bind_rows(
  eur_nodes %>% left_join(eur_y, by = "sampleset_id") %>% transmute(x = 0, y = .data$y_eur, label = paste(.data$sampleset_id, .data$study, sep = "\n"), study),
  sc_y %>% transmute(x = 1, y = .data$y_score, label = .data$pgs_id, study = "Score"),
  tgt_nodes %>% left_join(tgt_y, by = "sampleset_id") %>% transmute(x = 2, y = .data$y_tgt, label = paste(.data$sampleset_id, .data$study, sep = "\n"), study)
)
p_who <- ggplot() +
  geom_segment(
    data = links,
    aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend, colour = .data$study),
    linewidth = 0.6, alpha = 0.85
  ) +
  geom_label(
    data = who_nodes %>% filter(.data$study != "Score"),
    aes(x = .data$x, y = .data$y, label = .data$label, fill = .data$study),
    colour = "white", size = 2.5, linewidth = 0, label.padding = unit(0.15, "lines"), fontface = "bold"
  ) +
  geom_label(
    data = who_nodes %>% filter(.data$study == "Score"),
    aes(x = .data$x, y = .data$y, label = .data$label),
    fill = "#1a1a1a", colour = "white", size = 2.6, linewidth = 0, label.padding = unit(0.15, "lines")
  ) +
  scale_colour_manual(values = pal, breaks = setdiff(names(pal), "Score")) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = c(-0.35, 2.35), breaks = c(0, 1, 2), labels = c("European sample set", "Score", "Non-European sample set")) +
  labs(title = "Who shares a sample set", x = NULL, y = NULL) +
  theme_lesson(12) +
  theme(axis.text.y = element_blank(), legend.position = "bottom")

matrix_order <- sa %>%
  mutate(ord = case_when(
    .data$tgt_cluster == "UK Biobank" & .data$eur_cluster == "All of Us" ~ 1L,
    .data$tgt_cluster == "UK Biobank" ~ 2L,
    TRUE ~ 3L
  )) %>%
  arrange(.data$ord, .data$pgs_id)
if (!identical(as.integer(table(matrix_order$ord)[c("1", "2", "3")]), c(11L, 1L, 1L))) {
  stop("The worked comparison does not split into the eleven, PGS001295, and PGS004106.", call. = FALSE)
}
ids <- matrix_order$pgs_id

pair_grid <- function(same) {
  n <- length(ids)
  rows <- vector("list", n * n)
  k <- 0L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      k <- k + 1L
      val <- if (i == j) "1" else if (isTRUE(same(i, j))) "0.6" else "0"
      rows[[k]] <- tibble(row = ids[i], col = ids[j], value = val)
    }
  }
  bind_rows(rows) %>%
    mutate(
      row = factor(.data$row, levels = rev(ids)),
      col = factor(.data$col, levels = ids)
    )
}

eur_tokens <- split_tokens(matrix_order$eur_display)
tgt_tokens <- split_tokens(matrix_order$tgt_display)
eur_pss <- split_tokens(matrix_order$eur_sampleset_ids)
tgt_pss <- split_tokens(matrix_order$tgt_sampleset_ids)
eur_exact <- pair_grid(function(i, j) matrix_order$eur_cluster[i] == matrix_order$eur_cluster[j])
eur_any <- pair_grid(function(i, j) length(intersect(eur_tokens[[i]], eur_tokens[[j]])) > 0)
eur_group <- pair_grid(function(i, j) matrix_order$eur_group[i] == matrix_order$eur_group[j])
tgt_exact <- pair_grid(function(i, j) matrix_order$tgt_cluster[i] == matrix_order$tgt_cluster[j])
tgt_any <- pair_grid(function(i, j) length(intersect(tgt_tokens[[i]], tgt_tokens[[j]])) > 0)
tgt_group <- pair_grid(function(i, j) matrix_order$tgt_group[i] == matrix_order$tgt_group[j])
if (!identical(eur_any$value, eur_group$value) || !identical(tgt_any$value, tgt_group$value)) {
  stop("Sharing any study and the chained groups differ in this comparison.", call. = FALSE)
}

pss_differ <- function(a, b) length(setdiff(a, b)) > 0 || length(setdiff(b, a)) > 0
mark_pairs <- function(exact, any, pss) {
  exact %>%
    transmute(row, col, exact = .data$value) %>%
    left_join(any %>% transmute(row, col, any = .data$value), by = c("row", "col")) %>%
    mutate(
      i = match(as.character(.data$row), ids),
      j = match(as.character(.data$col), ids),
      differs = .data$exact != .data$any & .data$i != .data$j,
      same_label_other_set = .data$exact == "0.6" & .data$i != .data$j & mapply(pss_differ, pss[.data$i], pss[.data$j])
    )
}
eur_mark <- mark_pairs(eur_exact, eur_any, eur_pss)
tgt_mark <- mark_pairs(tgt_exact, tgt_any, tgt_pss)
if (sum(tgt_mark$differs) != 24L || sum(eur_mark$differs) != 2L) {
  stop("The differing pairs are not the twelve non-European pairs and the one European pair.", call. = FALSE)
}
if (sum(tgt_mark$same_label_other_set) != 22L) {
  stop("PGS001295 does not agree with the eleven on the non-European label.", call. = FALSE)
}

grid_plot <- function(values, marks, title, show_green) {
  hi <- marks %>%
    filter(.data$differs | (show_green & .data$same_label_other_set)) %>%
    mutate(mark = case_when(
      .data$differs ~ "Disagrees",
      TRUE ~ "Same label, different sample set"
    ))
  ggplot(values, aes(x = .data$col, y = .data$row, fill = .data$value)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_tile(
      data = hi,
      aes(colour = .data$mark),
      fill = NA, linewidth = 0.9, show.legend = FALSE
    ) +
    scale_fill_manual(values = c("0" = "#F7F7F7", "0.6" = "#88CCEE", "1" = "#332288"), drop = FALSE) +
    scale_colour_manual(values = c("Disagrees" = "#D55E00", "Same label, different sample set" = "#009E73")) +
    labs(title = title, x = NULL, y = NULL) +
    theme_lesson(10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 6),
      legend.position = "none"
    )
}

p_grids <- (grid_plot(eur_exact, eur_mark, "European side, same exact label", FALSE) |
  grid_plot(eur_any, eur_mark, "European side, shares any study", FALSE)) /
  (grid_plot(tgt_exact, tgt_mark, "Non-European side, same exact label", TRUE) |
    grid_plot(tgt_any, tgt_mark, "Non-European side, shares any study", FALSE))

p_share <- p_who / p_grids +
  plot_layout(heights = c(1.05, 1.7)) +
  plot_annotation(
    title = "3.  From sample sets to correlation",
    subtitle = "Type 2 diabetes, South Asian participants, European-only training. A grid entry is the correlation vcalc uses: 1 on the diagonal, 0.6 for a shared label, 0 otherwise.",
    caption = wrap_cap(paste(
      "Colour in the top panel is the study. The model treats two sample-set ids of one colour as the same people. PSS011749 and PSS004315 are different sample sets with the label UK Biobank.",
      "Orange outlines mark pairs where the two grids on that side disagree. On the non-European side those pairs are PGS004106 (Genes & Health + UK Biobank) with each of the other twelve. On the European side the pair is PGS004106 (EB + FinnGen + HUNT + UK Biobank) and PGS001295 (UK Biobank).",
      "Green outlines on the non-European exact-label grid mark PGS001295 with the eleven: the labels agree, and the sample sets do not.",
      "The linked fit uses chained groups. For this comparison those groups are the same pairs as sharing any study.",
      sep = "\n"
    ), width = 150),
    theme = annotation_theme
  )
save_lesson(p_share, "03_labels_and_correlation", width = 13.8, height = 16.5)
write_numbers(
  bind_rows(
    eur_exact %>% mutate(side = "European", rule = "Same exact label"),
    eur_any %>% mutate(side = "European", rule = "Shares any study"),
    tgt_exact %>% mutate(side = "Non-European", rule = "Same exact label"),
    tgt_any %>% mutate(side = "Non-European", rule = "Shares any study")
  ) %>% mutate(row = as.character(.data$row), col = as.character(.data$col)),
  "03_labels_and_correlation"
)

# -----------------------------------------------------------------
# 4. Building V
# -----------------------------------------------------------------

cov_sa <- difference_covariance(sa$eur_se, sa$tgt_se, sa$eur_cluster, sa$tgt_cluster, sa$pgs_id, rho_main)
if (max(abs(cov_sa$V - cov_sa$V_eur - cov_sa$V_tgt)) > 1e-10) {
  stop("The two covariances do not add to V.", call. = FALSE)
}
melt_num <- function(M, ids, panel) {
  dimnames(M) <- list(ids, ids)
  as.data.frame(M) %>%
    mutate(row = ids) %>%
    pivot_longer(-"row", names_to = "col", values_to = "value") %>%
    mutate(
      panel = panel,
      row = factor(.data$row, levels = rev(ids)),
      col = factor(.data$col, levels = ids)
    )
}
# Matrix order, not alphabetical pgs order.
cov_ord <- sa %>% slice(match(ids, .data$pgs_id))
cov_m <- difference_covariance(cov_ord$eur_se, cov_ord$tgt_se, cov_ord$eur_cluster, cov_ord$tgt_cluster, cov_ord$pgs_id, rho_main)
v_long <- bind_rows(
  melt_num(cov_m$V_eur, ids, "European covariance"),
  melt_num(cov_m$V_tgt, ids, "Non-European covariance"),
  melt_num(cov_m$V, ids, "V, the sum")
)
write_numbers(v_long %>% mutate(row = as.character(.data$row), col = as.character(.data$col)), "04_building_v")
v_lim <- c(0, max(v_long$value))

heat_v <- function(panel) {
  dat <- v_long %>% filter(.data$panel == .env$panel)
  if (panel == "European covariance") {
    hit <- dat %>% filter(as.character(.data$row) == "PGS000330", as.character(.data$col) == "PGS001295")
    if (nrow(hit) != 1L || abs(hit$value) > 1e-12) {
      stop("The European covariance panel is not the European block.", call. = FALSE)
    }
  }
  ggplot(dat, aes(x = .data$col, y = .data$row, fill = .data$value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "#F7F7F7", high = "#332288", limits = v_lim, trans = "sqrt", na.value = "#F7F7F7") +
    labs(title = panel, x = NULL, y = NULL) +
    theme_lesson(10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 6),
      legend.position = "none"
    )
}
pattern_values <- bind_rows(
  eur_exact %>% mutate(panel = "European correlation"),
  tgt_exact %>% mutate(panel = "Non-European correlation")
)
heat_pat <- function(panel) {
  dat <- pattern_values %>% filter(.data$panel == .env$panel)
  if (panel == "European correlation") {
    hit <- dat %>% filter(as.character(.data$row) == "PGS000330", as.character(.data$col) == "PGS001295")
    if (!identical(hit$value, "0")) {
      stop("The European correlation panel is not the European pattern.", call. = FALSE)
    }
  }
  ggplot(dat, aes(x = .data$col, y = .data$row, fill = .data$value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("0" = "#F7F7F7", "0.6" = "#88CCEE", "1" = "#332288"), drop = FALSE) +
    labs(title = panel, x = NULL, y = NULL) +
    theme_lesson(10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 6),
      legend.position = "none"
    )
}
se_df <- cov_ord %>%
  transmute(
    pgs_id = factor(.data$pgs_id, levels = rev(ids)),
    European = .data$eur_se,
    `Non-European` = .data$tgt_se
  ) %>%
  pivot_longer(-"pgs_id", names_to = "side", values_to = "se")
se_plot <- function(side, title) {
  dat <- se_df %>% filter(.data$side == .env$side)
  if (side == "European") {
    se_hit <- dat$se[as.character(dat$pgs_id) == "PGS004106"]
    if (length(se_hit) != 1L || se_hit > 0.002) {
      stop("The European standard-error panel is not the European side.", call. = FALSE)
    }
  }
  ggplot(dat, aes(x = 1, y = .data$pgs_id, fill = .data$se)) +
    geom_tile(colour = "white", width = 0.9, height = 0.9) +
    scale_fill_gradient(low = "#F7F7F7", high = "#332288", limits = range(se_df$se)) +
    scale_x_continuous(limits = c(0.4, 1.6), breaks = NULL) +
    labs(title = title, x = NULL, y = NULL) +
    theme_lesson(10) +
    theme(axis.text.y = element_blank(), legend.position = "none")
}

p_build <- (heat_pat("European correlation") | se_plot("European", "European standard error") | heat_v("European covariance")) /
  (heat_pat("Non-European correlation") | se_plot("Non-European", "Non-European standard error") | heat_v("Non-European covariance")) /
  (plot_spacer() | plot_spacer() | heat_v("V, the sum")) +
  plot_layout(widths = c(1.4, 0.35, 1.4)) +
  plot_annotation(
    title = "4.  Building V",
    subtitle = "Type 2 diabetes, South Asian participants, European-only training. \u03c1 = 0.6.",
    caption = wrap_cap("The correlation pattern is who shares a label: 0, 0.6, or 1. The standard error is how large that share is in the covariance. The European standard errors are small, so that block is faint. V is the sum of the two covariances. Colour within the covariances is on a square-root scale and is shared, so the sizes can be compared.", width = 155),
    theme = annotation_theme
  )
save_lesson(p_build, "04_building_v", width = 14.2, height = 14.4)

# -----------------------------------------------------------------
# 5. What the correlation is worth
# -----------------------------------------------------------------

n_scores_curve <- 1:15
curve_df <- tibble(
  n_scores = n_scores_curve,
  n_eff = n_scores_curve / (1 + (n_scores_curve - 1) * rho_main)
)
eff <- shared_cells %>%
  transmute(
    trait_label, trained_bucket, target_ancestry, n_scores = .data$n_pairs,
    n_eff = .data$n_pairs / .data$width_ratio_shared^2,
    side = case_when(
      .data$share_side == "both" ~ "Both sides",
      .data$share_side == "European" ~ "European side",
      .data$share_side == "Target" ~ "Non-European side",
      TRUE ~ .data$share_side
    ),
    worked = .data$trait_label == "Type 2 diabetes" &
      .data$trained_bucket == "Train: European-only" &
      .data$target_ancestry == "South Asian"
  )
write_numbers(bind_rows(curve_df %>% mutate(kind = "curve"), eff %>% mutate(kind = "comparison")), "05_effective_scores")
side_cols <- c("Both sides" = "#332288", "European side" = "#E69F00", "Non-European side" = "#009E73")

worked_eff <- eff %>% filter(.data$worked)
p_neff <- ggplot() +
  geom_line(data = curve_df, aes(x = .data$n_scores, y = .data$n_eff), colour = "#222222", linewidth = 0.7) +
  geom_point(data = eff, aes(x = .data$n_scores, y = .data$n_eff, colour = .data$side), size = 2.6) +
  geom_point(data = worked_eff, aes(x = .data$n_scores, y = .data$n_eff), size = 4.2, shape = 1, stroke = 1.1, colour = "#1a1a1a") +
  geom_text(
    data = worked_eff,
    aes(x = .data$n_scores, y = .data$n_eff, label = sprintf("Worked comparison, %.1f", .data$n_eff)),
    nudge_y = 0.55, hjust = 1, size = 3.2, family = "DejaVu Sans"
  ) +
  scale_colour_manual(values = side_cols) +
  scale_x_continuous(breaks = seq(2, 14, by = 2)) +
  labs(
    title = "5.  What the correlation is worth",
    subtitle = "Effective number of scores = score count \u00d7 (plain-average standard error / correlation-only standard error)\u00b2.",
    x = "Scores",
    y = "Effective number of scores",
    caption = wrap_cap("The curve is m / (1 + (m \u2212 1)\u03c1) at \u03c1 = 0.6. It is the effective number when every score shares one label on both sides, or the shared side carries all of the variance. Comparisons shared on one side sit near the score count, and that is why they fall off the curve. Points above the curve contain scores without a shared exact label. The worked comparison is one of them: PGS004106 is uncorrelated with the others.", width = 100)
  ) +
  theme_lesson(12)
save_lesson(p_neff, "05_effective_scores", width = 8.8, height = 7.2)

# -----------------------------------------------------------------
# 6. Between-score variance
# -----------------------------------------------------------------

sigma_grid <- seq(0, 4e-4, length.out = 41)
profile_ll <- function(df) {
  vapply(sigma_grid, function(s) {
    fit_difference_model(
      df$delta, df$eur_se, df$tgt_se, df$eur_cluster, df$tgt_cluster, df$pgs_id,
      rho = rho_main, sigma2 = s
    )$reml_ll
  }, numeric(1))
}
prof <- bind_rows(
  tibble(
    comparison = "Type 2 diabetes, South Asian participants, European-only training",
    sigma2 = sigma_grid,
    reml_ll = profile_ll(sa),
    fitted = fit_full$sigma2,
    fitted_ll = fit_full$reml_ll
  ),
  tibble(
    comparison = "Breast cancer, Hispanic or Latin American participants, European-only training",
    sigma2 = sigma_grid,
    reml_ll = profile_ll(bc),
    fitted = fit_bc$sigma2,
    fitted_ll = fit_bc$reml_ll
  )
)
write_numbers(prof, "06_sigma2_profile")
prof <- prof %>%
  mutate(comparison = factor(.data$comparison, levels = c(
    "Type 2 diabetes, South Asian participants, European-only training",
    "Breast cancer, Hispanic or Latin American participants, European-only training"
  )))
bc_name <- "Breast cancer, Hispanic or Latin American participants, European-only training"
boundary_lab <- tibble(
  comparison = factor(bc_name, levels = levels(prof$comparison)),
  sigma2 = 1.5e-5,
  reml_ll = max(prof$reml_ll[prof$comparison == bc_name]) - 0.08,
  label = "\u03c3\u00b2 = 0"
)
p_sig <- ggplot(prof, aes(x = .data$sigma2, y = .data$reml_ll)) +
  geom_line(linewidth = 0.8, colour = "#332288") +
  geom_vline(aes(xintercept = .data$fitted), linetype = "dashed", colour = "#D55E00") +
  geom_text(
    data = boundary_lab,
    aes(label = .data$label),
    hjust = 0, vjust = 1, size = 3.4, colour = "#D55E00", family = "DejaVu Sans"
  ) +
  facet_wrap(~comparison, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    limits = c(0, 4e-4),
    breaks = c(0, 1e-4, 2e-4, 3e-4, 4e-4),
    labels = c("0", "1e-4", "2e-4", "3e-4", "4e-4"),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "6.  Between-score variance",
    subtitle = "The curve is the REML log-likelihood with \u03c3\u00b2 held fixed. The dashed line is the fitted value.",
    x = "\u03c3\u00b2",
    y = "REML log-likelihood",
    caption = wrap_cap(sprintf(
      "The worked comparison has a peak near %s. Breast cancer, Hispanic or Latin American participants, European-only training has two scores and the curve is stuck at 0. %d of the %d comparisons have 2 or 3 scores and are flagged unstable for this variance.",
      sigma_txt, n_unstable, n_shared
    ), width = 108)
  ) +
  theme_lesson(12)
save_lesson(p_sig, "06_sigma2_profile", width = 10.5, height = 8.2)

# -----------------------------------------------------------------
# 7. Weights
# -----------------------------------------------------------------

weight_df <- bind_rows(
  tibble(pgs_id = names(fit_plain_rho0$weight), weight = as.numeric(fit_plain_rho0$weight), model = "Plain average"),
  tibble(pgs_id = names(fit_corr$weight), weight = as.numeric(fit_corr$weight), model = "Correlation only"),
  tibble(pgs_id = names(fit_full$weight), weight = as.numeric(fit_full$weight), model = "Full model"),
  tibble(pgs_id = names(fit_linked$weight), weight = as.numeric(fit_linked$weight), model = "Linked fit")
) %>%
  mutate(
    model = factor(.data$model, levels = c("Plain average", "Correlation only", "Full model", "Linked fit")),
    estimate = c(
      "Plain average" = fit_plain_rho0$estimate,
      "Correlation only" = fit_corr$estimate,
      "Full model" = fit_full$estimate,
      "Linked fit" = fit_linked$estimate
    )[as.character(.data$model)]
  )
weight_order <- weight_df %>%
  filter(.data$model == "Correlation only") %>%
  arrange(desc(.data$weight)) %>%
  pull(.data$pgs_id)
weight_df <- weight_df %>% mutate(pgs_id = factor(.data$pgs_id, levels = rev(weight_order)))
write_numbers(weight_df, "07_weights")
model_cols <- c(
  "Plain average" = "#222222",
  "Correlation only" = "#332288",
  "Full model" = "#01665E",
  "Linked fit" = "#E69F00"
)
p_w <- ggplot(weight_df, aes(x = .data$weight, y = .data$pgs_id, fill = .data$model)) +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_col(position = position_dodge(width = 0.8, reverse = TRUE), width = 0.72) +
  scale_fill_manual(values = model_cols) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    title = "7.  Why the estimate moves",
    subtitle = sprintf(
      "Plain average %s. Correlation only (\u03c1 = 0.6, \u03c3\u00b2 = 0) %s. Full model %s. Linked fit %s.",
      fmt_est(fit_plain_rho0$estimate), fmt_est(fit_corr$estimate),
      fmt_est(fit_full$estimate), fmt_est(fit_linked$estimate)
    ),
    x = "Share of the estimate",
    y = NULL,
    caption = wrap_cap("Shares are the row sums of the marginal precision matrix. They can be negative. With \u03c3\u00b2 = 0, PGS004106 does not share an exact label, so it takes most of the weight. PGS001295 shares the non-European label UK Biobank with the eleven and is pulled the other way. Estimating \u03c3\u00b2 brings PGS004106 back toward the plain average. The linked fit uses the chained groups.", width = 120)
  ) +
  theme_lesson(12)
save_lesson(p_w, "07_weights", width = 11.4, height = 8.8)

# -----------------------------------------------------------------
# 8. Where the interval width comes from
# -----------------------------------------------------------------

width_df <- shared_cells %>%
  transmute(
    label = sprintf("%s \u2014 %s, %s", .data$trait_label, .data$target_ancestry, train_short(.data$trained_bucket)),
    trait_label, trained_bucket, target_ancestry,
    plain = .data$se_fixed,
    correlation = .data$se_shared,
    full = .data$se_cov,
    worked = .data$trait_label == "Type 2 diabetes" &
      .data$trained_bucket == "Train: European-only" &
      .data$target_ancestry == "South Asian"
  ) %>%
  arrange(.data$full - .data$plain, .data$label) %>%
  mutate(label = factor(.data$label, levels = .data$label))
width_long <- width_df %>%
  pivot_longer(c("plain", "correlation", "full"), names_to = "step", values_to = "se") %>%
  mutate(step = factor(
    recode(.data$step, plain = "Plain average", correlation = "Correlation only", full = "Full model"),
    levels = c("Plain average", "Correlation only", "Full model")
  ))
write_numbers(width_long %>% mutate(label = as.character(.data$label)), "08_interval_width")
p_width <- ggplot(width_long, aes(x = .data$se, y = .data$label, colour = .data$step, group = .data$label)) +
  geom_line(colour = "grey70", linewidth = 0.4) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = model_cols[1:3]) +
  labs(
    title = "8.  Where the interval width comes from",
    subtitle = "One row is one comparison whose scores share people. Rows run from the largest increase in the standard error to the smallest.",
    x = "Standard error of the ancestry difference",
    y = NULL,
    caption = wrap_cap(sprintf(
      "Shared samples widen the interval %.1f\u2013%.1f times from the plain average to correlation only. The rest of the widening is the between-score variance. The worked comparison is type 2 diabetes, South Asian participants, European-only training.",
      min(shared_cells$width_ratio_shared), max(shared_cells$width_ratio_shared)
    ), width = 125)
  ) +
  theme_lesson(11) +
  theme(axis.text.y = element_text(size = 7.5))
save_lesson(p_width, "08_interval_width", width = 12.4, height = 11.8)

# -----------------------------------------------------------------
# 9. Rho for every comparison
# -----------------------------------------------------------------

rho_wide <- rho_tbl %>%
  mutate(rho_name = paste0("r", .data$rho)) %>%
  select("trait_label", "trained_bucket", "target_ancestry", "rho_name", "sig") %>%
  pivot_wider(names_from = "rho_name", values_from = "sig")
only_low <- rho_wide %>% filter(.data$r0 != .data$r0.3, .data$r0.3 == .data$r0.6, .data$r0.6 == .data$r0.9)
still_high <- rho_wide %>% filter(!(.data$r0.3 == .data$r0.6 & .data$r0.6 == .data$r0.9))
if (nrow(only_low) != 3L || nrow(still_high) != 2L) {
  stop("The zero-call changes across \u03c1 are not the expected five comparisons.", call. = FALSE)
}
east_eur <- still_high %>%
  filter(.data$trait_label == "Type 2 diabetes", .data$target_ancestry == "East Asian", grepl("European-only", .data$trained_bucket))
if (nrow(east_eur) != 1L || !isFALSE(east_eur$r0.3) || !isFALSE(east_eur$r0.6) || !isTRUE(east_eur$r0.9)) {
  stop("Type 2 diabetes, East Asian participants, European-only training does not return to excluding zero at \u03c1 = 0.9.", call. = FALSE)
}

zero_words <- function(sig03, sig06, sig09) {
  phrase <- function(s) if (isTRUE(s)) "excludes zero" else "includes zero"
  if (identical(sig03, sig06)) {
    sprintf("%s at \u03c1 = 0.3 and 0.6, and %s at \u03c1 = 0.9", phrase(sig03), phrase(sig09))
  } else {
    sprintf("%s at \u03c1 = 0.3, %s at \u03c1 = 0.6, and %s at \u03c1 = 0.9", phrase(sig03), phrase(sig06), phrase(sig09))
  }
}
list_phrase <- function(df) {
  paste(mapply(comp_phrase, df$trait_label, df$target_ancestry, df$trained_bucket), collapse = "; ")
}
still_bits <- mapply(function(tr, an, bu, s03, s06, s09) {
  sprintf("%s (%s)", comp_phrase(tr, an, bu), zero_words(s03, s06, s09))
}, still_high$trait_label, still_high$target_ancestry, still_high$trained_bucket,
still_high$r0.3, still_high$r0.6, still_high$r0.9)
east_rho <- rho_tbl %>%
  filter(
    .data$trait_label == "Type 2 diabetes",
    .data$target_ancestry == "East Asian",
    grepl("European-only", .data$trained_bucket),
    .data$rho %in% c(0, 0.9)
  )
east_0 <- east_rho$delta_hat[east_rho$rho == 0]
east_9 <- east_rho$delta_hat[east_rho$rho == 0.9]
if (length(east_0) != 1L || length(east_9) != 1L || east_0 <= 0 || east_9 >= 0) {
  stop("Type 2 diabetes, East Asian participants, European-only training does not change sign between \u03c1 = 0 and \u03c1 = 0.9.", call. = FALSE)
}
rho_caption <- paste(
  sprintf("Five comparisons change whether the interval excludes zero from \u03c1 = 0 to \u03c1 = 0.9. Three change only between 0 and 0.3: %s.", list_phrase(only_low)),
  sprintf("Two still change between 0.3 and 0.9: %s.", paste(still_bits, collapse = "; ")),
  sprintf(
    "Type 2 diabetes, East Asian participants, European-only training excludes zero at \u03c1 = 0 with estimate %s, and excludes zero again at \u03c1 = 0.9 with estimate %s. The sign reverses. Large \u03c1 makes the weights extreme as V becomes nearly singular. It is not a narrower interval around the same estimate, and it is not stronger evidence. Figure 11 shows the two comparisons whose intervals narrow as \u03c1 rises.",
    fmt_est(east_0), fmt_est(east_9)
  ),
  sep = "\n"
)

flip_keys <- bind_rows(only_low, still_high) %>%
  transmute(trait_label, trained_bucket, target_ancestry, flip = TRUE)
rho_plot_df <- rho_tbl %>%
  left_join(flip_keys, by = c("trait_label", "trained_bucket", "target_ancestry")) %>%
  mutate(
    flip = replace_na(.data$flip, FALSE),
    group = case_when(
      .data$trait_label == "Type 2 diabetes" & .data$target_ancestry == "East Asian" & grepl("European-only", .data$trained_bucket) ~ 1L,
      .data$flip & .data$rho %in% c(0.3, 0.6, 0.9) ~ 2L,
      .data$flip ~ 3L,
      TRUE ~ 4L
    ),
    label = sprintf("%s%s", if_else(.data$flip, "\u25b8 ", ""), comp_phrase(.data$trait_label, .data$target_ancestry, .data$trained_bucket)),
    rho_lab = factor(paste0("\u03c1 = ", .data$rho), levels = paste0("\u03c1 = ", c(0, 0.3, 0.6, 0.9)))
  )
label_order <- rho_plot_df %>%
  distinct(.data$label, .data$group, .data$trait_label) %>%
  arrange(.data$group, .data$trait_label) %>%
  pull(.data$label)
rho_plot_df <- rho_plot_df %>% mutate(label = factor(.data$label, levels = rev(unique(label_order))))
write_numbers(
  rho_plot_df %>%
    transmute(trait_label, trained_bucket, target_ancestry, rho, delta_hat, se_hat, lo, hi, sig, flip),
  "09_rho_grid"
)
p_rho <- ggplot(rho_plot_df, aes(x = .data$delta_hat, y = .data$label, colour = .data$rho_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(xmin = .data$lo, xmax = .data$hi),
    orientation = "y", width = 0.35, linewidth = 0.45,
    position = position_dodge(width = 0.7)
  ) +
  geom_point(size = 1.7, position = position_dodge(width = 0.7)) +
  scale_colour_manual(values = c("\u03c1 = 0" = "#B0B0B0", "\u03c1 = 0.3" = "#88CCEE", "\u03c1 = 0.6" = "#332288", "\u03c1 = 0.9" = "#D55E00")) +
  labs(
    title = "9.  The same comparisons at four values of \u03c1",
    subtitle = "A triangle in the label marks a comparison whose zero call changes. The reported intervals use \u03c1 = 0.6.",
    x = "Non-European AUC \u2212 European AUC, with 95% interval",
    y = NULL,
    caption = wrap_cap(rho_caption, width = 145)
  ) +
  theme_lesson(11) +
  theme(axis.text.y = element_text(size = 7.2))
save_lesson(p_rho, "09_rho_grid", width = 13.4, height = 13.6)

# The robust degrees of freedom are the box in figure 1.
# -----------------------------------------------------------------
# 10. What rho means
# -----------------------------------------------------------------

pair_correlation <- function(se_e, se_t, rho_e, rho_t) {
  se_d <- sqrt(se_e^2 + se_t^2)
  cov_d <- rho_e * se_e[1] * se_e[2] + rho_t * se_t[1] * se_t[2]
  c(correlation = unname(cov_d / (se_d[1] * se_d[2])), share = unname(cov_d / (rho_main * se_d[1] * se_d[2])))
}
eleven <- sa %>% filter(.data$tgt_cluster == "UK Biobank", .data$eur_cluster == "All of Us")
best <- NULL
for (i in seq_len(nrow(eleven) - 1L)) {
  for (j in (i + 1L):nrow(eleven)) {
    got <- pair_correlation(
      c(eleven$eur_se[i], eleven$eur_se[j]),
      c(eleven$tgt_se[i], eleven$tgt_se[j]),
      rho_main, rho_main
    )
    if (is.null(best) || abs(got[["correlation"]] - rho_main) < abs(best$correlation - rho_main)) {
      best <- list(
        correlation = got[["correlation"]], share = got[["share"]],
        scores = c(eleven$pgs_id[i], eleven$pgs_id[j])
      )
    }
  }
}
ea <- subset_comp(scores, "Type 2 diabetes", "Train: Multi incl. EUR", "East Asian") %>% arrange(.data$pgs_id)
if (nrow(ea) != 2L || length(unique(ea$eur_cluster)) != 1L || length(unique(ea$tgt_cluster)) != 2L) {
  stop("The East Asian comparison is not shared on the European side only.", call. = FALSE)
}
ea_got <- pair_correlation(ea$eur_se, ea$tgt_se, rho_main, 0)
if (abs(best$correlation - rho_main) > 0.005 || abs(ea_got[["correlation"]] - 0.009) > 0.0005) {
  stop("The two real pairs are not at correlations 0.60 and 0.009.", call. = FALSE)
}
if (abs(best$correlation - rho_main * best$share) > 1e-8 || abs(ea_got[["correlation"]] - rho_main * ea_got[["share"]]) > 1e-8) {
  stop("A real pair does not sit on the line.", call. = FALSE)
}

set.seed(20261102)
sim_shares <- seq(0, 1, by = 0.1)
n_rep <- 40L
n_draw <- 120L
sim_rows <- lapply(sim_shares, function(f) {
  cors <- replicate(n_rep, {
    z1 <- rnorm(n_draw)
    z2 <- rho_main * z1 + sqrt(1 - rho_main^2) * rnorm(n_draw)
    d1 <- sqrt(f) * z1 + sqrt(1 - f) * rnorm(n_draw)
    d2 <- sqrt(f) * z2 + sqrt(1 - f) * rnorm(n_draw)
    cor(d1, d2)
  })
  tibble(share = f, correlation = cors, mean_correlation = mean(cors))
})
sim_df <- bind_rows(sim_rows)
sim_gap <- sim_df %>%
  distinct(.data$share, .data$mean_correlation) %>%
  mutate(gap = abs(.data$mean_correlation - rho_main * .data$share))
if (max(sim_gap$gap) > 0.05) {
  stop(
    "The simulated correlations leave the line. Largest gap is ",
    formatC(max(sim_gap$gap), format = "f", digits = 3), ".",
    call. = FALSE
  )
}
share_line <- seq(0, 1, length.out = 201)
line_df <- tibble(share = share_line, correlation = rho_main * share_line)
real_df <- tibble(
  share = c(best$share, ea_got[["share"]]),
  correlation = c(best$correlation, ea_got[["correlation"]]),
  label = c(
    sprintf("Two of the eleven, both sides\n%s and %s, %.2f", best$scores[1], best$scores[2], best$correlation),
    sprintf("East Asian, European side only\n%.3f", ea_got[["correlation"]])
  )
)
write_numbers(
  bind_rows(
    line_df %>% mutate(kind = "line"),
    sim_df %>% mutate(kind = "simulation"),
    real_df %>% select("share", "correlation") %>% mutate(kind = "real pair")
  ),
  "10_what_rho_means"
)
p_rho_mean <- ggplot() +
  geom_line(data = line_df, aes(x = .data$share, y = .data$correlation), linewidth = 0.9, colour = "#222222") +
  geom_point(data = sim_df, aes(x = .data$share, y = .data$correlation), alpha = 0.35, size = 1.5, colour = "#332288") +
  geom_point(data = real_df, aes(x = .data$share, y = .data$correlation), size = 3.2, colour = "#D55E00") +
  geom_text(
    data = tibble(
      x = c(0.98, 0.04),
      y = c(0.82, 0.42),
      hjust = c(1, 0),
      label = real_df$label
    ),
    aes(x = .data$x, y = .data$y, hjust = .data$hjust, label = .data$label),
    size = 3.2, lineheight = 1.05, family = "DejaVu Sans"
  ) +
  scale_x_continuous(limits = c(-0.02, 1.05), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.35, 1.05), breaks = seq(-0.2, 1, by = 0.2)) +
  labs(
    title = "10.  What \u03c1 means",
    subtitle = "\u03c1 is the correlation of two scores' AUC errors on one side with the same label.",
    x = "Share of the variance on the shared side",
    y = "Correlation of the ancestry differences",
    caption = wrap_cap(paste(
      "For two scores with equal standard errors, the correlation of the ancestry differences is exactly \u03c1 times the shared side's share of the variance. The line is that rule. The blue points are sampling scatter around it, from a fixed seed.",
      sprintf("Two of the eleven share a label on both sides (%s and %s), so the correlation of their differences is %.2f.", best$scores[1], best$scores[2], best$correlation),
      "Type 2 diabetes, East Asian participants, multi-ancestry training shares on the European side only, and the correlation of the differences is 0.009, because the European standard errors are tiny. Both real pairs sit on the line.",
      "This is a sampling model. It is not a correlation estimated from these participants, and it does not change a reported interval.",
      sep = "\n"
    ), width = 118)
  ) +
  theme_lesson(12)
save_lesson(p_rho_mean, "10_what_rho_means", width = 11.2, height = 8.4)

# -----------------------------------------------------------------
# 11. Earlier checks
# -----------------------------------------------------------------

narrow <- width_narrowed(rho_tbl) %>% arrange(.data$trait_label, .data$trained_bucket, .data$target_ancestry)
expect_narrow <- tibble(
  trait_label = c("Type 2 diabetes", "Type 2 diabetes"),
  trained_bucket = c("Train: European-only", "Train: Multi incl. EUR"),
  target_ancestry = c("East Asian", "Hispanic or Latin American")
) %>% arrange(.data$trait_label, .data$trained_bucket, .data$target_ancestry)
if (!identical(
  narrow %>% select("trait_label", "trained_bucket", "target_ancestry"),
  expect_narrow
)) {
  stop("The comparisons whose intervals narrow as \u03c1 rises are not the two expected ones.", call. = FALSE)
}

min_se <- scores %>%
  group_by(.data$trait_label, .data$trained_bucket, .data$target_ancestry) %>%
  summarise(min_se = min(.data$delta_se), .groups = "drop")
earlier <- shared_cells %>%
  left_join(min_se, by = c("trait_label", "trained_bucket", "target_ancestry")) %>%
  transmute(
    trait_label, trained_bucket, target_ancestry,
    label = comp_phrase(.data$trait_label, .data$target_ancestry, .data$trained_bucket),
    narrows = paste(.data$trait_label, .data$trained_bucket, .data$target_ancestry) %in%
      paste(narrow$trait_label, narrow$trained_bucket, narrow$target_ancestry),
    avg = .data$delta_fixed, avg_lo = .data$lo_fixed, avg_hi = .data$hi_fixed,
    small = .data$delta_fixed,
    small_lo = .data$delta_fixed - 1.96 * .data$min_se,
    small_hi = .data$delta_fixed + 1.96 * .data$min_se,
    one = .data$delta_one_per_pss, one_lo = .data$lo_one_per_pss, one_hi = .data$hi_one_per_pss,
    cov = .data$delta_cov, cov_lo = .data$lo_cov, cov_hi = .data$hi_cov
  )
if (any(!is.finite(c(earlier$small_lo, earlier$one, earlier$cov)))) {
  stop("An earlier-check estimate is missing.", call. = FALSE)
}

earlier_long <- bind_rows(
  earlier %>% transmute(label, narrows, estimate = .data$avg, lo = .data$avg_lo, hi = .data$avg_hi, which = "Average of all scores"),
  earlier %>% transmute(label, narrows, estimate = .data$small, lo = .data$small_lo, hi = .data$small_hi, which = "Smallest standard error"),
  earlier %>% transmute(label, narrows, estimate = .data$one, lo = .data$one_lo, hi = .data$one_hi, which = "One score per sample set"),
  earlier %>% transmute(label, narrows, estimate = .data$cov, lo = .data$cov_lo, hi = .data$cov_hi, which = "Covariance model")
) %>%
  mutate(
    which = factor(.data$which, levels = c("Average of all scores", "Smallest standard error", "One score per sample set", "Covariance model")),
    label = factor(.data$label, levels = rev(sort(unique(.data$label))))
  )
write_numbers(earlier_long %>% mutate(label = as.character(.data$label)), "11_earlier_checks")
check_cols <- c(
  "Average of all scores" = "#222222",
  "Smallest standard error" = "#E41A1C",
  "One score per sample set" = "#6A3D9A",
  "Covariance model" = "#01665E"
)
p_earlier <- ggplot(earlier_long, aes(x = .data$estimate, y = .data$label, colour = .data$which)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(xmin = .data$lo, xmax = .data$hi),
    orientation = "y", width = 0.4, linewidth = 0.45,
    position = position_dodge(width = 0.72)
  ) +
  geom_point(size = 1.6, position = position_dodge(width = 0.72)) +
  scale_colour_manual(values = check_cols) +
  labs(
    title = "11.  Earlier checks beside the covariance model",
    subtitle = "The 24 comparisons whose scores share people.",
    x = "Non-European AUC \u2212 European AUC, with 95% interval",
    y = NULL,
    caption = wrap_cap(paste(
      "The smallest standard error keeps the plain average and uses the smallest score standard error as its standard error. That is roughly the plain average's own standard error if every error were perfectly correlated and the plain-average weights were kept. It is not the refit at \u03c1 = 1.",
      "As \u03c1 approaches 1, V becomes nearly singular. The fit can use differences between highly correlated scores to cancel error, so the interval can shrink and the weights can go extreme.",
      sprintf(
        "Two comparisons already narrow as \u03c1 rises: %s. Figure 9 points here for the first. At \u03c1 = 0.9 its estimate has changed sign and excludes zero on the other side. That is extreme weights, not stronger evidence.",
        paste(mapply(comp_phrase, narrow$trait_label, narrow$target_ancestry, narrow$trained_bucket), collapse = "; and ")
      ),
      "The sample-set term is a different model and is not drawn.",
      sep = "\n"
    ), width = 140)
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_lesson(11) +
  theme(axis.text.y = element_text(size = 7.2))
save_lesson(p_earlier, "11_earlier_checks", width = 13.2, height = 14.2)

message("Walk-through written to ", root_dir)
