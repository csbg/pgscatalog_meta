#!/usr/bin/env Rscript

# Why the weighted lm for European-only type 2 diabetes has sigma = 11.5.
# One score, PGS000805, contributes two European UK Biobank evaluations on
# the same 423,729 people. The model can fit only one value for both, so the
# 0.17 gap drives the residual scale. Stage 2 is unaffected: it pairs within
# one covariate class.
#
# Run from the repository root, after analysis/10_diabetes_cohort_models.R:
#   Rscript tests/q_weighted_lm_pgs000805.R
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
diabetes_dir <- file.path(base_dir, "diabetes_cohort_models")
table_dir <- file.path(base_dir, "tables")
root_dir <- file.path(base_dir, "q_weighted_lm_pgs000805")
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
old_files <- list.files(root_dir, pattern = "\\.(png|pdf|csv)$", full.names = TRUE)
if (length(old_files)) unlink(old_files)

read_need <- function(path) {
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run analysis/10_diabetes_cohort_models.R first.", call. = FALSE)
  }
  read_csv(path, show_col_types = FALSE)
}

eval_rows <- read_need(file.path(diabetes_dir, "diabetes_evaluation_rows.csv"))
coef_csv <- read_need(file.path(diabetes_dir, "diabetes_model_coefficients.csv"))
cells <- read_need(file.path(table_dir, "09_dependence_cells.csv"))

slide_width <- 13.333
slide_height <- 7.5
z95 <- qnorm(0.975)

fmt1 <- function(x) formatC(x, format = "f", digits = 1)
fmt3 <- function(x) formatC(x, format = "f", digits = 3)
fmt_n <- function(x) format(as.integer(round(x)), big.mark = ",")

theme_lesson <- function(base = 12) {
  theme_minimal(base_size = base, base_family = "DejaVu Sans") +
    theme(
      plot.title = element_text(face = "bold", size = base + 1.5, colour = "#1a1a1a"),
      plot.subtitle = element_text(size = base - 2, colour = "#333333", lineheight = 1.12),
      plot.caption = element_text(size = base - 3.5, colour = "#444444", hjust = 0, lineheight = 1.15),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title = element_text(size = base - 1.5, colour = "#222222"),
      axis.text = element_text(colour = "#222222"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base - 2.5),
      strip.text = element_text(face = "bold", hjust = 0, size = base - 1.5, colour = "#1a1a1a"),
      strip.background = element_rect(fill = "#F3F3F3", colour = NA),
      plot.margin = margin(6, 10, 4, 8)
    )
}

wrap_cap <- function(text, width = 150) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  paste(vapply(lines, function(line) {
    paste(strwrap(line, width = width), collapse = "\n")
  }, character(1)), collapse = "\n")
}

# Panel subtitles sit in a third of the slide; wrap them so they cannot run
# into the next panel.
wrap_sub <- function(text, width = 52) wrap_cap(text, width = width)

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

# -----------------------------------------------------------------
# European-only rows and the weighted lm Nikolaus fit
# -----------------------------------------------------------------

eur <- eval_rows %>%
  filter(.data$train_bucket == "Train: European-only") %>%
  mutate(
    ancestry = factor(
      .data$ancestry_display,
      levels = c("European", setdiff(evaluation_display_levels, "European"))
    ),
    PGS = factor(.data$pgs_id),
    weight_auc = 1 / (.data$se_auc^2)
  )
stopifnot(nrow(eur) == 56L, n_distinct(eur$pgs_id) == 13L)

fit_orig <- lm(auc ~ ancestry + PGS, data = eur, weights = weight_auc)
sm_orig <- summary(fit_orig)
sigma_orig <- sm_orig$sigma

anc_keep <- c("East Asian", "Hispanic or Latin American", "South Asian")
coef_orig <- sm_orig$coefficients
anc_rows <- grepl("^ancestry", rownames(coef_orig))
est_orig <- setNames(
  coef_orig[anc_rows, 1],
  sub("^ancestry", "", rownames(coef_orig)[anc_rows])
)
se_orig <- setNames(
  coef_orig[anc_rows, 2],
  sub("^ancestry", "", rownames(coef_orig)[anc_rows])
)

expected_csv <- coef_csv %>%
  filter(
    .data$train_bucket == "Train: European-only",
    grepl("^lm\\(AUC ~ ancestry \\+ PGS, weights", .data$model),
    .data$ancestry %in% anc_keep
  )
stopifnot(nrow(expected_csv) == 3L)
for (anc in anc_keep) {
  row <- expected_csv %>% filter(.data$ancestry == .env$anc)
  if (abs(est_orig[[anc]] - row$estimate) > 1e-8 || abs(se_orig[[anc]] - row$se) > 1e-8) {
    stop(
      "Refit of the weighted lm does not match diabetes_model_coefficients.csv for ",
      anc, ".",
      call. = FALSE
    )
  }
}
# The three numbers called out in the note for Georg.
if (
  abs(est_orig[["East Asian"]] - 0.011883751017775643) > 1e-8 ||
    abs(est_orig[["Hispanic or Latin American"]] - 0.04531639285060013) > 1e-8 ||
    abs(est_orig[["South Asian"]] - (-0.0033809643158825326)) > 1e-8
) {
  stop("Weighted lm ancestry coefficients drifted from 0.0119, 0.0453, -0.0034.", call. = FALSE)
}

dup_ppm <- c("PPM002120", "PPM002122")
dup_rows <- eur %>% filter(.data$performance_id %in% dup_ppm)
stopifnot(nrow(dup_rows) == 2L)
stopifnot(all(dup_rows$pgs_id == "PGS000805"))
stopifnot(all(dup_rows$sampleset_id == "PSS001046"))
stopifnot(n_distinct(dup_rows$n_individuals) == 1L)
n_dup_people <- dup_rows$n_individuals[[1]]

# Residual in units of the reported SE. If SE^2 were the true variance these
# would be near 1; the two UK Biobank rows are near ±50.
w_fit <- weights(fit_orig)
e_fit <- residuals(fit_orig)
resid_se <- e_fit / eur$se_auc
wss <- w_fit * e_fit^2
dup_idx <- which(eur$performance_id %in% dup_ppm)
rss_share <- sum(wss[dup_idx]) / sum(wss)
if (!is.finite(rss_share) || abs(rss_share - sum(wss[dup_idx]) / sum(wss)) > 0) {
  stop("Residual sum of squares share was not recomputed from the fit.", call. = FALSE)
}
stopifnot(rss_share > 0.88, rss_share < 0.92)
stopifnot(abs(sigma_orig - 11.5) < 0.1)

fitted_eur_ukb <- unique(fitted(fit_orig)[dup_idx])
stopifnot(length(fitted_eur_ukb) == 1L)

# Drop the age/sex/BMI/PCs European UK Biobank row; keep score only, as the
# Stage 1 class priority does when more than one class is present.
eur_drop <- eur %>% filter(.data$performance_id != "PPM002120")
fit_drop <- lm(auc ~ ancestry + PGS, data = eur_drop, weights = weight_auc)
sm_drop <- summary(fit_drop)
sigma_drop <- sm_drop$sigma
coef_drop <- sm_drop$coefficients
anc_drop <- grepl("^ancestry", rownames(coef_drop))
est_drop <- setNames(
  coef_drop[anc_drop, 1],
  sub("^ancestry", "", rownames(coef_drop)[anc_drop])
)
se_drop <- setNames(
  coef_drop[anc_drop, 2],
  sub("^ancestry", "", rownames(coef_drop)[anc_drop])
)

# Same point estimates as the original weighted lm; SE as if sigma = 1.
est_sigma1 <- est_orig
se_sigma1 <- se_orig / sigma_orig

cov_cells <- cells %>%
  filter(
    .data$trait_label == "Type 2 diabetes",
    .data$trained_bucket == "Train: European-only",
    .data$target_ancestry %in% anc_keep
  )
stopifnot(nrow(cov_cells) == 3L)

# Multi-ancestry footnote: PGS000804 has the same two-model pattern on three
# sample sets.
multi_dup <- eval_rows %>%
  filter(
    .data$train_bucket == "Train: Multi incl. EUR",
    .data$pgs_id == "PGS000804"
  ) %>%
  count(.data$ancestry_display, .data$sampleset_id, name = "n") %>%
  filter(.data$n == 2L)
stopifnot(nrow(multi_dup) == 3L)
stopifnot(setequal(
  multi_dup$ancestry_display,
  c("African", "European", "Hispanic or Latin American")
))

# -----------------------------------------------------------------
# Panel A: what the model sees for PGS000805
# -----------------------------------------------------------------

class_short <- c(
  "None (PRS only)" = "Score only",
  "Age/sex/PCs/technical" = "Age/sex/PCs",
  "Clinical or other covariates" = "Age/sex/BMI/PCs"
)
class_cols <- c(
  "Score only" = "#332288",
  "Age/sex/PCs" = "#E69F00",
  "Age/sex/BMI/PCs" = "#D55E00"
)

# Put the two UK Biobank European rows adjacent for the bracket.
order805 <- c(
  "PPM021603", # East Asian UKB
  "PPM021629", # South Asian UKB
  "PPM021497", # Hispanic All of Us
  "PPM021422", # European All of Us
  "PPM002120", # European UKB clinical
  "PPM002122"  # European UKB score only
)

pgs805 <- eur %>%
  filter(.data$pgs_id == "PGS000805") %>%
  mutate(
    class_label = unname(class_short[.data$covariate_class]),
    row_label = paste0(
      .data$ancestry_display, ", ", .data$study_name, "\n", .data$class_label
    ),
    is_dup = .data$performance_id %in% dup_ppm,
    performance_id = factor(.data$performance_id, levels = rev(order805))
  ) %>%
  arrange(.data$performance_id)
stopifnot(nrow(pgs805) == 6L)
stopifnot(!anyDuplicated(as.character(pgs805$row_label)))
pgs805 <- pgs805 %>%
  mutate(row_label = factor(.data$row_label, levels = .data$row_label))

dup_y <- as.numeric(pgs805$performance_id[match(dup_ppm, pgs805$performance_id)])
dup_y_mid <- mean(dup_y)
auc_clinical <- pgs805$auc[pgs805$performance_id == "PPM002120"]
auc_score_only <- pgs805$auc[pgs805$performance_id == "PPM002122"]
bracket_lab <- paste0(
  "The same ", fmt_n(n_dup_people), " people.\n",
  "Score with age, sex, BMI, PCs: ", fmt3(auc_clinical), "\n",
  "Score alone: ", fmt3(auc_score_only)
)

p_a <- ggplot(pgs805, aes(y = .data$performance_id)) +
  geom_errorbar(
    aes(xmin = .data$estimate_ci_lower, xmax = .data$estimate_ci_upper, colour = .data$class_label),
    orientation = "y", width = 0.2, linewidth = 0.55
  ) +
  geom_point(aes(x = .data$auc, colour = .data$class_label), size = 2.4) +
  # Fitted value tick between the two UK Biobank European rows.
  annotate(
    "segment",
    x = fitted_eur_ukb, xend = fitted_eur_ukb,
    y = min(dup_y) - 0.28, yend = max(dup_y) + 0.28,
    colour = "#222222", linewidth = 0.9
  ) +
  annotate(
    "point",
    x = fitted_eur_ukb, y = dup_y_mid,
    shape = 18, size = 3.2, colour = "#222222"
  ) +
  annotate(
    "text",
    x = fitted_eur_ukb + 0.012, y = dup_y_mid,
    label = paste0("Fitted\n", fmt3(fitted_eur_ukb)),
    size = 2.6, family = "DejaVu Sans", colour = "#222222", hjust = 0, lineheight = 0.95
  ) +
  # Bracket to the right of the two UK Biobank European rows.
  annotate(
    "segment",
    x = 0.855, xend = 0.855,
    y = min(dup_y), yend = max(dup_y),
    colour = "#444444", linewidth = 0.55
  ) +
  annotate(
    "segment",
    x = 0.84, xend = 0.855,
    y = min(dup_y), yend = min(dup_y),
    colour = "#444444", linewidth = 0.55
  ) +
  annotate(
    "segment",
    x = 0.84, xend = 0.855,
    y = max(dup_y), yend = max(dup_y),
    colour = "#444444", linewidth = 0.55
  ) +
  annotate(
    "text",
    x = 0.865, y = dup_y_mid,
    label = bracket_lab,
    hjust = 0, size = 2.45, family = "DejaVu Sans", colour = "#333333", lineheight = 1.05
  ) +
  scale_y_discrete(
    breaks = levels(pgs805$performance_id),
    labels = levels(pgs805$row_label)
  ) +
  scale_colour_manual(values = class_cols, breaks = names(class_cols)) +
  coord_cartesian(xlim = c(0.62, 0.98), clip = "off") +
  scale_x_continuous(breaks = seq(0.65, 0.85, by = 0.1)) +
  labs(
    title = "A. The six evaluations of PGS000805",
    subtitle = wrap_sub(paste0(
      "Two are Europeans in UK Biobank, once with BMI in the model and once without. ",
      "The weighted lm can fit only one value for both."
    )),
    x = "AUC with 95% interval",
    y = NULL
  ) +
  theme_lesson(11) +
  theme(
    axis.text.y = element_text(size = 8, lineheight = 0.95),
    legend.position = "bottom",
    plot.margin = margin(6, 64, 4, 8),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3)
  )

# -----------------------------------------------------------------
# Panel B: where sigma = 11.5 comes from
# -----------------------------------------------------------------

resid_df <- eur %>%
  mutate(
    resid_se = as.numeric(resid_se),
    is_dup = .data$performance_id %in% dup_ppm,
    label = if_else(
      .data$is_dup,
      paste0("European, UK Biobank\n", class_short[.data$covariate_class]),
      NA_character_
    )
  ) %>%
  arrange(.data$resid_se) %>%
  mutate(rank = row_number())

p_b <- ggplot(resid_df, aes(x = .data$rank, y = .data$resid_se)) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "grey55", linewidth = 0.45) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(
    data = resid_df %>% filter(!.data$is_dup),
    colour = "#B0B0B0", size = 1.7
  ) +
  geom_point(
    data = resid_df %>% filter(.data$is_dup),
    aes(colour = unname(class_short[.data$covariate_class])),
    size = 3.1
  ) +
  ggrepel::geom_text_repel(
    data = resid_df %>% filter(.data$is_dup),
    aes(label = .data$label, colour = unname(class_short[.data$covariate_class])),
    size = 2.5, family = "DejaVu Sans", lineheight = 0.95,
    min.segment.length = 0, box.padding = 0.35, show.legend = FALSE
  ) +
  scale_colour_manual(values = class_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.10))) +
  annotate(
    "label",
    x = 28, y = 38,
    label = paste0(
      "These 2 of ", nrow(eur), " evaluations give\n",
      formatC(100 * rss_share, format = "f", digits = 0),
      "% of the residual sum of squares."
    ),
    size = 2.8, family = "DejaVu Sans", hjust = 0.5, vjust = 1,
    fill = alpha("white", 0.85), colour = "#222222", lineheight = 1.05
  ) +
  labs(
    title = paste0("B. Where \u03c3 = ", fmt1(sigma_orig), " comes from"),
    subtitle = wrap_sub(paste0(
      "Distance of each evaluation from its fitted value, in standard errors. ",
      "If the standard errors were right, most would fall between the dashed lines."
    )),
    x = paste0(nrow(eur), " evaluations, European-only training, sorted"),
    y = "(AUC \u2212 fitted) / SE"
  ) +
  theme_lesson(11) +
  theme(panel.grid.major.x = element_blank())

# -----------------------------------------------------------------
# Panel C: ancestry estimates under four fits
# -----------------------------------------------------------------

# Look the values up before tibble(): inside tibble(), a later argument that
# says `estimate` sees the new one-value column, not the named vector.
ci_tbl <- function(estimate, se, ancestry, fit_label, sigma_val) {
  est <- unname(estimate[ancestry])
  s <- unname(se[ancestry])
  stopifnot(is.finite(est), is.finite(s))
  tibble(
    ancestry = ancestry,
    fit = fit_label,
    estimate = est,
    se = s,
    lo = est - z95 * s,
    hi = est + z95 * s,
    sigma = sigma_val
  )
}

fit_orig_lab <- paste0("Weighted lm as fitted (\u03c3 = ", fmt1(sigma_orig), ")")
fit_drop_lab <- paste0("Without the BMI evaluation (\u03c3 = ", fmt1(sigma_drop), ")")
fit_sig1_lab <- "Standard errors taken as known (\u03c3 = 1)"
fit_cov_lab <- "Covariance model (\u03c1 = 0.6)"
fit_levels <- c(fit_orig_lab, fit_drop_lab, fit_sig1_lab, fit_cov_lab)
fit_cols <- setNames(c("#7A7A7A", "#332288", "#009E73", "#8C2D04"), fit_levels)

panel_c <- bind_rows(
  lapply(anc_keep, function(a) ci_tbl(est_orig, se_orig, a, fit_orig_lab, sigma_orig)),
  lapply(anc_keep, function(a) ci_tbl(est_drop, se_drop, a, fit_drop_lab, sigma_drop)),
  lapply(anc_keep, function(a) ci_tbl(est_sigma1, se_sigma1, a, fit_sig1_lab, 1)),
  cov_cells %>%
    transmute(
      ancestry = .data$target_ancestry,
      fit = fit_cov_lab,
      estimate = .data$delta_cov,
      se = .data$se_cov,
      lo = .data$lo_cov,
      hi = .data$hi_cov,
      sigma = NA_real_
    )
) %>%
  mutate(
    ancestry = factor(.data$ancestry, levels = anc_keep),
    # First listed fit at the top of each facet.
    fit = factor(.data$fit, levels = rev(fit_levels))
  )

# Drawn covariance values must match the comparison table.
for (anc in anc_keep) {
  drawn <- panel_c %>% filter(.data$ancestry == .env$anc, .data$fit == fit_cov_lab)
  cell <- cov_cells %>% filter(.data$target_ancestry == .env$anc)
  if (
    abs(drawn$estimate - cell$delta_cov) > 1e-8 ||
      abs(drawn$lo - cell$lo_cov) > 1e-8 ||
      abs(drawn$hi - cell$hi_cov) > 1e-8
  ) {
    stop("Drawn covariance values do not match 09_dependence_cells.csv for ", anc, ".", call. = FALSE)
  }
}

stopifnot(all(is.finite(panel_c$lo)), all(is.finite(panel_c$hi)))

p_c <- ggplot(panel_c, aes(y = .data$fit, colour = .data$fit)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(xmin = .data$lo, xmax = .data$hi),
    orientation = "y", width = 0.2, linewidth = 0.55
  ) +
  geom_point(aes(x = .data$estimate), size = 2.2) +
  facet_wrap(~ancestry, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = fit_cols, breaks = fit_levels) +
  scale_x_continuous(breaks = seq(-0.2, 0.2, by = 0.1)) +
  coord_cartesian(xlim = c(-0.22, 0.22)) +
  labs(
    title = "C. The ancestry differences",
    subtitle = wrap_sub(paste0(
      "Coefficient and 95% interval from each weighted lm, with the covariance ",
      "model for reference. Grey and green share the point; only \u03c3 differs."
    )),
    x = "\u0394AUC (non-European \u2212 European)",
    y = NULL
  ) +
  theme_lesson(11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.box.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 9),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3)
  ) +
  guides(colour = guide_legend(ncol = 1))

# -----------------------------------------------------------------
# Compose, write numbers, save
# -----------------------------------------------------------------

caption <- wrap_cap(paste0(
  "Type 2 diabetes, European-only training. In lm, the weights 1/SE\u00b2 only say how much each evaluation counts ",
  "relative to the others. lm then estimates its own error scale \u03c3 and multiplies every standard error by it. ",
  "Removing the BMI evaluation lowers \u03c3 from ", fmt1(sigma_orig), " to ", fmt1(sigma_drop), ". ",
  "Even then, the weighted lm treats evaluations on the same people as independent: this figure is about the weights, not about shared samples.\n",
  "PGS000804 (multi-ancestry training) has the same pattern on African, European, and Hispanic or Latin American sample sets. ",
  "The main \u0394AUCs are not affected, because Stage 2 pairs European and non-European evaluations only within one covariate set."
), width = 175)

fig <- (p_a | p_b | p_c) +
  plot_layout(widths = c(1.1, 0.9, 1.1)) +
  plot_annotation(
    title = paste0(
      "Why the weighted lm has \u03c3 = ", fmt1(sigma_orig),
      ": one score evaluated twice on the same people"
    ),
    caption = caption,
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 15, family = "DejaVu Sans", colour = "#1a1a1a"
      ),
      plot.caption = element_text(
        size = 8.5, family = "DejaVu Sans", colour = "#444444", hjust = 0, lineheight = 1.15
      ),
      plot.margin = margin(8, 12, 6, 8)
    )
  )

save_slide(fig, "01_weighted_lm_pgs000805")

numbers <- bind_rows(
  pgs805 %>%
    transmute(
      panel = "A",
      performance_id, pgs_id, sampleset_id, ancestry_display, study_name,
      covariate_class, class_label, auc,
      lo = .data$estimate_ci_lower, hi = .data$estimate_ci_upper,
      fitted_eur_ukb = ifelse(.data$is_dup, fitted_eur_ukb, NA_real_),
      resid_se = NA_real_, rss_share = NA_real_, sigma = NA_real_,
      fit = NA_character_, estimate = NA_real_, se = NA_real_
    ),
  resid_df %>%
    transmute(
      panel = "B",
      performance_id, pgs_id, sampleset_id, ancestry_display, study_name,
      covariate_class, class_label = unname(class_short[.data$covariate_class]),
      auc, lo = NA_real_, hi = NA_real_,
      fitted_eur_ukb = NA_real_,
      resid_se, rss_share = .env$rss_share, sigma = .env$sigma_orig,
      fit = NA_character_, estimate = NA_real_, se = NA_real_
    ),
  panel_c %>%
    transmute(
      panel = "C",
      performance_id = NA_character_, pgs_id = NA_character_,
      sampleset_id = NA_character_, ancestry_display = as.character(.data$ancestry),
      study_name = NA_character_, covariate_class = NA_character_,
      class_label = NA_character_, auc = NA_real_,
      lo, hi, fitted_eur_ukb = NA_real_, resid_se = NA_real_,
      rss_share = NA_real_, sigma, fit = as.character(.data$fit),
      estimate, se
    )
)

write_csv(numbers, file.path(root_dir, "01_weighted_lm_pgs000805.csv"))

message(glue::glue(
  "Weighted-lm checks passed: coefficients match CSV; ",
  "2 of 56 rows = {formatC(100 * rss_share, format = 'f', digits = 1)}% of residual sum of squares; ",
  "sigma = {fmt1(sigma_orig)}."
))
