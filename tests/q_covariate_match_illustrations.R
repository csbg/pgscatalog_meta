#!/usr/bin/env Rscript

# ============================================================
# Covariate class on both sides of a delta.
# A delta is one score and one covariate class. A comparison has one
# average. Its colour is score-only, covariate-adjusted, or both kinds.
# Five scores change value relative to the pairing that ignored covariate class.
# Seven more scores enter because each class is now its own side.
#
# Run from the repository root, after analysis/09_review_sensitivity.R:
#   Rscript tests/q_covariate_match_illustrations.R
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
root_dir <- file.path("results", paste0("review_sensitivity_", stamp), "q_covariate_match_illustrations")
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

pair_path <- file.path(table_dir, "10_covariates_pairs.csv")
colour_path <- file.path(table_dir, "10_covariates_comparison_colour.csv")
change_path <- file.path(table_dir, "10_covariates_rule_changes.csv")
stopifnot(file.exists(pair_path), file.exists(colour_path), file.exists(change_path))

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

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

pairs <- read_csv(pair_path, show_col_types = FALSE)
colour <- read_csv(colour_path, show_col_types = FALSE)
changes <- read_csv(change_path, show_col_types = FALSE)

stopifnot(nrow(pairs) == 127L, nrow(colour) == 32L)
stopifnot(!anyDuplicated(pairs[c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")]))
stopifnot(all(pairs$covariate_class %in% covariate_class_levels))
stopifnot(all(colour$covariate_colour %in% covariate_colour_levels))
stopifnot(!any(colour$covariate_colour == "Both, or mixed within a score"))

chd_trait <- "Coronary heart disease (incident and prevalent)"
chd_bucket <- "Train: Multi incl. EUR"
chd_colour <- colour %>%
  filter(trait_label == chd_trait, trained_bucket == chd_bucket)
stopifnot(nrow(chd_colour) == 2L)
stopifnot(all(chd_colour$covariate_colour == "Uses covariates"))

expected <- tibble(
  trait_label = c(rep(chd_trait, 4), "Type 2 diabetes"),
  trained_bucket = chd_bucket,
  target_ancestry = c(
    "African", "Hispanic or Latin American",
    "African", "Hispanic or Latin American",
    "Middle Eastern or North African"
  ),
  pgs_id = c("PGS000013", "PGS000013", "PGS000018", "PGS000018", "PGS000804")
)
value_changes <- changes %>% filter(change_type == "value")
stopifnot(nrow(inner_join(
  value_changes, expected,
  by = c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")
)) == 5L)
stopifnot(nrow(value_changes) == 5L)
stopifnot(sum(changes$change_type == "admitted") == 7L)
stopifnot(!any(changes$change_type == "dropped"))

changes <- value_changes %>%
  mutate(
    label = paste0(pgs_id, "\n", target_ancestry),
    label = factor(label, levels = rev(label))
  )

long <- bind_rows(
  changes %>% transmute(label, delta = delta_pre_rule, rule = "Before the class rule"),
  changes %>% transmute(label, delta = delta_new, rule = "One class on both sides")
) %>%
  mutate(rule = factor(rule, levels = c("Before the class rule", "One class on both sides")))

p <- ggplot(long, aes(x = delta, y = label, colour = rule)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = 2.8, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("Before the class rule" = "#B0B0B0", "One class on both sides" = "#332288")) +
  labs(
    title = "Five scores whose delta changes",
    subtitle = "Each score now uses one covariate class on the European side and the non-European side.",
    x = "ΔAUC (non-European − European)",
    y = NULL,
    caption = paste0(
      "Coronary heart disease, PGS000013 and PGS000018: the European side is the clinical evaluation. ",
      "Type 2 diabetes, PGS000804, Middle Eastern or North African participants: the clinical European evaluation ",
      "is paired with the clinical evaluation, and the delta is ",
      fmt(changes$delta_new[changes$pgs_id == "PGS000804"]),
      "."
    )
  ) +
  theme_lesson(12)
save_lesson(p, "01_five_changed_deltas", width = 9.2, height = 5.2)

write_csv(
  changes %>%
    select(
      trait_label, trained_bucket, target_ancestry, pgs_id, covariate_class,
      delta_pre_rule, delta_new, delta_diff
    ),
  file.path(root_dir, "five_changed_deltas.csv")
)

message("Covariate-class checks passed: 127 pairs, 32 comparisons, five changed scores, seven newly paired scores.")
