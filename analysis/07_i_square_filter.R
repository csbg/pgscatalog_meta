#!/usr/bin/env Rscript

# ============================================================
# Script: 07_i_square_filter.R
# Project: pgscatalog_meta
# Purpose: Apply heterogeneity filtering to Stage-1 pooled AUC estimates using an I² threshold.
# Inputs: meta_roadmap_two_stages/<DATE>/tables/stage1_pooled_cells_<DATE>.csv
# Outputs: meta_roadmap_two_stages/<DATE>/tables/stage1_pooled_cells_I2filtered_<DATE>.csv; I² summary tables and plots
# Run after: 06_meta_ivw.R
# Run before: 08_generate_publication_figures.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr); library(glue); library(stringr)
})

if (!file.exists("R/load.R")) {
  stop("R/load.R not found. Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

# ---------------------------
# 0) CONFIG
# ---------------------------
stamp <- pipeline_stamp
root_pw    <- file.path("meta_roadmap_two_stages", stamp)
tables_dir <- file.path(root_pw, "tables")
plots_dir  <- file.path(root_pw, "plots", "heterogeneity")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

slugify <- function(x) stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")

# --- Ancestry color palette (consistent with other scripts) ---
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

# ---------------------------
# 1) LOAD STAGE-1 TABLE
# ---------------------------
fp_stage1 <- file.path(tables_dir, paste0("stage1_pooled_cells_", stamp, ".csv"))
if (!file.exists(fp_stage1)) {
  stop(glue("Stage-1 file not found: {fp_stage1}. Please run 06_meta_ivw.R first."))
}

message(glue(">> Loading Stage-1 data from: {fp_stage1}"))
df <- read_csv(fp_stage1, show_col_types = FALSE)

required_cols <- c("trait_label", "pgs_id", "ancestry_display", "I2", "k_eval")
missing_cols  <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(glue("Stage-1 file is missing required columns: {paste(missing_cols, collapse = ', ')}"))
}

# ---------------------------
# 2) SUMMARY OF I² PER CELL
# ---------------------------
df_i2 <- df %>%
  filter(!is.na(I2), k_eval >= 2) %>%
  mutate(
    trait_label      = as.character(trait_label),
    ancestry_display = as.character(ancestry_display)
  )

if (nrow(df_i2) == 0) {
  stop("No replicated cells with non-NA I² (k_eval >= 2). Cannot perform heterogeneity filter.")
}

i2_summary <- df_i2 %>%
  arrange(trait_label, ancestry_display, pgs_id) %>%
  select(trait_label, pgs_id, ancestry_display, I2, k_eval, dplyr::any_of(c("Q", "auc_pooled")))

out_summary <- file.path(tables_dir, paste0("stage1_I2_summary_", stamp, ".csv"))
write_csv(i2_summary, out_summary)
message(glue(">> Wrote I² summary: {out_summary}"))

# ---------------------------
# 3) FILTER OUT HIGH I² CELLS (> 80%)
# ---------------------------
df_flagged <- df %>%
  mutate(
    flag_high_I2 = flag_high_i2(I2, k_eval)
  )

n_high <- sum(df_flagged$flag_high_I2, na.rm = TRUE)
n_tot  <- nrow(df_flagged)
message(glue(">> Flagging {n_high} / {n_tot} rows with I² > 80 and k_eval >= 2."))

stage1_filtered <- df_flagged %>%
  filter(!flag_high_I2) %>%
  select(-flag_high_I2)

out_filtered <- file.path(tables_dir, paste0("stage1_pooled_cells_I2filtered_", stamp, ".csv"))
write_csv(stage1_filtered, out_filtered)
message(glue(">> Wrote filtered Stage-1 table: {out_filtered}"))

# ---------------------------
# 4) BAR PLOTS OF I² PER TRAIT
# ---------------------------
df_i2_plot <- df_i2 %>%
  mutate(
    ancestry_display = factor(
      ancestry_display,
      levels = intersect(names(pal_ancestry), unique(ancestry_display)),
      ordered = TRUE
    )
  )

traits <- sort(unique(df_i2_plot$trait_label))
pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf

for (tr in traits) {
  d_tr <- df_i2_plot %>% filter(trait_label == tr)
  if (nrow(d_tr) == 0) next

  # Order PGS by European pooled AUC if present, otherwise by median I²
  if ("auc_pooled" %in% names(d_tr)) {
    # First compute a global fallback from the full data
    global_fallback <- median(d_tr$auc_pooled, na.rm = TRUE)
    if (is.na(global_fallback)) {
      global_fallback <- 0
    }
    
    pgs_order_tbl <- d_tr %>%
      group_by(pgs_id) %>%
      summarise(
        m = median(auc_pooled[ancestry_display == "European"], na.rm = TRUE),
        .groups = "drop"
      )
    
    pgs_order <- pgs_order_tbl %>%
      mutate(m = ifelse(is.na(m), global_fallback, m)) %>%
      arrange(m) %>%
      pull(pgs_id)
  } else {
    pgs_order <- d_tr %>%
      group_by(pgs_id) %>%
      summarise(m = median(I2, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(m)) %>%
      pull(pgs_id)
  }

  d_tr <- d_tr %>%
    mutate(pgs_id = factor(pgs_id, levels = pgs_order))

  g <- ggplot(d_tr, aes(x = pgs_id, y = I2, fill = ancestry_display)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 80, linetype = "dashed", linewidth = 0.5, color = "grey40") +
    scale_fill_manual(values = pal_ancestry, name = "Ancestry") +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title    = glue("Heterogeneity (I²) by PGS and ancestry — {tr}"),
      subtitle = "Dashed line at 80%: evaluations above this threshold are excluded from downstream analyses.",
      x        = "PGS ID",
      y        = "I² (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank()
    )

  # Dynamic width/height
  n_pgs <- length(unique(d_tr$pgs_id))
  w <- max(7, 3 + 0.25 * n_pgs)
  h <- max(4, 3 + 0.25 * n_pgs)

  out_pdf <- file.path(plots_dir, paste0("I2_bar_", slugify(tr), "_", stamp, ".pdf"))
  ggsave(out_pdf, g, width = w, height = h, device = pdf_device, limitsize = FALSE)
  message(glue("   - Saved I² bar plot for {tr}: {out_pdf}"))
}

message("\n✅ Heterogeneity filter (I² > 80%) complete.")
