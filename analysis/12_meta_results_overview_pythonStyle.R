#!/usr/bin/env Rscript
## ============================================================
## 12_results_plots_extra_forest_traitY_faceted.R
##
## Complementary ΔAUC forest:
##   - Uses the SAME data prep as Plot 3 in 12_results_plots.R
##   - Uses a faceted layout with facet_grid to ensure proportional
##     panel heights (space = "free_y") and reduce empty spaces.
##   - Hides y-axis labels (ancestries) to reduce clutter.
##   - Trait names are displayed as left-aligned facet strips
##     acting as the new y-axis.
##
## Does NOT modify or replace the original three plots.
## Run in the SAME R session after 12_results_plots.R.
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(forcats)
})

message(">> Extra: Faceted ΔAUC forest (proportional panel heights, clean y-axis)")

## ------------------------------------------------------------------
## 0) Validate environment & paths
## ------------------------------------------------------------------

needed <- c("stage1_tbl", "keep_bucket_pattern", "keep_traits",
            "display_levels", "pal_ancestry", "wrap_trait",
            "plots_pw", "pdf_device", "stamp")
missing <- needed[!vapply(needed, exists, logical(1))]
if (length(missing) > 0) {
  stop(glue(
    "Missing objects: {paste(missing, collapse = ', ')}.\n",
    "Run 12_results_plots.R in this session before this script."
  ))
}

root_results <- dirname(plots_pw)
results_dir  <- file.path(root_results, "plots", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
message(glue(">> Results plots will be saved to: {results_dir}"))

## ------------------------------------------------------------------
## 1) REBUILD paired_df EXACTLY as in 12_results_plots.R
## ------------------------------------------------------------------

scat_df <- stage1_tbl %>%
  transmute(
    trait_label,
    ancestry_display = as.character(ancestry_display),
    pgs_id,
    auc_pooled,
    se,
    k_eval,
    trained_bucket
  ) %>%
  filter(
    ancestry_display != "Not reported",
    grepl(paste(keep_bucket_pattern, collapse = "|"),
          trained_bucket, ignore.case = TRUE)
  ) %>%
  mutate(trait_wrapped_check = wrap_trait(trait_label)) %>%
  filter(trait_wrapped_check %in% keep_traits) %>%
  select(-trait_wrapped_check)

eur_df <- scat_df %>%
  filter(ancestry_display == "European") %>%
  transmute(
    trait_label, pgs_id, trained_bucket,
    eur_auc = auc_pooled,
    eur_se  = se,
    eur_k   = k_eval,
    eur_lo  = auc_pooled - 1.96 * se,
    eur_hi  = auc_pooled + 1.96 * se
  )

tgt_df <- scat_df %>%
  filter(ancestry_display != "European") %>%
  transmute(
    trait_label, pgs_id, trained_bucket,
    target_ancestry = ancestry_display,
    tgt_auc = auc_pooled,
    tgt_se  = se,
    tgt_k   = k_eval,
    tgt_lo  = auc_pooled - 1.96 * se,
    tgt_hi  = auc_pooled + 1.96 * se
  )

paired_df <- tgt_df %>%
  inner_join(eur_df, by = c("trait_label", "pgs_id", "trained_bucket")) %>%
  mutate(
    is_pooled = (eur_k >= 2) & (tgt_k >= 2),
    target_ancestry = factor(
      target_ancestry,
      levels = intersect(display_levels, unique(target_ancestry)),
      ordered = TRUE
    )
  ) %>%
  filter(is.finite(eur_auc), is.finite(tgt_auc))

message(glue(">> Extra forest: {nrow(paired_df)} paired dots across {n_distinct(paired_df$trait_label)} traits"))

if (nrow(paired_df) == 0) {
  stop("No paired EUR/target data — cannot build ΔAUC forest.")
}

## ------------------------------------------------------------------
## 2) REBUILD delta_df EXACTLY as in Plot 3 (no extra filters)
## ------------------------------------------------------------------

delta_df <- paired_df %>%
  mutate(
    delta           = tgt_auc - eur_auc,
    delta_se        = sqrt(tgt_se^2 + eur_se^2),
    flag_not_pooled = !is_pooled
  ) %>%
  filter(is.finite(delta), is.finite(delta_se), delta_se > 0) %>%
  group_by(trait_label, trained_bucket, target_ancestry) %>%
  summarise(
    n_pairs        = n(),
    w              = 1 / (delta_se^2),
    delta_hat      = sum(w * delta) / sum(w),
    se_hat         = sqrt(1 / sum(w)),
    all_not_pooled = all(flag_not_pooled, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  mutate(
    lo             = delta_hat - 1.96 * se_hat,
    hi             = delta_hat + 1.96 * se_hat,
    pooled_group   = all_not_pooled,
    target_ancestry = factor(
      target_ancestry,
      levels = intersect(display_levels, unique(target_ancestry)),
      ordered = TRUE
    ),
    trained_bucket = factor(
      trained_bucket,
      levels = sort(unique(as.character(trained_bucket)))
    ),
    # Pre-wrap traits for cleaner facet headers
    trait_clean = wrap_trait(trait_label)
  )

if (nrow(delta_df) == 0) {
  stop("No valid ΔAUC data after pooling — cannot build extra forest.")
}

## ------------------------------------------------------------------
## 3) Axis limits and colors
## ------------------------------------------------------------------

rng     <- range(c(delta_df$lo, delta_df$hi), na.rm = TRUE)
max_abs <- max(abs(rng))
lim_x   <- ceiling(max_abs * 20) / 20

anc_targets <- levels(delta_df$target_ancestry)
anc_colors  <- pal_ancestry[anc_targets]

## ------------------------------------------------------------------
## 4) Faceted ΔAUC forest
## ------------------------------------------------------------------

# Use standard dodge to separate training buckets within the same ancestry
pos_dodge <- position_dodge(width = 0.5)

g_delta_faceted <- ggplot(
  delta_df,
  aes(
    x        = delta_hat,
    y        = target_ancestry,
    color    = target_ancestry,
    shape    = pooled_group,
    linetype = trained_bucket
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_errorbarh(
    aes(xmin = lo, xmax = hi, group = interaction(target_ancestry, trained_bucket)),
    linewidth = 0.3,
    height    = 0,
    alpha     = 0.9,
    position  = pos_dodge
  ) +
  geom_point(
    aes(group = interaction(target_ancestry, trained_bucket)),
    size     = 1.4,
    stroke   = 0.3,
    position = pos_dodge
  ) +
  scale_color_manual(
    values = anc_colors,
    name   = "Target ancestry"
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 1, `FALSE` = 16),
    labels = c(
      `TRUE`  = "All single-eval (no pooled cells)",
      `FALSE` = "Includes pooled cells"
    ),
    name   = "Pooling status"
  ) +
  scale_linetype_discrete(name = "Training weights") +
  scale_x_continuous(
    limits = c(-lim_x, lim_x),
    breaks = pretty(c(-lim_x, lim_x), n = 5),
    name   = expression(Delta * AUC ~ "(target - EUR)")
  ) +
  # This section keeps your traits on the y-axis
  facet_grid(
    trait_clean ~ ., 
    scales = "free_y", 
    space = "free_y",
    switch = "y" # Moves the trait labels to the left side
  ) +
  labs(
    y = NULL, 
    caption = "Points = pooled ΔAUC across PGS; hollow = all single-eval (no pooled cells); filled = includes pooled cells; 95% CI from inverse-variance-weighted ΔAUC."
  ) +
  theme_CrossAncestryGenPhen(
    base_size    = 6,
    show_borders = TRUE,
    show_grid    = FALSE
  ) +
  theme(
    legend.position   = "right",
    legend.box        = "vertical",
    plot.caption      = element_text(size = 6, hjust = 0),
    panel.spacing     = unit(0.2, "lines"),
    
    # This removes the repeating ancestries from the y-axis
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    
    # This styles your traits so they look like normal y-axis labels
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6, margin = margin(r = 5)) 
  )

## ------------------------------------------------------------------
## 5) Save the faceted plot
## ------------------------------------------------------------------

# Calculate height dynamically based on the number of unique y-axis groups
n_traits     <- n_distinct(delta_df$trait_clean)
n_total_rows <- n_distinct(paste(delta_df$trait_clean, delta_df$target_ancestry))
h_delta      <- max(4.0, (n_total_rows * 0.15) + (n_traits * 0.2) + 1.5)

out_delta_faceted <- file.path(
  results_dir,
  paste0("deltaAUC_forest_faceted_", stamp, ".pdf")
)

ggsave(
  out_delta_faceted,
  g_delta_faceted,
  width  = 6.0, 
  height = h_delta,
  device = pdf_device,
  limitsize = FALSE
)

message(glue("   Saved faceted forest: {out_delta_faceted}"))
message("\n✅ Faceted ΔAUC forest complete.")