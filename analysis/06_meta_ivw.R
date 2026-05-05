#!/usr/bin/env Rscript

# ============================================================
# Script: 06_meta_ivw.R
# Project: pgscatalog_meta
# Purpose: Perform Stage-1 inverse-variance weighted meta-analysis of AUC within each trait × PGS × ancestry cell.
# Inputs: results/pgs_auc_ci_audit_<DATE>/eval_df_final_auc_ci.csv
# Outputs: meta_roadmap_two_stages/<DATE>/tables/stage1_pooled_cells_<DATE>.csv
# Run after: 05_pgs_auc_ci_audit.R
# Run before: 07_i_square_filter.R
# ============================================================

setwd("~/pgscatalog")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(forcats); library(janitor); library(purrr)
  library(ggplot2); library(scales); library(jsonlite); library(glue)
  library(tibble)
})

# -------------------------
# 0) CONFIG
# -------------------------
if (!exists("audit_dir"))  audit_dir  <- "/home/people/nnunes/pgscatalog/results/pgs_auc_ci_audit_20251006"  # where eval_df_final_auc_ci.csv lives
if (!exists("stamp") || !is.character(stamp) || length(stamp) != 1) {
  stamp <- format(Sys.Date(), "%Y%m%d")
}

# drop "Not reported" by default from plots
drop_label <- "Not reported"

# keep this FALSE for a simpler run (no k_min ancestry gate)
if (!exists("apply_min_ancestries_filter")) apply_min_ancestries_filter <- FALSE
if (!exists("min_ancestries"))              min_ancestries <- 2L
if (!exists("anc_exclude_from_count"))      anc_exclude_from_count <- c("Not reported")

display_levels <- c(
  "European","African","East Asian","South Asian",
  "Hispanic or Latin American","Middle Eastern or North African",
  "Other/Mixed","Not reported",
  "Multi-ancestry including European","Multi-ancestry excluding European"
)

# Palette
if (!exists("pal_ancestry")) {
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
}

pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
`%||%` <- function(a, b) if (!is.null(a)) a else b
slugify  <- function(x) stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")

theme_clean <- function(){
  theme_minimal(base_size = 11) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      strip.text       = element_text(face = "bold"),
      panel.spacing.x  = unit(1.5, "mm"),
      panel.spacing.y  = unit(2.0, "mm"),
      plot.margin      = margin(t = 6, r = 10, b = 28, l = 12)
    )
}

# -------------------------
# 1) HELPERS (ancestry map + IVW)
# -------------------------
make_ancestry_display <- function(x){
  x <- as.character(x)
  g <- function(p) grepl(p, x, ignore.case = TRUE, perl = TRUE)
  dplyr::case_when(
    is.na(x) | x == "" | g("Unknown|Not reported|\\bNR\\b|unspecified") ~ "Not reported",
    g("Multi-ancestry including European") ~ "Multi-ancestry including European",
    g("Multi-ancestry excluding European") ~ "Multi-ancestry excluding European",
    g("^European(,|$)|^EUR(,|$)|^European\\b") ~ "European",
    g("Sub-?Saharan|\\bSSA\\b") ~ "African",
    g("African\\s*American|Afro-?Caribbean|Black or African American|Black/African American") ~ "African",
    (g("^African") & !g("^African\\s*American")) ~ "African",
    g("Middle|Greater Middle|North African|\\bMENA\\b") ~ "Middle Eastern or North African",
    g("East\\s*Asian")  ~ "East Asian",
    g("South\\s*Asian") ~ "South Asian",
    g("Hispanic|Lat(in|inx|ino)|Latin American") ~ "Hispanic or Latin American",
    g("Other|Mixed|Admixed") ~ "Other/Mixed",
    TRUE ~ x
  )
}

clamp01 <- function(a, eps = 1e-6) pmin(pmax(a, eps), 1 - eps)

auc_ci_to_logit <- function(auc, lo, hi, eps_se = 1e-8){
  a    <- clamp01(as.numeric(auc))
  lo1 <- clamp01(as.numeric(lo))
  hi1 <- clamp01(as.numeric(hi))
  eta <- qlogis(a)
  se  <- (qlogis(hi1) - qlogis(lo1)) / (2*1.96)
  # guard against zero/negative/NA se
  se  <- ifelse(is.finite(se) & se > 0, se, eps_se)
  tibble(eta = eta, se = se)
}

ivw_pool_logit <- function(eta, se){
  # keep only finite, positive SE
  ok <- is.finite(eta) & is.finite(se) & (se > 0)
  eta <- eta[ok]; se <- se[ok]
  k <- length(eta)
  if (k == 0) return(tibble(eta = NA_real_, se = NA_real_, Q = NA_real_, I2 = NA_real_, k_eval = 0L))
  if (k == 1) return(tibble(eta = eta,       se = se,       Q = NA_real_, I2 = NA_real_, k_eval = 1L))
  
  w <- 1 / (se^2)
  w <- ifelse(is.finite(w), w, 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0)
    return(tibble(eta = NA_real_, se = NA_real_, Q = NA_real_, I2 = NA_real_, k_eval = k))
  
  eta_hat <- sum(w * eta) / sw
  se_hat  <- sqrt(1 / sw)
  
  Q  <- sum(w * (eta - eta_hat)^2)
  df <- k - 1
  I2 <- if (isTRUE(Q > 0)) max(0, (Q - df) / Q) * 100 else 0
  
  tibble(eta = eta_hat, se = se_hat, Q = Q, I2 = I2, k_eval = k)
}

# Optional ancestry-count filter (off by default)
filter_min_ancestries <- function(df){
  df %>%
    mutate(anc_for_count = ifelse(ancestry_display %in% anc_exclude_from_count, NA, ancestry_display)) %>%
    group_by(trait_label, pgs_id) %>%
    mutate(n_anc_distinct = n_distinct(anc_for_count, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(n_anc_distinct >= min_ancestries) %>%
    select(-anc_for_count, -n_anc_distinct)
}

# -------------------------
# 2) OUTPUT DIRS
# -------------------------
root <- file.path("meta_roadmap_two_stages", stamp)
dirs <- list(
  base   = root,
  plots  = file.path(root, "plots"),
  raw    = file.path(root, "plots", "raw"),
  stage1 = file.path(root, "plots", "stage1"),
  stage2 = file.path(root, "plots", "stage2"),
  tables = file.path(root, "tables"),
  logs   = file.path(root, "logs")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))



# -------------------------
# 3) INPUT & PREP
# -------------------------
message(">> Looking for eval_df_final_auc_ci.csv under: ", normalizePath(audit_dir, mustWork = FALSE))
cand <- list.files(audit_dir, pattern = "eval_df_final_auc_ci\\.csv$", full.names = TRUE, recursive = TRUE)
stopifnot("No eval_df_final_auc_ci.csv found under audit_dir" = length(cand) > 0)
fp_eval_final <- cand[order(file.info(cand)$mtime, decreasing = TRUE)][1]
message("  * Using: ", fp_eval_final)

eval_final <- readr::read_csv(fp_eval_final, show_col_types = FALSE) %>% janitor::clean_names()

needed <- c("pgs_id","performance_id","sampleset_id","ancestry_eval",
            "auc","estimate_ci_lower","estimate_ci_upper",
            "reported_trait","efo_label")
missing <- setdiff(needed, names(eval_final))
if (length(missing) > 0) for (nm in missing) eval_final[[nm]] <- NA

# ---------- ORIGINAL METHOD: use EFO first, then reported_trait ----------
# Fields expected (any missing ones are created as NA so code is robust):
for (nm in c("efo_id","efo_label","reported_trait")) if (!nm %in% names(eval_final)) eval_final[[nm]] <- NA_character_

na_empty <- function(x) { x <- as.character(x); ifelse(!is.na(x) & trimws(x) == "", NA_character_, x) }
slug_lc  <- function(x) stringr::str_replace_all(tolower(x), "[^a-z0-9]+", "_")

# ======================================================================
# TRAIT HARMONIZATION (The Elegant Fix)
# ======================================================================

# 1. Define Dictionary of Aliases
trait_dictionary <- c(
  # Fix Capitalization Duplicates
  "Type 2 Diabetes"    = "Type 2 diabetes",
  "Breast Cancer"      = "Breast cancer",
  "Prostate Cancer"    = "Prostate cancer",
  
  # Fix Abbreviations / Synonyms
  "T1D"                = "Type 1 diabetes",
  "Cancer of prostate" = "Prostate cancer",
  "Breast cancer [female]" = "Breast cancer",
  
  # Standardize naming
  "Coronary artery disease" = "Coronary heart disease"
)

# 2. Helper to clean text
clean_text <- function(x) {
  trimws(x)
}

# 3. Apply Logic
eval_final <- eval_final %>%
  mutate(
    efo_id      = na_empty(efo_id),
    efo_label   = na_empty(efo_label),
    reported_trait = na_empty(reported_trait),
    
    # Step A: Get raw label (EFO preferred)
    raw_label = dplyr::coalesce(efo_label, reported_trait),
    
    # Step B: Basic Clean
    temp_label = clean_text(raw_label),
    
    # Step C: Apply Dictionary Remapping (The Magic Step)
    trait_label = dplyr::recode(temp_label, !!!trait_dictionary),
    
    # Step D: Re-generate ID based on the CLEAN label
    # This ensures "Type 2 Diabetes" and "Type 2 diabetes" get merged into ONE ID.
    trait_id = paste0("custom:", slug_lc(trait_label))
  )
# Proceed with processing; keep both id & label for grouping/printing
eval_proc <- eval_final %>%
  mutate(
    ancestry_display = make_ancestry_display(ancestry_eval)
  ) %>%
  mutate(across(c(auc, estimate_ci_lower, estimate_ci_upper), as.numeric)) %>%
  filter(is.finite(auc), is.finite(estimate_ci_lower), is.finite(estimate_ci_upper))

# Audit: show how rows map to EFO vs reported
dir.create(dirs$tables, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  eval_proc %>%
    transmute(trait_id, trait_label, efo_id, efo_label, reported_trait) %>%
    count(trait_id, trait_label, efo_id, efo_label, reported_trait, name = "n_rows") %>%
    arrange(desc(n_rows)),
  file.path(dirs$tables, paste0("trait_mapping_audit_original_method_", stamp, ".csv"))
)

mapping_check <- eval_proc %>%
  transmute(ancestry_raw = ancestry_eval, ancestry_display) %>%
  count(ancestry_raw, ancestry_display, sort = TRUE)

# Save ancestry mapping check
write_csv(mapping_check, file.path(dirs$tables, paste0("ancestry_display_mapping_", stamp, ".csv")))
# Optional ancestry filter
eval_for_stage <- if (isTRUE(apply_min_ancestries_filter)) filter_min_ancestries(eval_proc) else eval_proc

# -------------------------
# 4) PLOT 1 — RAW (flipped)
# -------------------------
plot_raw_flipped <- function(df){
  message(">> Plot 1/3: RAW (flipped)")
  df <- df %>%
    mutate(ancestry_display = as.character(ancestry_display)) %>%
    filter(ancestry_display != drop_label)
  
  if (nrow(df) == 0) return(invisible(NULL))
  
  # keep observed ancestry order but constrained by display_levels
  lvl_anc <- intersect(display_levels, unique(df$ancestry_display))
  df <- df %>%
    mutate(ancestry_display = factor(ancestry_display, levels = lvl_anc, ordered = TRUE))
  
  traits <- sort(unique(df$trait_label))
  purrr::walk(traits, function(tr){
    d <- df %>% filter(trait_label == tr)
    if (nrow(d) == 0) return(invisible(NULL))
    
    # stable PGS order
    pgs_levels <- d %>%
      group_by(pgs_id) %>% summarise(m = median(auc), .groups="drop") %>%
      arrange(m) %>% pull(pgs_id)
    
    # how many facets and rows per facet (for dynamic height)
    n_anc <- dplyr::n_distinct(d$ancestry_display)
    max_rows_per_anc <- d %>% distinct(ancestry_display, pgs_id) %>% count(ancestry_display) %>% pull(n) %>% max()
    
    g <- ggplot(d, aes(y = factor(pgs_id, levels = pgs_levels), x = auc, color = ancestry_display)) +
      geom_errorbarh(aes(xmin = estimate_ci_lower, xmax = estimate_ci_upper), height = 0, alpha = 0.9) +
      geom_point(size = 1.9) +
      geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.25, alpha = 0.5) +
      facet_grid(rows = vars(ancestry_display), scales = "free_y") +   # <<< rows, not columns
      scale_color_manual(values = pal_ancestry, limits = lvl_anc, guide = "none") +
      coord_cartesian(xlim = c(0.45, 1.0)) +
      labs(title = glue("RAW AUC — {tr}"), x = "AUC (95% CI)", y = "PGS") +
      theme_clean() +
      theme(strip.text.y = element_text(angle = 0, hjust = 0.5))
    
    # dynamic size: more ancestries/rows => taller
    width_this  <- 9.5
    height_this <- max(5.0, 2.0 + 1.1 * n_anc + 0.16 * max_rows_per_anc)
    fn <- file.path(dirs$raw, paste0("raw_flipped_", slugify(tr), "_", stamp, ".pdf"))
    ggsave(fn, g, width = width_this, height = height_this, device = pdf_device, limitsize = FALSE)
    
    write_csv(d %>% arrange(ancestry_display, pgs_id),
              file.path(dirs$tables, paste0("raw_rows_", slugify(tr), "_", stamp, ".csv")))
  })
}

# -------------------------
# 5) STAGE-1 POOLING (within ancestry per PGS) + plot
# -------------------------
pool_stage1 <- function(df){
  df %>%
    mutate(tmp = pmap(list(auc, estimate_ci_lower, estimate_ci_upper), auc_ci_to_logit)) %>%
    tidyr::unnest_wider(tmp) %>%
    # extra safety: drop any non-finite logits/SE before grouping
    filter(is.finite(eta), is.finite(se), se > 0) %>%
    group_by(pgs_id, trait_label, ancestry_display) %>%
    summarise(pooled = list(ivw_pool_logit(eta, se)), .groups = "drop") %>%
    tidyr::unnest_wider(pooled) %>%
    mutate(
      auc_pooled = plogis(eta),
      lo_pooled  = plogis(eta - 1.96*se),
      hi_pooled  = plogis(eta + 1.96*se)
    )
}

# -------------------------
# 5) STAGE-1 POOLING: FINAL (Grouped by Ancestry, One Eval per Line)
# -------------------------
plot_stage1_final <- function(pooled){
  message(">> Plot 2/3: Stage-1 Final (Grouped by Ancestry)")
  
  # 1. Clean Data
  df <- pooled %>%
    mutate(ancestry_display = as.character(ancestry_display)) %>%
    filter(ancestry_display != drop_label)
  
  if (nrow(df) == 0) return(invisible(NULL))
  
  # 2. Loop Per Trait
  traits <- sort(unique(df$trait_label))
  
  purrr::walk(traits, function(tr){
    d_trait <- df %>% filter(trait_label == tr)
    if (nrow(d_trait) == 0) return(invisible(NULL))
    
    # 3. Sort PGS IDs by European Performance (for consistency)
    pgs_order <- d_trait %>%
      group_by(pgs_id) %>%
      summarise(m = median(auc_pooled[ancestry_display == "European"], na.rm = TRUE)) %>%
      # Fallback if no European data
      mutate(m = ifelse(is.na(m), median(d_trait$auc_pooled), m)) %>%
      arrange(m) %>% 
      pull(pgs_id)
    
    d_trait <- d_trait %>%
      mutate(
        pgs_id = factor(pgs_id, levels = pgs_order),
        # Ensure Ancestry Factor order (European on Top)
        ancestry_display = factor(ancestry_display, levels = display_levels, ordered = TRUE)
      )
    
    # 4. Dynamic Height
    # We count the total number of Evaluations (rows) across all ancestries
    total_evals <- nrow(d_trait)
    
    # Height: ~0.2 inches per line is enough for "one eval per line"
    h_dynamic <- 1.5 + (total_evals * 0.25)
    h_dynamic <- max(5, h_dynamic)
    
    # 5. Plot
    g <- ggplot(d_trait, aes(x = auc_pooled, y = pgs_id, color = ancestry_display)) +
      
      # Reference Line
      geom_vline(xintercept = 0.5, linetype = "solid", color = "gray85", linewidth = 0.5) +
      
      # Error Bars & Dots
      geom_errorbarh(aes(xmin = lo_pooled, xmax = hi_pooled), height = 0, linewidth = 0.6, alpha = 0.8) +
      geom_point(size = 2.0) +
      
      # FACET BY ANCESTRY
      # scales="free_y" removes empty rows. space="free_y" resizes panels to fit data.
      facet_grid(ancestry_display ~ ., scales = "free_y", space = "free_y", switch = "y") +
      
      # Colors
      scale_color_manual(values = pal_ancestry, guide = "none") +
      scale_x_continuous(breaks = seq(0.5, 1.0, 0.1), limits = c(0.45, 1.0)) +
      
      labs(
        title = glue("Stage 1: {tr}"),
        subtitle = glue("Pooled AUC grouped by Ancestry (Total {total_evals} evaluations)."),
        x = "Pooled AUC (95% CI)", y = NULL
      ) +
      
      theme_minimal(base_size = 12) +
      theme(
        # Headers (Ancestry Names) on the LEFT side (switch=y)
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 11, hjust = 1),
        strip.background = element_rect(fill = "gray96", color = NA),
        strip.placement = "outside", # Moves ancestry label outside axis text
        
        # Panel spacing
        panel.spacing = unit(0.5, "lines"),
        panel.border = element_rect(color = "gray80", fill = NA),
        
        # Axis Text
        axis.text.y = element_text(size = 9, color = "black"),
        panel.grid.major.y = element_line(color = "gray95", linetype = "dotted")
      )
    
    # 6. Save (NEW FILENAME to distinguish it)
    fn <- file.path(dirs$stage1, paste0("stage1_grouped_ancestries_", slugify(tr), "_", stamp, ".pdf"))
    
    ggsave(fn, g, width = 8.5, height = h_dynamic, device = pdf_device, limitsize = FALSE)
    
    write_csv(d_trait, file.path(dirs$tables, paste0("stage1_rows_", slugify(tr), "_", stamp, ".csv")))
  })
}


# -------------------------
# 5b) Combined RAW + Stage-1 pooling per trait (MODIFIED)
# -------------------------
plot_pooling_combined <- function(raw_df, stage1_tbl){
  message(">> Plot: Combined RAW evaluations + Stage-1 pooled (one file per trait)")
  message("   * STRICT LOGIC: Raw dots and Pooled bars must BOTH exist to be plotted.")
  
  # Drop ancestries that we do not want to display
  raw_df <- raw_df %>%
    mutate(ancestry_display = as.character(ancestry_display)) %>%
    filter(ancestry_display != drop_label)
  
  stage1_df <- stage1_tbl %>%
    mutate(ancestry_display = as.character(ancestry_display)) %>%
    filter(ancestry_display != drop_label)
  
  if (nrow(raw_df) == 0 || nrow(stage1_df) == 0) {
    message("No data left after filtering for combined pooling plot.")
    return(invisible(NULL))
  }
  
  # Traits present after filtering
  traits <- sort(unique(stage1_df$trait_label))
  
  purrr::walk(traits, function(tr){
    d_raw    <- raw_df    %>% filter(trait_label == tr)
    d_pooled <- stage1_df %>% filter(trait_label == tr)
    
    if (nrow(d_raw) == 0 || nrow(d_pooled) == 0) return(invisible(NULL))
    
    # =================================================================
    # START MODIFICATION: STRICT INTERSECTION (A AND B)
    # =================================================================
    # 1. Identify keys available in RAW data for this trait
    keys_raw <- d_raw %>% 
      select(pgs_id, ancestry_display) %>% 
      distinct()
    
    # 2. Identify keys available in POOLED data (that are valid)
    keys_pooled <- d_pooled %>% 
      filter(!is.na(auc_pooled)) %>% 
      select(pgs_id, ancestry_display) %>% 
      distinct()
    
    # 3. Find the Intersection: Keys present in BOTH
    valid_keys <- dplyr::inner_join(keys_raw, keys_pooled, by = c("pgs_id", "ancestry_display"))
    
    # 4. Filter BOTH datasets to only include these shared keys
    # This prevents "Raw without Pooled" AND "Pooled without Raw"
    d_raw <- d_raw %>% semi_join(valid_keys, by = c("pgs_id", "ancestry_display"))
    d_pooled <- d_pooled %>% semi_join(valid_keys, by = c("pgs_id", "ancestry_display"))
    
    # =================================================================
    # END MODIFICATION
    # =================================================================
    
    # If filtering emptied the dataframes
    if(nrow(d_raw) == 0 || nrow(d_pooled) == 0) return(invisible(NULL))
    
    # Ancestry order restricted to those present
    lvl_anc <- intersect(display_levels, unique(c(d_raw$ancestry_display, d_pooled$ancestry_display)))
    d_raw <- d_raw %>%
      mutate(ancestry_display = factor(ancestry_display, levels = lvl_anc, ordered = TRUE))
    d_pooled <- d_pooled %>%
      mutate(ancestry_display = factor(ancestry_display, levels = lvl_anc, ordered = TRUE))
    
    # Order PGS IDs by Stage-1 pooled performance (prefer European, otherwise overall)
    pgs_order_tbl <- d_pooled %>%
      group_by(pgs_id) %>%
      summarise(
        m = median(auc_pooled[ancestry_display == "European"], na.rm = TRUE),
        .groups = "drop"
      )
    
    global_fallback <- median(d_pooled$auc_pooled, na.rm = TRUE)
    if (is.na(global_fallback)) {
      global_fallback <- 0.6
    }
    
    pgs_order <- pgs_order_tbl %>%
      mutate(m = ifelse(is.na(m), global_fallback, m)) %>%
      arrange(m) %>%
      pull(pgs_id)
    
    d_raw <- d_raw %>%
      mutate(pgs_id = factor(pgs_id, levels = pgs_order))
    d_pooled <- d_pooled %>%
      mutate(pgs_id = factor(pgs_id, levels = pgs_order))
    
    # Dynamic height based on valid PGS × ancestry combinations
    total_lines <- d_pooled %>%
      distinct(ancestry_display, pgs_id) %>%
      nrow()
    n_anc <- dplyr::n_distinct(d_pooled$ancestry_display)
    h_dynamic <- max(5, 2.0 + 1.1 * n_anc + 0.16 * total_lines)
    
    g <- ggplot() +
      # Reference line
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray80", linewidth = 0.4) +
      
      # --- RAW evaluations (transparent) ---
      geom_errorbarh(
        data = d_raw,
        aes(y = pgs_id, x = auc,
            xmin = estimate_ci_lower, xmax = estimate_ci_upper,
            color = ancestry_display),
        height   = 0,
        linewidth = 0.4,
        alpha    = 0.35
      ) +
      geom_point(
        data = d_raw,
        aes(y = pgs_id, x = auc, color = ancestry_display),
        size  = 1.4,
        alpha = 0.35
      ) +
      
      # --- Stage-1 pooled (solid, no transparency) ---
      geom_errorbarh(
        data = d_pooled,
        aes(y = pgs_id, x = auc_pooled,
            xmin = lo_pooled, xmax = hi_pooled,
            color = ancestry_display),
        height   = 0,
        linewidth = 0.8
      ) +
      geom_point(
        data = d_pooled,
        aes(y = pgs_id, x = auc_pooled, color = ancestry_display),
        size = 2.0
      ) +
      
      facet_grid(ancestry_display ~ ., scales = "free_y", space = "free_y", switch = "y") +
      scale_color_manual(values = pal_ancestry, limits = lvl_anc, guide = "none") +
      coord_cartesian(xlim = c(0.45, 1.0)) +
      labs(
        title    = glue("Raw and Stage-1 pooled AUC — {tr}"),
        subtitle = "Transparent: raw cohort evaluations; Solid: Stage-1 pooled (Matched)",
        x        = "AUC (95% CI)",
        y        = "PGS"
      ) +
      theme_clean() +
      theme(
        strip.text.y   = element_text(angle = 0, hjust = 0.5),
        panel.spacing  = unit(0.5, "lines"),
        panel.border   = element_rect(color = "gray80", fill = NA),
        axis.text.y    = element_text(size = 9, color = "black"),
        panel.grid.major.y = element_line(color = "gray95", linetype = "dotted")
      )
    
    fn <- file.path(dirs$stage1, paste0("combined_raw_and_stage1_", slugify(tr), "_", stamp, ".pdf"))
    ggsave(fn, g, width = 9.5, height = h_dynamic, device = pdf_device, limitsize = FALSE)
    
    message(glue("  Saved combined pooling plot for {tr}: {fn}"))
  })
}

# -------------------------
# 6) STAGE-2 POOLING (across PGS per trait × ancestry) + plot
# -------------------------
pool_stage2 <- function(stage1_tbl){
  stage1_tbl %>%
    select(trait_label, ancestry_display, pgs_id, eta, se, I2, k_eval) %>%
    group_by(trait_label, ancestry_display) %>%
    summarise(pooled = list(ivw_pool_logit(eta, se)),
              median_I2_pgs = median(I2, na.rm = TRUE),
              .groups = "drop") %>%
    unnest_wider(pooled) %>%
    mutate(
      auc_pool = plogis(eta),
      lo_pool  = plogis(eta - 1.96*se),
      hi_pool  = plogis(eta + 1.96*se)
    )
}

# -------------------------
# 6b) IMPROVED STAGE-2 PLOT (Individual Files per Trait)
# -------------------------
plot_stage2_individual_traits <- function(stage2_tbl){
  message(">> Plot 3/3: Stage-2 Comparative (Individual Files)")
  
  # 1. Clean Data
  d <- stage2_tbl %>%
    mutate(ancestry_display = as.character(ancestry_display)) %>%
    filter(ancestry_display != drop_label)
  
  if (nrow(d) == 0) return(invisible(NULL))
  
  # 2. STRICT FILTER (Must have EUR + Non-EUR)
  traits_stats <- d %>%
    group_by(trait_label) %>%
    summarise(
      has_eur = "European" %in% ancestry_display,
      has_non_eur = any(ancestry_display != "European"),
      n_ancestries = n_distinct(ancestry_display),
      total_k = sum(k_eval), 
      .groups = "drop"
    ) %>%
    filter(has_eur & has_non_eur)
  
  n_kept <- nrow(traits_stats)
  message(glue("   * Generating plots for {n_kept} traits..."))
  
  if (n_kept == 0) return(invisible(NULL))
  
  # 3. PREPARE PLOT DATA
  d_plot <- d %>%
    filter(trait_label %in% traits_stats$trait_label) %>%
    left_join(traits_stats %>% select(trait_label, n_ancestries), by = "trait_label") %>%
    mutate(
      anc_label = glue("{ancestry_display} (k={k_eval})")
    )
  
  # Calculate European References
  eur_refs <- d_plot %>%
    filter(ancestry_display == "European") %>%
    select(trait_label, eur_auc = auc_pool)
  
  d_plot <- d_plot %>% left_join(eur_refs, by = "trait_label")
  
  # 4. LOOP & SAVE INDIVIDUALLY
  # Base directory
  base_out <- file.path(dirs$stage2, "by_ancestry_count")
  dir.create(base_out, recursive = TRUE, showWarnings = FALSE)
  
  # We iterate through every single trait
  traits_list <- sort(unique(d_plot$trait_label))
  
  purrr::walk(traits_list, function(tr){
    
    # Subset
    d_sub <- d_plot %>% filter(trait_label == tr)
    
    # Determine which sub-folder it belongs to (e.g. "3_ancestries")
    N <- unique(d_sub$n_ancestries)
    sub_dir <- file.path(base_out, paste0(N, "_ancestries"))
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Plot
    g <- ggplot(d_sub, aes(x = auc_pool, y = ancestry_display, color = ancestry_display)) +
      
      # Reference Line
      geom_vline(aes(xintercept = eur_auc), linetype = "dashed", color = "black", alpha = 0.6) +
      
      # Error Bar
      geom_errorbarh(aes(xmin = lo_pool, xmax = hi_pool), height = 0.3, linewidth = 0.8) +
      
      # Dot
      geom_point(aes(size = k_eval), alpha = 0.9) +
      
      # Facet (Just for the title strip style, essentially)
      facet_wrap(~trait_label, scales = "free") +
      
      # Scales
      scale_color_manual(values = pal_ancestry, guide = "none") +
      scale_size_continuous(range = c(3, 8), guide = "none") + # Slightly larger dots for single plot
      scale_x_continuous(breaks = seq(0.4, 1.0, 0.05)) +
      
      labs(
        title = glue("Comparative Accuracy: {tr}"),
        subtitle = "Vertical Dashed Line = European Mean.",
        x = "Pooled AUC (95% CI)", y = NULL
      ) +
      
      theme_minimal(base_size = 12) +
      theme(
        panel.border = element_rect(color = "gray80", fill = NA),
        strip.text = element_text(face = "bold", size = 12, hjust = 0),
        strip.background = element_rect(fill = "gray97"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 11, color = "gray20")
      )
    
    # Save as Individual File
    # Clean filename: "Gout.pdf"
    clean_name <- slugify(tr)
    fn <- file.path(sub_dir, paste0(clean_name, "_compare_", stamp, ".pdf"))
    
    ggsave(fn, g, width = 8, height = 5, device = pdf_device)
  })
  
  message(glue("✅ Saved {length(traits_list)} individual plots in: {base_out}"))
}

# -------------------------
# 7) RUN
# -------------------------
# Log params
write(jsonlite::toJSON(list(
  timestamp   = as.character(Sys.time()),
  stamp       = stamp,
  input_csv   = fp_eval_final,
  audit_dir   = audit_dir,
  drop_label  = drop_label,
  apply_min_ancestries_filter = apply_min_ancestries_filter,
  min_ancestries = min_ancestries
), pretty = TRUE, auto_unbox = TRUE),
file = file.path(dirs$logs, "params.json"))

# 1) Stage-1 pooling (within ancestry × PGS)
stage1_tbl <- pool_stage1(eval_for_stage)

# Save Stage-1 pooled cells for downstream I² + Wilcoxon scripts
out_stage1 <- file.path(dirs$tables, paste0("stage1_pooled_cells_", stamp, ".csv"))
readr::write_csv(stage1_tbl, out_stage1)
message(glue(">> Wrote Stage-1 pooled cells: {out_stage1}"))

# 2) Combined RAW + pooled plot per trait (boss feedback: single pooling plot per phenotype)
plot_pooling_combined(eval_for_stage, stage1_tbl)

# 3) Stage-2 pooling across PGS (kept as in original script; comment out if not needed)
stage2_tbl <- pool_stage2(stage1_tbl)
plot_stage2_individual_traits(stage2_tbl)

message(glue("\n✅ Done. Outputs in: {normalizePath(root, mustWork = FALSE)}"))
