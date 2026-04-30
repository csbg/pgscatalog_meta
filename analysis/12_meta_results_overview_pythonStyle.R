#!/usr/bin/env Rscript
## ============================================================
## 12_meta_results_overview_selfcontained_repro.R
##
## Self-contained final Step 12 plotting script, with reproduction filters matching the older overview workflow.
##
## Intended order:
##   1) source("12_meta_roadmap_two_stages.R")
##   2) source("12_i_square_filter.R")
##   3) source("12_meta_results_overview_selfcontained_repro.R")
##
## What this script does:
##   - Loads the I²-filtered Stage-1 table produced by 12_i_square_filter.R.
##   - Falls back to the unfiltered Stage-1 table only if the filtered table is absent.
##   - Rebuilds the paired weighted ΔAUC results from the Stage-1 table.
##   - Rebuilds the plotting objects that older scripts expected to already exist.
##   - Saves a faceted ΔAUC forest plot without requiring the wrapper scripts.
##   - Reproduces the old trait-selection gate:
##       selected buckets -> heatmap rows -> traits with any n_pairs >= 2 -> forest.
##
## Main input:
##   meta_roadmap_two_stages/<STAMP>/tables/stage1_pooled_cells_I2filtered_<STAMP>.csv
##
## Main outputs:
##   meta_roadmap_two_stages/<STAMP>/tables/paired_weightedtest_by_bucket_rebuilt_<STAMP>.csv
##   meta_roadmap_two_stages/<STAMP>/tables/reproduction_gate_heatmap_df_<STAMP>.csv
##   meta_roadmap_two_stages/<STAMP>/tables/reproduction_gate_keep_traits_<STAMP>.csv
##   meta_roadmap_two_stages/<STAMP>/tables/deltaAUC_forest_faceted_data_<STAMP>.csv
##   meta_roadmap_two_stages/<STAMP>/plots/results/deltaAUC_forest_faceted_<STAMP>.pdf
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(glue)
  library(stringr)
  library(forcats)
  library(scales)
  library(purrr)
  library(grid)
})

message(">> Self-contained Step 12 overview: starting")

# ============================================================
# 0) Self-contained setup
# ============================================================

if (!exists("stamp") || !is.character(stamp) || length(stamp) != 1) {
  stamp <- format(Sys.Date(), "%Y%m%d")
}

root_pw     <- file.path("meta_roadmap_two_stages", stamp)
tables_pw   <- file.path(root_pw, "tables")
plots_pw    <- file.path(root_pw, "plots")
results_dir <- file.path(plots_pw, "results")

invisible(lapply(
  list(tables_pw, plots_pw, results_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf

`%||%` <- function(a, b) if (!is.null(a)) a else b

slugify <- function(x) stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
wrap_trait <- function(x, width = 40) stringr::str_wrap(x, width = width)

# Evaluation-ancestry order used throughout the manuscript figures.
display_levels <- c(
  "European", "African", "East Asian", "South Asian",
  "Hispanic or Latin American", "Middle Eastern or North African",
  "Other/Mixed", "Not reported",
  "Multi-ancestry including European",
  "Multi-ancestry excluding European"
)

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

bucket_panel_labels <- function(x) {
  dplyr::case_when(
    grepl("European-only|EUR-only", x, ignore.case = TRUE) ~ "Weights: EUR-only",
    grepl("Multi.*incl|incl.*EUR", x, ignore.case = TRUE)  ~ "Weights: Multi incl. EUR",
    grepl("Multi.*excl|excl.*EUR", x, ignore.case = TRUE)  ~ "Weights: Multi excl. EUR",
    TRUE ~ paste("Weights:", x)
  )
}

# Lightweight version of your manuscript theme.
theme_CrossAncestryGenPhen <- function(
    base_size    = 6,
    show_axis    = TRUE,
    show_facets  = TRUE,
    show_grid    = FALSE,
    show_borders = FALSE,
    rotate       = NULL,
    base_family  = "Arial",
    ...
) {
  p <- ggplot2::theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = NA),
    strip.placement  = "inside",
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(1, "mm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = if (show_grid) element_line(color = "grey80", linewidth = 0.3) else element_blank(),
    axis.text      = element_text(size = base_size, family = base_family),
    axis.title     = element_text(size = base_size, family = base_family),
    plot.title     = element_text(size = base_size, hjust = 0.5, family = base_family),
    plot.subtitle  = element_text(size = base_size, hjust = 0.5, family = base_family),
    strip.text     = element_text(size = base_size, family = base_family),
    plot.caption   = element_text(size = base_size, hjust = 0, family = base_family),
    legend.title   = element_text(size = base_size, family = base_family),
    legend.text    = element_text(size = base_size, family = base_family),
    legend.background = element_rect(fill = NA, color = NA),
    ...
  )
  
  if (show_borders) {
    p <- p + theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line    = element_blank()
    )
  }
  
  if (!show_axis) {
    p <- p + theme(axis.line = element_blank(), axis.ticks = element_blank())
  }
  
  if (!show_facets) {
    p <- p + theme(strip.text = element_blank())
  }
  
  if (!is.null(rotate)) {
    p <- p + theme(axis.text.x = element_text(angle = rotate, hjust = 1, vjust = 1))
  }
  
  p
}

# ============================================================
# 1) Load I²-filtered Stage-1 table
# ============================================================

filtered_csv <- file.path(tables_pw, paste0("stage1_pooled_cells_I2filtered_", stamp, ".csv"))
unfiltered_csv <- file.path(tables_pw, paste0("stage1_pooled_cells_", stamp, ".csv"))

if (file.exists(filtered_csv)) {
  input_csv <- filtered_csv
  message(">> Using I²-filtered Stage-1 table: ", input_csv)
} else if (file.exists(unfiltered_csv)) {
  input_csv <- unfiltered_csv
  warning("I²-filtered table not found. Using unfiltered Stage-1 table: ", input_csv)
} else {
  stop(
    "Could not find Stage-1 input table.\nTried:\n  ",
    filtered_csv, "\n  ", unfiltered_csv,
    "\nRun 12_meta_roadmap_two_stages.R and then 12_i_square_filter.R first."
  )
}

stage1_tbl <- readr::read_csv(input_csv, show_col_types = FALSE)

# ------------------------------------------------------------
# Rescue training-bucket metadata if Stage-1 does not contain it.
# This is expected for Stage-1 tables produced by 12_meta_roadmap_two_stages.R,
# because those tables are grouped by trait × PGS × evaluation ancestry and may
# not carry train_bucket/trained_bucket forward.
# ------------------------------------------------------------
if ("train_bucket" %in% names(stage1_tbl) && !"trained_bucket" %in% names(stage1_tbl)) {
  stage1_tbl <- stage1_tbl %>% dplyr::rename(trained_bucket = train_bucket)
}

if (!"trained_bucket" %in% names(stage1_tbl)) {
  message(">> Stage-1 table has no trained_bucket column; rescuing from eval_df_final_auc_ci.csv")
  
  # You can override this before sourcing the script if needed, e.g.:
  # audit_dir <- "results/pgs_auc_ci_audit_20251006"
  if (!exists("audit_dir") || !is.character(audit_dir) || length(audit_dir) != 1) {
    audit_dir <- file.path("results", "pgs_auc_ci_audit_20251006")
  }
  
  candidate_roots <- unique(c(
    audit_dir,
    "results",
    ".",
    path.expand("~/pgscatalog/results"),
    path.expand("~/pgscatalog")
  ))
  candidate_roots <- candidate_roots[dir.exists(candidate_roots)]
  
  cand_meta <- unique(unlist(lapply(candidate_roots, function(root) {
    list.files(root, pattern = "eval_df_final_auc_ci\\.csv$", full.names = TRUE, recursive = TRUE)
  })))
  
  if (length(cand_meta) == 0) {
    stop(
      "Stage-1 table does not contain trained_bucket, and I could not find eval_df_final_auc_ci.csv to rescue it.\n",
      "Tried searching under: ", paste(candidate_roots, collapse = ", "), "\n",
      "Fix options:\n",
      "  1) Set audit_dir <- 'path/to/results/pgs_auc_ci_audit_YYYYMMDD' before sourcing this script; or\n",
      "  2) Copy eval_df_final_auc_ci.csv somewhere under results/; or\n",
      "  3) Add train_bucket/trained_bucket to the Stage-1 table upstream."
    )
  }
  
  fp_meta <- cand_meta[order(file.info(cand_meta)$mtime, decreasing = TRUE)][1]
  message(">> Using training metadata: ", fp_meta)
  
  meta_names <- names(readr::read_csv(fp_meta, n_max = 0, show_col_types = FALSE))
  if (!all(c("pgs_id", "train_bucket") %in% meta_names)) {
    stop(
      "Found eval_df_final_auc_ci.csv, but it does not contain both pgs_id and train_bucket: ", fp_meta
    )
  }
  
  bucket_map <- readr::read_csv(
    fp_meta,
    col_select = c("pgs_id", "train_bucket"),
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(
      pgs_id = as.character(pgs_id),
      trained_bucket = as.character(train_bucket)
    ) %>%
    dplyr::filter(!is.na(pgs_id), !is.na(trained_bucket), trained_bucket != "") %>%
    dplyr::distinct(pgs_id, trained_bucket)
  
  bucket_conflicts <- bucket_map %>%
    dplyr::count(pgs_id, name = "n_buckets") %>%
    dplyr::filter(n_buckets > 1)
  
  if (nrow(bucket_conflicts) > 0) {
    conflict_fp <- file.path(tables_pw, paste0("training_bucket_conflicts_", stamp, ".csv"))
    readr::write_csv(
      bucket_map %>% dplyr::semi_join(bucket_conflicts, by = "pgs_id"),
      conflict_fp
    )
    
    # For reproducibility, do not silently choose one bucket if the metadata
    # maps the same PGS to multiple different training buckets. The previous
    # wrapper could accidentally duplicate rows in this situation; stopping is
    # safer and makes the source of non-reproducibility visible.
    if (!exists("allow_bucket_conflicts_first") || !isTRUE(allow_bucket_conflicts_first)) {
      stop(
        "Some PGS IDs map to multiple training buckets. Wrote conflict table: ", conflict_fp, "\n",
        "Resolve the conflicts upstream, or set allow_bucket_conflicts_first <- TRUE before sourcing this script to keep the first bucket per PGS."
      )
    }
    
    warning(
      "Some PGS IDs have more than one train_bucket in metadata. Because allow_bucket_conflicts_first=TRUE, keeping the first bucket per PGS. Conflicting PGS count: ",
      nrow(bucket_conflicts)
    )
  }
  
  bucket_map <- bucket_map %>%
    dplyr::group_by(pgs_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(pgs_id, trained_bucket)
  stage1_tbl <- stage1_tbl %>%
    dplyr::mutate(pgs_id = as.character(pgs_id)) %>%
    dplyr::left_join(bucket_map, by = "pgs_id")
  
  n_missing_bucket <- sum(is.na(stage1_tbl$trained_bucket) | stage1_tbl$trained_bucket == "")
  if (n_missing_bucket > 0) {
    warning(
      "Training bucket could not be rescued for ", n_missing_bucket,
      " Stage-1 rows. They will be labelled as 'Train: Unknown'."
    )
    stage1_tbl <- stage1_tbl %>%
      dplyr::mutate(trained_bucket = dplyr::if_else(
        is.na(trained_bucket) | trained_bucket == "",
        "Train: Unknown",
        trained_bucket
      ))
  }
}

required_cols <- c(
  "trait_label", "ancestry_display", "pgs_id",
  "auc_pooled", "se", "trained_bucket", "k_eval"
)
missing_cols <- setdiff(required_cols, names(stage1_tbl))

if (length(missing_cols) > 0) {
  stop(
    "The Stage-1 table is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    "\nThis script cannot be self-contained until those columns are present in the Stage-1 table."
  )
}

stage1_tbl <- stage1_tbl %>%
  mutate(
    trait_label      = as.character(trait_label),
    ancestry_display = as.character(ancestry_display),
    pgs_id           = as.character(pgs_id),
    trained_bucket   = as.character(trained_bucket),
    auc_pooled       = as.numeric(auc_pooled),
    se               = as.numeric(se),
    k_eval           = as.numeric(k_eval)
  )

message(glue(">> Loaded {nrow(stage1_tbl)} Stage-1 rows across {n_distinct(stage1_tbl$trait_label)} traits"))

# ============================================================
# 2) Rebuild weighted paired ΔAUC results from Stage-1 table
# ============================================================

safe_weighted_paired_test <- function(eur_auc, tgt_auc, eur_se, tgt_se) {
  ok <- is.finite(eur_auc) & is.finite(tgt_auc) &
    is.finite(eur_se) & is.finite(tgt_se) &
    eur_se > 0 & tgt_se > 0
  
  n_pairs <- sum(ok)
  
  if (n_pairs < 3) {
    return(list(estimate = NA_real_, p.value = NA_real_, n_pairs = n_pairs))
  }
  
  diff    <- tgt_auc[ok] - eur_auc[ok]
  diff_se <- sqrt(eur_se[ok]^2 + tgt_se[ok]^2)
  weights <- 1 / diff_se^2
  
  model <- stats::lm(diff ~ 1, weights = weights)
  
  list(
    estimate = as.numeric(stats::coef(model)[1]),
    p.value  = summary(model)$coefficients[1, 4],
    n_pairs  = n_pairs
  )
}

dat <- stage1_tbl %>%
  transmute(
    trait    = trait_label,
    ancestry = as.character(ancestry_display),
    bucket   = trained_bucket,
    pgs_id,
    auc      = auc_pooled,
    se,
    k_eval
  ) %>%
  filter(ancestry != "Not reported") %>%
  mutate(
    ancestry = factor(
      ancestry,
      levels = intersect(display_levels, unique(ancestry)),
      ordered = TRUE
    )
  )

if (!("European" %in% dat$ancestry)) {
  stop("No European ancestry present after Stage-1 pooling. Cannot compute target - EUR ΔAUC.")
}

bucket_levels <- unique(dat$bucket)
bucket_levels <- bucket_levels[!is.na(bucket_levels)]
anc_all_targets <- setdiff(levels(dat$ancestry), "European")

results <- list()

for (b in bucket_levels) {
  dat_b <- dat %>% filter(bucket == b)
  
  for (anc_t in anc_all_targets) {
    if (!(anc_t %in% dat_b$ancestry)) next
    
    dat_bt <- dat_b %>%
      filter(ancestry %in% c("European", anc_t))
    
    wide_bt <- dat_bt %>%
      select(trait, ancestry, pgs_id, auc, se, k_eval) %>%
      pivot_wider(
        names_from  = ancestry,
        values_from = c(auc, se, k_eval),
        names_sep   = "_"
      )
    
    auc_eur_col <- "auc_European"
    auc_tgt_col <- paste0("auc_", anc_t)
    se_eur_col  <- "se_European"
    se_tgt_col  <- paste0("se_", anc_t)
    k_eur_col   <- "k_eval_European"
    k_tgt_col   <- paste0("k_eval_", anc_t)
    
    needed_here <- c(auc_eur_col, auc_tgt_col, se_eur_col, se_tgt_col)
    if (!all(needed_here %in% names(wide_bt))) next
    
    res_bt <- wide_bt %>%
      group_by(trait) %>%
      group_modify(~ {
        eur_auc <- .x[[auc_eur_col]]
        tgt_auc <- .x[[auc_tgt_col]]
        eur_se  <- .x[[se_eur_col]]
        tgt_se  <- .x[[se_tgt_col]]
        
        ok <- is.finite(eur_auc) & is.finite(tgt_auc) &
          is.finite(eur_se) & is.finite(tgt_se) &
          eur_se > 0 & tgt_se > 0
        
        n_pairs <- sum(ok)
        
        if (n_pairs == 0) {
          return(tibble(
            median_delta    = NA_real_,
            estimate        = NA_real_,
            p.value         = NA_real_,
            n_pairs         = 0L,
            flag_not_pooled = FALSE
          ))
        }
        
        delta_vec    <- tgt_auc[ok] - eur_auc[ok]
        median_delta <- median(delta_vec, na.rm = TRUE)
        
        wout <- safe_weighted_paired_test(
          eur_auc = eur_auc,
          tgt_auc = tgt_auc,
          eur_se  = eur_se,
          tgt_se  = tgt_se
        )
        
        flag_not_pooled <- FALSE
        if (all(c(k_eur_col, k_tgt_col) %in% names(.x))) {
          eur_k <- .x[[k_eur_col]][ok]
          tgt_k <- .x[[k_tgt_col]][ok]
          flag_not_pooled <- all(eur_k < 2 | tgt_k < 2, na.rm = TRUE)
        }
        
        tibble(
          median_delta    = median_delta,
          estimate        = wout$estimate,
          p.value         = wout$p.value,
          n_pairs         = as.integer(wout$n_pairs),
          flag_not_pooled = flag_not_pooled
        )
      }) %>%
      ungroup() %>%
      mutate(
        ancestry_target = anc_t,
        bucket = b,
        .before = 1
      )
    
    results[[length(results) + 1]] <- res_bt
  }
}

if (length(results) == 0) {
  stop("No valid trait × bucket × ancestry combinations for ΔAUC calculation.")
}

res_pw <- bind_rows(results) %>%
  mutate(
    p.adj = {
      pvec <- p.value
      mask <- is.finite(pvec)
      padj <- rep(NA_real_, length(pvec))
      if (any(mask)) padj[mask] <- p.adjust(pvec[mask], method = "fdr")
      padj
    },
    sig = !is.na(p.adj) & p.adj < 0.05,
    flag_not_pooled = coalesce(flag_not_pooled, FALSE)
  )

out_res_pw <- file.path(tables_pw, paste0("paired_weightedtest_by_bucket_rebuilt_", stamp, ".csv"))
readr::write_csv(res_pw, out_res_pw)
message(glue(">> Wrote rebuilt paired weighted-test table: {out_res_pw}"))

# ============================================================
# 3) Build plotting object expected by overview-style plots
# ============================================================

anc_pretty <- c(
  "African"                         = "African",
  "East Asian"                      = "East Asian",
  "South Asian"                     = "South Asian",
  "Hispanic or Latin American"      = "Hispanic/LatAm",
  "Middle Eastern or North African" = "MENA",
  "Other/Mixed"                     = "Other/Mixed"
)

anc_present <- intersect(names(anc_pretty), unique(res_pw$ancestry_target))
if (length(anc_present) == 0) {
  anc_present <- unique(res_pw$ancestry_target)
}

trait_order <- res_pw %>%
  group_by(trait) %>%
  summarise(median_delta_trait = median(median_delta, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_delta_trait) %>%
  pull(trait)

trait_levels_wrapped <- unique(wrap_trait(trait_order))

lim_max <- max(abs(res_pw$median_delta), na.rm = TRUE)
if (!is.finite(lim_max) || lim_max == 0) lim_max <- 0.05
lim_max <- ceiling(lim_max * 20) / 20

res_pw2 <- res_pw %>%
  filter(ancestry_target %in% anc_present) %>%
  mutate(
    ancestry_target = factor(ancestry_target, levels = anc_present, ordered = TRUE),
    trait_wrapped = factor(wrap_trait(trait), levels = trait_levels_wrapped, ordered = TRUE),
    txt_col = ifelse(abs(median_delta) > 0.35 * lim_max, "white", "black"),
    anc_lab = forcats::fct_relabel(
      ancestry_target,
      ~ ifelse(.x %in% names(anc_pretty), unname(anc_pretty[.x]), .x)
    ),
    bucket_panel = factor(
      bucket_panel_labels(bucket),
      levels = unique(bucket_panel_labels(bucket))
    )
  )

res_pw2_nonempty <- res_pw2 %>%
  filter(!is.na(median_delta) | !is.na(estimate))

if (nrow(res_pw2_nonempty) == 0) {
  stop("No rows available for plotting after rebuilding res_pw2_nonempty.")
}

message(glue(">> Rebuilt res_pw2_nonempty with {nrow(res_pw2_nonempty)} rows"))

# ============================================================
# 4) Reproduce old overview trait-selection gate
# ============================================================

# IMPORTANT: this is the part that prevents the forest from becoming
# over-inclusive. The older overview script did not plot every pairable trait.
# It first selected the manuscript buckets, then kept only traits where the
# heatmap had at least one tile based on >=2 PGS pairs.
keep_bucket_pattern <- c("European-only", "Multi incl")

heatmap_df <- res_pw2_nonempty %>%
  dplyr::filter(grepl(paste(keep_bucket_pattern, collapse = "|"), bucket, ignore.case = TRUE))

if (nrow(heatmap_df) == 0) {
  heatmap_df <- res_pw2_nonempty %>%
    dplyr::filter(grepl(paste(keep_bucket_pattern, collapse = "|"), bucket_panel, ignore.case = TRUE))
}

if (nrow(heatmap_df) == 0) {
  stop(
    "No rows survived the manuscript bucket filter.\n",
    "keep_bucket_pattern = ", paste(keep_bucket_pattern, collapse = " | "), "\n",
    "Available buckets: ", paste(sort(unique(res_pw2_nonempty$bucket)), collapse = " | ")
  )
}

message(glue(">> Buckets selected for reproduction gate: {paste(unique(heatmap_df$bucket_panel), collapse = ', ')}"))
message(glue(">> Traits before min-2-PGS gate: {dplyr::n_distinct(heatmap_df$trait_wrapped)}"))

heatmap_df <- heatmap_df %>%
  dplyr::group_by(trait_wrapped) %>%
  dplyr::filter(any(n_pairs >= 2, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    trait_wrapped = forcats::fct_drop(trait_wrapped),
    bucket_panel  = forcats::fct_drop(bucket_panel),
    anc_lab       = forcats::fct_drop(anc_lab)
  )

if (nrow(heatmap_df) == 0) {
  stop("No traits survived the old min-2-PGS gate: any(n_pairs >= 2) within selected buckets.")
}

keep_traits <- unique(as.character(heatmap_df$trait_wrapped))

selected_traits_tbl <- heatmap_df %>%
  dplyr::distinct(trait_wrapped) %>%
  dplyr::arrange(trait_wrapped)

out_heatmap_gate <- file.path(tables_pw, paste0("reproduction_gate_heatmap_df_", stamp, ".csv"))
out_keep_traits  <- file.path(tables_pw, paste0("reproduction_gate_keep_traits_", stamp, ".csv"))
readr::write_csv(heatmap_df, out_heatmap_gate)
readr::write_csv(selected_traits_tbl, out_keep_traits)

message(glue(">> Traits after min-2-PGS gate: {length(keep_traits)}"))
message(glue(">> Wrote reproduction gate table: {out_heatmap_gate}"))
message(glue(">> Wrote selected trait list: {out_keep_traits}"))

# ============================================================
# 5) Rebuild paired_df and delta_df for the faceted forest
# ============================================================

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
    grepl(paste(keep_bucket_pattern, collapse = "|"), trained_bucket, ignore.case = TRUE)
  ) %>%
  mutate(trait_wrapped_check = wrap_trait(trait_label)) %>%
  filter(trait_wrapped_check %in% keep_traits) %>%
  select(-trait_wrapped_check)

if (nrow(scat_df) == 0) {
  stop("No Stage-1 rows left for the selected training buckets and traits. Check keep_bucket_pattern.")
}

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

message(glue(">> Forest paired data: {nrow(paired_df)} paired rows across {n_distinct(paired_df$trait_label)} traits"))

if (nrow(paired_df) == 0) {
  stop("No paired EUR/target data available for the faceted forest.")
}

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
    lo              = delta_hat - 1.96 * se_hat,
    hi              = delta_hat + 1.96 * se_hat,
    pooled_group    = all_not_pooled,
    trait_clean     = wrap_trait(trait_label),
    target_ancestry = factor(
      target_ancestry,
      levels = intersect(display_levels, unique(target_ancestry)),
      ordered = TRUE
    ),
    trained_bucket = factor(
      trained_bucket,
      levels = sort(unique(as.character(trained_bucket)))
    )
  )

if (nrow(delta_df) == 0) {
  stop("No valid delta_df rows available for the faceted forest.")
}

out_delta_data <- file.path(tables_pw, paste0("deltaAUC_forest_faceted_data_", stamp, ".csv"))
readr::write_csv(delta_df, out_delta_data)
message(glue(">> Wrote faceted forest plotting data: {out_delta_data}"))

# ============================================================
# 6) Axis limits and colors
# ============================================================

rng <- range(c(delta_df$lo, delta_df$hi), na.rm = TRUE)
max_abs <- max(abs(rng))
if (!is.finite(max_abs) || max_abs == 0) max_abs <- 0.05
lim_x <- ceiling(max_abs * 20) / 20

anc_targets <- levels(delta_df$target_ancestry)
anc_targets <- anc_targets[anc_targets %in% unique(as.character(delta_df$target_ancestry))]
anc_colors <- pal_ancestry[anc_targets]
anc_colors <- anc_colors[!is.na(anc_colors)]

# ============================================================
# 7) Faceted ΔAUC forest
# ============================================================

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
    name = "Pooling status"
  ) +
  scale_linetype_discrete(name = "Training weights") +
  scale_x_continuous(
    limits = c(-lim_x, lim_x),
    breaks = pretty(c(-lim_x, lim_x), n = 5),
    name   = expression(Delta * AUC ~ "(target - EUR)")
  ) +
  facet_grid(
    trait_clean ~ .,
    scales = "free_y",
    space  = "free_y",
    switch = "y"
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
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6, margin = margin(r = 5))
  )

# ============================================================
# 8) Save the faceted plot
# ============================================================

n_traits     <- n_distinct(delta_df$trait_clean)
n_total_rows <- n_distinct(paste(delta_df$trait_clean, delta_df$target_ancestry))
h_delta      <- max(4.0, (n_total_rows * 0.15) + (n_traits * 0.2) + 1.5)

out_delta_faceted <- file.path(results_dir, paste0("deltaAUC_forest_faceted_", stamp, ".pdf"))

ggsave(
  out_delta_faceted,
  g_delta_faceted,
  width  = 6.0,
  height = h_delta,
  device = pdf_device,
  limitsize = FALSE
)

message(glue(">> Saved faceted forest: {out_delta_faceted}"))
message("\n✅ Self-contained faceted ΔAUC forest complete.")
