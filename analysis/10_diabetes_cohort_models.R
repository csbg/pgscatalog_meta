#!/usr/bin/env Rscript

# ============================================================
# Script: 10_diabetes_cohort_models.R
# Project: pgscatalog_meta
# Purpose: Type 2 diabetes side analysis. Fit ancestry models that
#          mark shared sample sets, write an HTML page for Georg, and leave
#          script 09 / DerSimonian-Laird unchanged.
#
# Reader-facing words, matching tests/q_random_effects_t2d_illustrations.R:
#   PGS, evaluation, sample set, pair, comparison, training set.
#   In the model formula, (1|cohort) means the sample set. Say that once.
#   Do not write cell, code, cohorts left, key, bucket, or k on the page.
# Inputs:  eval_df_final_auc_ci, Stage-1 I2-filtered cells, sample sets,
#          09_dependence_cells.csv
# Outputs: results/review_sensitivity_<STAMP>/diabetes_cohort_models/
# Run from the repository root after analysis/09_review_sensitivity.R.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(janitor)
  library(lme4)
  library(metafor)
  library(jsonlite)
})

if (!file.exists("R/load.R")) {
  stop("Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

stamp <- pipeline_stamp
out_dir <- file.path(
  "results", paste0("review_sensitivity_", stamp), "diabetes_cohort_models"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmt3 <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 3), "NA")
fmt4 <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 4), "NA")
html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

normalize_study <- function(x) {
  x0 <- trimws(as.character(x))
  x0[is.na(x0) | !nzchar(x0)] <- "(blank)"
  xl <- tolower(x0)
  out <- x0
  out[xl %in% c("ukb", "uk biobank") | grepl("^uk\\s*biobank", xl)] <- "UK Biobank"
  out[xl %in% c("g&h", "g and h") | grepl("genes\\s*&\\s*health|genes and health", xl)] <- "Genes & Health"
  out[xl %in% c("allofus", "all of us") | grepl("all\\s*of\\s*us", xl)] <- "All of Us"
  out
}

pool_delta_check <- function(delta, se, model = c("fixed", "random_dl")) {
  model <- match.arg(model)
  w <- 1 / se^2
  sw <- sum(w)
  dh <- sum(w * delta) / sw
  Q <- sum(w * (delta - dh)^2)
  k <- length(delta)
  df <- k - 1L
  tau2 <- 0
  if (model == "random_dl" && k >= 2) {
    C <- sw - sum(w^2) / sw
    tau2 <- if (is.finite(C) && C > 0) max(0, (Q - df) / C) else 0
  }
  ws <- 1 / (se^2 + tau2)
  dh2 <- sum(ws * delta) / sum(ws)
  se_hat <- sqrt(1 / sum(ws))
  I2 <- if (k >= 2 && Q > 0) max(0, (Q - df) / Q) * 100 else NA_real_
  tibble(
    delta_hat = dh2, se_hat = se_hat,
    lo = dh2 - 1.96 * se_hat, hi = dh2 + 1.96 * se_hat,
    tau2 = tau2, Q = Q, df = df, I2 = I2, k = k
  )
}

safe_coef_ci <- function(est, se, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  tibble(
    estimate = as.numeric(est),
    se = as.numeric(se),
    lo = as.numeric(est) - z * as.numeric(se),
    hi = as.numeric(est) + z * as.numeric(se)
  )
}

extract_ancestry_lm <- function(fit, model_label, cohort_def, note = NA_character_) {
  sm <- summary(fit)$coefficients
  rn <- rownames(sm)
  keep <- grepl("^ancestry", rn)
  if (!any(keep)) {
    return(tibble(
      model = model_label, cohort_def = cohort_def, ancestry = character(),
      estimate = numeric(), se = numeric(), lo = numeric(), hi = numeric(),
      excludes_zero = logical(), note = character()
    ))
  }
  anc <- sub("^ancestry", "", rn[keep])
  bind_cols(
    tibble(model = model_label, cohort_def = cohort_def, ancestry = anc),
    safe_coef_ci(sm[keep, 1], sm[keep, 2])
  ) %>%
    mutate(
      excludes_zero = lo > 0 | hi < 0,
      note = note
    )
}

extract_ancestry_lmer <- function(fit, model_label, cohort_def, note = NA_character_) {
  sm <- summary(fit)$coefficients
  rn <- rownames(sm)
  keep <- grepl("^ancestry", rn)
  if (!any(keep)) {
    return(tibble(
      model = model_label, cohort_def = cohort_def, ancestry = character(),
      estimate = numeric(), se = numeric(), lo = numeric(), hi = numeric(),
      excludes_zero = logical(), note = character()
    ))
  }
  anc <- sub("^ancestry", "", rn[keep])
  bind_cols(
    tibble(model = model_label, cohort_def = cohort_def, ancestry = anc),
    safe_coef_ci(sm[keep, 1], sm[keep, 2])
  ) %>%
    mutate(
      excludes_zero = lo > 0 | hi < 0,
      note = note
    )
}

extract_ancestry_rma <- function(fit, model_label, cohort_def, note = NA_character_) {
  if (inherits(fit, "try-error") || is.null(fit)) {
    return(tibble(
      model = model_label, cohort_def = cohort_def, ancestry = character(),
      estimate = numeric(), se = numeric(), lo = numeric(), hi = numeric(),
      excludes_zero = logical(), note = note
    ))
  }
  b <- as.numeric(fit$beta)
  se <- as.numeric(fit$se)
  rn <- rownames(fit$beta)
  if (is.null(rn)) rn <- names(b)
  keep <- grepl("^ancestry", rn)
  if (!any(keep)) {
    return(tibble(
      model = model_label, cohort_def = cohort_def, ancestry = character(),
      estimate = numeric(), se = numeric(), lo = numeric(), hi = numeric(),
      excludes_zero = logical(), note = note
    ))
  }
  anc <- sub("^ancestry", "", rn[keep])
  tibble(
    model = model_label,
    cohort_def = cohort_def,
    ancestry = anc,
    estimate = b[keep],
    se = se[keep],
    lo = fit$ci.lb[keep],
    hi = fit$ci.ub[keep],
    excludes_zero = fit$ci.lb[keep] > 0 | fit$ci.ub[keep] < 0,
    note = note
  )
}

lmer_note <- function(fit) {
  if (inherits(fit, "try-error")) return(as.character(fit))
  msgs <- character()
  warns <- attr(fit, "fit_warnings")
  conv_msg <- fit@optinfo$conv$lme4$messages
  noisy <- c(warns, conv_msg)
  if (length(noisy) && any(grepl("converge|gradient|Hessian", noisy, ignore.case = TRUE))) {
    msgs <- c(msgs, "lmer convergence warning")
  }
  if (isTRUE(lme4::isSingular(fit, tol = 1e-4))) {
    msgs <- c(msgs, "singular fit")
  }
  vc <- as.data.frame(lme4::VarCorr(fit))
  if (any(vc$vcov < 1e-10, na.rm = TRUE)) {
    msgs <- c(msgs, "variance on boundary")
  }
  if (!length(msgs)) NA_character_ else paste(unique(msgs), collapse = "; ")
}

run_lmer <- function(formula, dat) {
  warns <- character()
  fit <- withCallingHandlers(
    try(
      lmer(
        formula,
        data = dat,
        weights = weight_auc,
        control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))
      ),
      silent = TRUE
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (!inherits(fit, "try-error")) attr(fit, "fit_warnings") <- warns
  fit
}

# ---------------------------------------------------------------------------
# 0) Load inputs
# ---------------------------------------------------------------------------

message(glue::glue(">> Diabetes cohort models on Catalog stamp {stamp}"))

fp_eval  <- file.path("results", paste0("pgs_auc_ci_audit_", stamp), "eval_df_final_auc_ci.csv")
fp_s1f   <- file.path(catalog_results_dir(stamp), "stage1", paste0("stage1_pooled_cells_I2filtered_", stamp, ".csv"))
fp_pss   <- catalog_bulk_file("evaluation_sample_sets", stamp)
fp_ppm   <- catalog_bulk_file("performance_metrics", stamp)
fp_dep   <- file.path("results", paste0("review_sensitivity_", stamp), "tables", "09_dependence_cells.csv")
fp_pairs <- file.path("results", paste0("review_sensitivity_", stamp), "tables", "10_covariates_pairs.csv")
fp_ids   <- file.path("results", paste0("review_sensitivity_", stamp), "tables", "13_cell_identifiers_long.csv")

stopifnot(
  file.exists(fp_eval), file.exists(fp_s1f), file.exists(fp_pss),
  file.exists(fp_ppm), file.exists(fp_dep), file.exists(fp_pairs),
  file.exists(fp_ids)
)

eval_raw <- read_csv(fp_eval, show_col_types = FALSE)
stage1_f <- read_csv(fp_s1f, show_col_types = FALSE)
dep_cells <- read_csv(fp_dep, show_col_types = FALSE) %>%
  filter(trait_label == "Type 2 diabetes")
pairs_t2d <- read_csv(fp_pairs, show_col_types = FALSE) %>%
  filter(trait_label == "Type 2 diabetes")
ids_t2d <- read_csv(fp_ids, show_col_types = FALSE) %>%
  filter(trait_label == "Type 2 diabetes")

pss <- read_csv(fp_pss, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
  clean_names() %>%
  transmute(
    sampleset_id  = as.character(pgs_sample_set_pss),
    n_individuals = suppressWarnings(as.numeric(number_of_individuals)),
    n_cases       = suppressWarnings(as.numeric(number_of_cases)),
    n_controls    = suppressWarnings(as.numeric(number_of_controls)),
    cohort_txt    = if ("cohort_s" %in% names(.)) as.character(cohort_s) else NA_character_
  ) %>%
  distinct(sampleset_id, .keep_all = TRUE)

ppm_cov <- read_csv(fp_ppm, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
  clean_names() %>%
  transmute(
    performance_id = as.character(pgs_performance_metric_ppm_id),
    covariates = covariates_included_in_the_model
  ) %>%
  distinct(performance_id, .keep_all = TRUE)

# ---------------------------------------------------------------------------
# 1) Diabetes evaluation table: same filters Stage 1 uses, then keep
#    scores that survive I2 filtering and have a European + target pair.
# ---------------------------------------------------------------------------

eval_proc <- eval_raw %>%
  mutate(
    pgs_id = as.character(pgs_id),
    sampleset_id = as.character(sampleset_id),
    performance_id = as.character(performance_id),
    reported_trait = as.character(reported_trait),
    train_bucket = as.character(train_bucket),
    ancestry_eval = as.character(ancestry_eval),
    raw_label = dplyr::coalesce(na_if(trimws(reported_trait), ""), NA_character_),
    trait_label = harmonize_trait_label(NA_character_, raw_label),
    ancestry_display = ancestry_eval
  ) %>%
  mutate(across(c(auc, estimate_ci_lower, estimate_ci_upper), as.numeric)) %>%
  filter(ancestry_display %in% evaluation_display_levels) %>%
  filter(is.finite(auc), is.finite(estimate_ci_lower), is.finite(estimate_ci_upper)) %>%
  left_join(
    pss %>% select(sampleset_id, n = n_individuals, cases = n_cases, controls = n_controls),
    by = "sampleset_id"
  ) %>%
  left_join(ppm_cov, by = "performance_id") %>%
  mutate(covariate_class = classify_covariates(covariates))

interval_issue <- auc_interval_issue(eval_proc$estimate_ci_lower, eval_proc$estimate_ci_upper)
if (any(!is.na(interval_issue))) {
  stop("Bad AUC intervals reached the diabetes table.", call. = FALSE)
}

one_eval <- keep_one_evaluation(
  eval_proc,
  class_col = "covariate_class",
  return_dropped = FALSE
)

diabetes_stage1_cells <- stage1_f %>%
  filter(trait_label == "Type 2 diabetes") %>%
  transmute(
    pgs_id = as.character(pgs_id),
    trait_label,
    ancestry_display = as.character(ancestry_display),
    covariate_class = as.character(covariate_class),
    eta_stage1 = eta,
    se_stage1 = se,
    auc_pooled_stage1 = auc_pooled,
    k_eval_stage1 = k_eval
  )

keep_bucket_pattern <- c("European-only", "Multi incl")

diabetes_rows <- one_eval %>%
  filter(trait_label == "Type 2 diabetes") %>%
  filter(grepl(paste(keep_bucket_pattern, collapse = "|"), train_bucket, ignore.case = TRUE)) %>%
  inner_join(
    diabetes_stage1_cells %>% select(pgs_id, trait_label, ancestry_display, covariate_class),
    by = c("pgs_id", "trait_label", "ancestry_display", "covariate_class")
  ) %>%
  left_join(pss %>% select(sampleset_id, cohort_txt, n_individuals), by = "sampleset_id") %>%
  mutate(
    tmp = pmap(list(auc, estimate_ci_lower, estimate_ci_upper), auc_ci_to_logit)
  ) %>%
  unnest_wider(tmp) %>%
  filter(is.finite(eta), is.finite(se), se > 0) %>%
  mutate(
    se_auc = auc_scale_se(auc, se, "delta_method"),
    study_name = normalize_study(cohort_txt),
    weight_auc = 1 / (se_auc^2),
    weight_logit = 1 / (se^2)
  )

# Scores that have European and at least one target ancestry in the same bucket.
pair_keys <- diabetes_rows %>%
  group_by(train_bucket, pgs_id) %>%
  summarise(
    has_eur = any(ancestry_display == "European"),
    has_tgt = any(ancestry_display != "European"),
    .groups = "drop"
  ) %>%
  filter(has_eur, has_tgt)

diabetes_eval <- diabetes_rows %>%
  inner_join(pair_keys %>% select(train_bucket, pgs_id), by = c("train_bucket", "pgs_id")) %>%
  arrange(train_bucket, ancestry_display, pgs_id, sampleset_id, performance_id)

# Stage-1 logit check: within-cell IVW of the evaluation logits must match.
stage1_check <- diabetes_eval %>%
  group_by(pgs_id, trait_label, ancestry_display, covariate_class) %>%
  summarise(
    pooled = list(ivw_pool_logit(eta, se)),
    .groups = "drop"
  ) %>%
  unnest_wider(pooled) %>%
  inner_join(
    diabetes_stage1_cells,
    by = c("pgs_id", "trait_label", "ancestry_display", "covariate_class")
  ) %>%
  mutate(
    abs_eta_diff = abs(eta - eta_stage1),
    abs_se_diff = abs(se - se_stage1)
  )

bad_check <- stage1_check %>% filter(abs_eta_diff > 1e-8 | abs_se_diff > 1e-8)
if (nrow(bad_check) > 0) {
  print(bad_check %>% select(pgs_id, ancestry_display, covariate_class, eta, eta_stage1, se, se_stage1))
  stop("Diabetes evaluation logits do not match Stage-1 I2-filtered cells.", call. = FALSE)
}
message(glue::glue(
  ">> Evaluation table: {nrow(diabetes_eval)} rows, ",
  "{n_distinct(diabetes_eval$pgs_id)} scores, ",
  "{n_distinct(diabetes_eval$sampleset_id)} sample sets. Stage-1 check OK."
))

write_csv(
  diabetes_eval %>%
    select(
      performance_id, pgs_id, sampleset_id, cohort_txt, study_name,
      ancestry_display, train_bucket, trait_label,
      auc, estimate_ci_lower, estimate_ci_upper,
      eta, se, se_auc, n_individuals, covariates, covariate_class
    ),
  file.path(out_dir, "diabetes_evaluation_rows.csv")
)

# ---------------------------------------------------------------------------
# 2) Verify the South Asian DerSimonian-Laird numbers from the notes
# ---------------------------------------------------------------------------

sa_pairs <- pairs_t2d %>%
  filter(
    trained_bucket == "Train: European-only",
    target_ancestry == "South Asian"
  )

sa_meta <- ids_t2d %>%
  filter(
    trained_bucket == "Train: European-only",
    target_ancestry == "South Asian",
    side == "Target"
  ) %>%
  group_by(pgs_id, sampleset_id) %>%
  summarise(
    n_individuals = max(n_individuals, na.rm = TRUE),
    cohort_txt = paste(unique(trimws(cohort_txt)), collapse = "+"),
    .groups = "drop"
  ) %>%
  group_by(pgs_id) %>%
  summarise(
    tgt_key = paste(sort(unique(sampleset_id)), collapse = ";"),
    .groups = "drop"
  )

sa_kept <- sa_pairs %>%
  left_join(sa_meta, by = "pgs_id") %>%
  group_by(tgt_key) %>%
  arrange(pgs_id, .by_group = TRUE) %>%
  mutate(kept = row_number() == 1L) %>%
  ungroup()

stopifnot(nrow(sa_kept) == 13L, sum(sa_kept$kept) == 3L)
sa_all <- pool_delta_check(sa_kept$delta, sa_kept$delta_se, "fixed")
sa_one <- pool_delta_check(sa_kept$delta[sa_kept$kept], sa_kept$delta_se[sa_kept$kept], "fixed")
sa_re  <- pool_delta_check(sa_kept$delta[sa_kept$kept], sa_kept$delta_se[sa_kept$kept], "random_dl")

sa_pub <- dep_cells %>%
  filter(trained_bucket == "Train: European-only", target_ancestry == "South Asian")
stopifnot(nrow(sa_pub) == 1L)
stopifnot(abs(sa_all$delta_hat - sa_pub$delta_fixed) < 1e-8)
stopifnot(abs(sa_one$delta_hat - sa_pub$delta_one_per_pss) < 1e-8)
stopifnot(abs(sa_re$delta_hat - sa_pub$delta_random) < 1e-8)
stopifnot(abs(sa_re$I2 - sa_pub$I2_stage2) < 1e-6)

sa_check <- tibble(
  estimate = c("IVW, all 13 scores", "IVW, 3 independent scores", "Random effects (DL), 3 independent scores"),
  delta_hat = c(sa_all$delta_hat, sa_one$delta_hat, sa_re$delta_hat),
  lo = c(sa_all$lo, sa_one$lo, sa_re$lo),
  hi = c(sa_all$hi, sa_one$hi, sa_re$hi),
  Q = c(sa_all$Q, sa_one$Q, sa_re$Q),
  I2 = c(sa_all$I2, sa_one$I2, sa_re$I2),
  k = c(sa_all$k, sa_one$k, sa_re$k),
  note = "Current DerSimonian-Laird on three DeltaAUC values. Not one of the new ancestry models."
)
write_csv(sa_check, file.path(out_dir, "south_asian_dl_check.csv"))
message(glue::glue(
  ">> South Asian check: IVW13={fmt4(sa_all$delta_hat)}, IVW3={fmt4(sa_one$delta_hat)}, ",
  "RE={fmt4(sa_re$delta_hat)}, Q={fmt3(sa_re$Q)}, I2={fmt3(sa_re$I2)}"
))

# ---------------------------------------------------------------------------
# 3) Fit models within each training bucket
# ---------------------------------------------------------------------------

fit_bucket <- function(df, cohort_col, cohort_label, include_lm = FALSE) {
  d <- df %>%
    mutate(
      ancestry = factor(ancestry_display, levels = c(
        "European",
        setdiff(evaluation_display_levels, "European")
      )),
      PGS = factor(pgs_id),
      cohort = factor(.data[[cohort_col]])
    ) %>%
    filter(is.finite(auc), is.finite(se_auc), se_auc > 0, is.finite(eta), is.finite(se), se > 0)

  n_pgs <- n_distinct(d$PGS)
  n_cohort <- n_distinct(d$cohort)
  n_rows <- nrow(d)

  out <- list()

  if (isTRUE(include_lm)) {
    # lm models have no cohort term; fit once per bucket.
    fit1 <- lm(auc ~ ancestry + PGS, data = d)
    out[[length(out) + 1]] <- extract_ancestry_lm(
      fit1, "lm(AUC ~ ancestry + PGS)", "— (not used)",
      note = glue::glue("{n_rows} rows, {n_pgs} scores")
    )
    fit2 <- lm(auc ~ ancestry + PGS, data = d, weights = weight_auc)
    out[[length(out) + 1]] <- extract_ancestry_lm(
      fit2, "lm(AUC ~ ancestry + PGS, weights = 1/SE^2)", "— (not used)",
      note = "Joint weighted regression. Not the cell-wise Figure 1C pool."
    )
  }

  # 3) lmer with fixed PGS + random cohort.
  # Pass the weight vector, not the column name: lmer would otherwise look
  # the name up in the formula environment.
  fit3 <- run_lmer(auc ~ ancestry + PGS + (1 | cohort), d)
  n3 <- if (inherits(fit3, "try-error")) as.character(fit3) else lmer_note(fit3)
  if (!inherits(fit3, "try-error")) {
    out[[length(out) + 1]] <- extract_ancestry_lmer(
      fit3, "lmer(AUC ~ ancestry + PGS + (1|cohort), weights = 1/SE^2)", cohort_label, note = n3
    )
  } else {
    out[[length(out) + 1]] <- tibble(
      model = "lmer(AUC ~ ancestry + PGS + (1|cohort), weights = 1/SE^2)",
      cohort_def = cohort_label, ancestry = NA_character_,
      estimate = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_,
      excludes_zero = NA, note = n3
    )
  }

  # 4) lmer with random PGS + random cohort
  fit4 <- run_lmer(auc ~ ancestry + (1 | PGS) + (1 | cohort), d)
  if (!inherits(fit4, "try-error")) {
    note4_parts <- c(glue::glue("{n_pgs} scores as random intercepts"), lmer_note(fit4))
    note4 <- paste(note4_parts[!is.na(note4_parts) & nzchar(note4_parts)], collapse = "; ")
    out[[length(out) + 1]] <- extract_ancestry_lmer(
      fit4, "lmer(AUC ~ ancestry + (1|PGS) + (1|cohort), weights = 1/SE^2)", cohort_label, note = note4
    )
  } else {
    out[[length(out) + 1]] <- tibble(
      model = "lmer(AUC ~ ancestry + (1|PGS) + (1|cohort), weights = 1/SE^2)",
      cohort_def = cohort_label, ancestry = NA_character_,
      estimate = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_,
      excludes_zero = NA, note = as.character(fit4)
    )
  }

  # metafor cross-check: fixed PGS + random cohort on AUC
  fit5 <- try(
    rma.mv(
      yi = d$auc, V = d$se_auc^2, mods = ~ ancestry + PGS,
      random = ~ 1 | cohort, data = d, sparse = TRUE, verbose = FALSE
    ),
    silent = TRUE
  )
  n5 <- if (inherits(fit5, "try-error")) as.character(fit5) else {
    msgs <- character()
    if (!is.null(fit5$sigma2) && any(fit5$sigma2 < 1e-10, na.rm = TRUE)) {
      msgs <- c(msgs, "variance on boundary")
    }
    if (!length(msgs)) "V = SE_AUC^2 fixed" else paste(c("V = SE_AUC^2 fixed", msgs), collapse = "; ")
  }
  out[[length(out) + 1]] <- extract_ancestry_rma(
    fit5, "rma.mv(AUC ~ ancestry + PGS + (1|cohort), V = SE^2)", cohort_label, note = n5
  )

  # metafor: random PGS + random cohort on AUC
  fit6 <- try(
    rma.mv(
      yi = d$auc, V = d$se_auc^2, mods = ~ ancestry,
      random = list(~ 1 | PGS, ~ 1 | cohort), data = d, sparse = TRUE, verbose = FALSE
    ),
    silent = TRUE
  )
  n6 <- if (inherits(fit6, "try-error")) as.character(fit6) else {
    msgs <- c(glue::glue("{n_pgs} scores as random intercepts"), "V = SE_AUC^2 fixed")
    if (!is.null(fit6$sigma2) && any(fit6$sigma2 < 1e-10, na.rm = TRUE)) {
      msgs <- c(msgs, "variance on boundary")
    }
    paste(msgs, collapse = "; ")
  }
  out[[length(out) + 1]] <- extract_ancestry_rma(
    fit6, "rma.mv(AUC ~ ancestry + (1|PGS) + (1|cohort), V = SE^2)", cohort_label, note = n6
  )

  # metafor on logit: fixed PGS + random cohort
  fit7 <- try(
    rma.mv(
      yi = d$eta, V = d$se^2, mods = ~ ancestry + PGS,
      random = ~ 1 | cohort, data = d, sparse = TRUE, verbose = FALSE
    ),
    silent = TRUE
  )
  n7 <- if (inherits(fit7, "try-error")) as.character(fit7) else {
    msgs <- "V = SE_logit^2 fixed; outcome is logit(AUC)"
    if (!is.null(fit7$sigma2) && any(fit7$sigma2 < 1e-10, na.rm = TRUE)) {
      msgs <- paste(msgs, "variance on boundary", sep = "; ")
    }
    msgs
  }
  out[[length(out) + 1]] <- extract_ancestry_rma(
    fit7, "rma.mv(logit(AUC) ~ ancestry + PGS + (1|cohort), V = SE^2)", cohort_label, note = n7
  )

  # metafor on logit: random PGS + random cohort
  fit8 <- try(
    rma.mv(
      yi = d$eta, V = d$se^2, mods = ~ ancestry,
      random = list(~ 1 | PGS, ~ 1 | cohort), data = d, sparse = TRUE, verbose = FALSE
    ),
    silent = TRUE
  )
  n8 <- if (inherits(fit8, "try-error")) as.character(fit8) else {
    msgs <- c(
      glue::glue("{n_pgs} scores as random intercepts"),
      "V = SE_logit^2 fixed; outcome is logit(AUC)"
    )
    if (!is.null(fit8$sigma2) && any(fit8$sigma2 < 1e-10, na.rm = TRUE)) {
      msgs <- c(msgs, "variance on boundary")
    }
    paste(msgs, collapse = "; ")
  }
  out[[length(out) + 1]] <- extract_ancestry_rma(
    fit8, "rma.mv(logit(AUC) ~ ancestry + (1|PGS) + (1|cohort), V = SE^2)", cohort_label, note = n8
  )

  bind_rows(out) %>%
    mutate(
      train_bucket = unique(df$train_bucket),
      n_rows = n_rows,
      n_pgs = n_pgs,
      n_cohorts = n_cohort
    )
}

coef_rows <- list()
for (bucket in sort(unique(diabetes_eval$train_bucket))) {
  d_b <- diabetes_eval %>% filter(train_bucket == bucket)
  message(glue::glue(">> Fitting models for {bucket} ({nrow(d_b)} rows)"))
  coef_rows[[length(coef_rows) + 1]] <- fit_bucket(
    d_b, "sampleset_id", "sample-set id", include_lm = TRUE
  )
  coef_rows[[length(coef_rows) + 1]] <- fit_bucket(
    d_b, "study_name", "study name", include_lm = FALSE
  )
}
coef_tbl <- bind_rows(coef_rows) %>%
  relocate(train_bucket, model, cohort_def, ancestry)

write_csv(coef_tbl, file.path(out_dir, "diabetes_model_coefficients.csv"))

# ---------------------------------------------------------------------------
# 4) Reduced-analysis row counts (for Georg's last question)
# ---------------------------------------------------------------------------

# One score per sample set (lowest PGS id within sample set × ancestry), as in script 09.
one_per_pss <- diabetes_eval %>%
  group_by(train_bucket, ancestry_display, sampleset_id) %>%
  arrange(pgs_id, .by_group = TRUE) %>%
  mutate(keep = row_number() == 1L) %>%
  ungroup()

# Largest cohort when one score is evaluated more than once within ancestry.
# Keep the row with the largest n_individuals within (pgs, ancestry); ties by sampleset_id.
largest_cohort <- diabetes_eval %>%
  group_by(train_bucket, ancestry_display, pgs_id, covariate_class) %>%
  arrange(desc(n_individuals), sampleset_id, .by_group = TRUE) %>%
  mutate(keep = row_number() == 1L) %>%
  ungroup()

drop_counts <- tibble(
  rule = c(
    "One score per sample set (lowest PGS ID within sample set × ancestry)",
    "Largest cohort when one score is evaluated more than once within ancestry"
  ),
  n_rows_kept = c(sum(one_per_pss$keep), sum(largest_cohort$keep)),
  n_rows_dropped = c(sum(!one_per_pss$keep), sum(!largest_cohort$keep)),
  n_rows_total = nrow(diabetes_eval)
)
write_csv(drop_counts, file.path(out_dir, "diabetes_reduced_rule_counts.csv"))

# ---------------------------------------------------------------------------
# 5) HTML for Georg
# ---------------------------------------------------------------------------

model_blurbs <- tribble(
  ~model_key, ~title, ~blurb, ~question,
  "lm_unweighted",
  "Unweighted regression: lm(AUC ~ ancestry + PGS)",
  "Each PGS has its own baseline. Every evaluation has the same weight. A shared sample set is not marked.",
  "Should the unweighted regression be the dependence model?",
  "lm_weighted",
  "Weighted regression: lm(AUC ~ ancestry + PGS, weights = 1/SE^2)",
  "Each PGS has its own baseline, and a precise AUC pulls harder. A shared sample set is not marked. For some ancestry groups the point matches All PGS IVW. The interval is the regression interval, not the All PGS IVW interval.",
  "Should the weighted regression be the dependence model?",
  "lmer_pgs_fixed",
  "Sample-set term, PGS fixed: lmer(AUC ~ ancestry + PGS + (1|cohort), weights = 1/SE^2)",
  "In this model, cohort is the sample set. Evaluations from the same sample set share an offset, and the PGS stays a fixed intercept. Eleven South Asian PGS share PSS011749, and all of them are South Asian, so the sample-set offset and the ancestry coefficient compete for the same shift.",
  "Should the model with a sample-set term and a fixed PGS be the dependence model?",
  "lmer_pgs_random",
  "PGS and sample set both random: lmer(AUC ~ ancestry + (1|PGS) + (1|cohort), weights = 1/SE^2)",
  "Cohort is again the sample set. The PGS is a random intercept as well. European-only training has 13 PGS. Multi-ancestry training has 6.",
  "Should the model with both the PGS and the sample set random be the dependence model?"
)

map_model_key <- function(model) {
  case_when(
    grepl("^lm\\(AUC ~ ancestry \\+ PGS\\)$", model) ~ "lm_unweighted",
    grepl("^lm\\(AUC ~ ancestry \\+ PGS, weights", model) ~ "lm_weighted",
    grepl("lmer\\(AUC ~ ancestry \\+ PGS \\+ \\(1\\|cohort\\)", model) ~ "lmer_pgs_fixed",
    grepl("lmer\\(AUC ~ ancestry \\+ \\(1\\|PGS\\) \\+ \\(1\\|cohort\\)", model) ~ "lmer_pgs_random",
    grepl("^rma\\.mv\\(AUC", model) ~ "rma_auc",
    grepl("^rma\\.mv\\(logit", model) ~ "rma_logit",
    TRUE ~ "other"
  )
}

coef_view <- coef_tbl %>%
  mutate(model_key = map_model_key(model)) %>%
  filter(!is.na(ancestry))

# The page shows Nikolaus's four models. The cohort in those mixed models is
# the sample-set id. Study-name fits and the metafor weight check stay in the CSV.
page_rows <- coef_view %>%
  filter(
    model_key %in% model_blurbs$model_key,
    cohort_def %in% c("— (not used)", "sample-set id")
  )

cohort_support <- diabetes_eval %>%
  group_by(train_bucket, sampleset_id) %>%
  summarise(n_in_set = n(), .groups = "drop") %>%
  group_by(train_bucket) %>%
  summarise(
    n_rows_bucket = sum(n_in_set),
    n_sets = n(),
    n_sets_repeated = sum(n_in_set > 1),
    .groups = "drop"
  )

page_note <- function(r) {
  bits <- character()
  if (identical(r$model_key, "lm_unweighted")) {
    bits <- c(bits, "A shared sample set is not marked.")
  }
  if (identical(r$model_key, "lm_weighted")) {
    ref <- ref_rows %>%
      filter(train_bucket == r$train_bucket, ancestry == r$ancestry)
    if (nrow(ref) == 1L && abs(r$estimate - ref$estimate) < 0.002) {
      bits <- c(bits, sprintf(
        "Point matches All PGS IVW (%s). The interval is the regression interval, not the All PGS IVW interval. A shared sample set is not marked.",
        fmt4(ref$estimate)
      ))
    } else if (nrow(ref) == 1L) {
      bits <- c(bits, sprintf(
        "Point is %s; All PGS IVW is %s. This joint regression is not the All PGS IVW pool for that comparison. A shared sample set is not marked.",
        fmt4(r$estimate), fmt4(ref$estimate)
      ))
    }
  }
  if (r$model_key %in% c("lmer_pgs_fixed", "lmer_pgs_random")) {
    if (isTRUE(grepl("convergence", r$note))) {
      bits <- c(bits, "lmer convergence warning.")
    }
    support <- cohort_support %>% filter(train_bucket == r$train_bucket)
    if (nrow(support) == 1L && grepl("Multi", r$train_bucket)) {
      bits <- c(bits, sprintf(
        "%d evaluations and %d sample sets; %d sample sets contain more than one PGS. The sample-set term is barely estimated.",
        support$n_rows_bucket, support$n_sets, support$n_sets_repeated
      ))
    }
    if (identical(r$model_key, "lmer_pgs_random") && !is.na(r$n_pgs)) {
      bits <- c(bits, sprintf("%d PGS as random intercepts.", r$n_pgs))
    }
  }
  if (!length(bits)) "" else paste(bits, collapse = " ")
}

ref_rows <- dep_cells %>%
  transmute(
    train_bucket = trained_bucket,
    ancestry = target_ancestry,
    estimate = delta_fixed,
    lo = lo_fixed,
    hi = hi_fixed,
    excludes_zero = sig_fixed,
    n_pairs,
    source = "All PGS IVW"
  )

# Cohort composition for European-only South Asian / Hispanic notes
sa_cohorts <- diabetes_eval %>%
  filter(
    train_bucket == "Train: European-only",
    ancestry_display == "South Asian"
  ) %>%
  count(sampleset_id, study_name, cohort_txt, sort = TRUE)

hisp_cohorts <- diabetes_eval %>%
  filter(
    train_bucket == "Train: European-only",
    ancestry_display == "Hispanic or Latin American"
  ) %>%
  count(sampleset_id, study_name, cohort_txt, sort = TRUE)

eas_ref <- ref_rows %>%
  filter(train_bucket == "Train: European-only", ancestry == "East Asian")
eas_study <- coef_tbl %>%
  filter(
    train_bucket == "Train: European-only",
    ancestry == "East Asian",
    cohort_def == "study name",
    grepl("^lmer\\(AUC ~ ancestry \\+ PGS \\+", model)
  )
sas_lmer <- coef_tbl %>%
  filter(
    train_bucket == "Train: European-only",
    ancestry == "South Asian",
    cohort_def == "sample-set id",
    grepl("^lmer\\(AUC ~ ancestry \\+ PGS \\+", model)
  )
sas_rma <- coef_tbl %>%
  filter(
    train_bucket == "Train: European-only",
    ancestry == "South Asian",
    cohort_def == "sample-set id",
    grepl("^rma\\.mv\\(AUC ~ ancestry \\+ PGS \\+", model)
  )

cell_table_html <- function(bucket, ancestry) {
  ref <- ref_rows %>% filter(train_bucket == bucket, ancestry == .env$ancestry)
  rows <- page_rows %>%
    filter(train_bucket == bucket, ancestry == .env$ancestry) %>%
    arrange(model_key, cohort_def, model)

  lines <- c(
    "<table>",
    "<thead><tr>",
    "<th>Source</th><th>Sample set in the model</th><th>Estimate</th><th>95% CI</th><th>Excludes 0</th><th>Note</th>",
    "</tr></thead><tbody>"
  )
  if (nrow(ref) == 1L) {
    lines <- c(lines, sprintf(
      "<tr class='ref'><td>%s</td><td>—</td><td>%s</td><td>%s to %s</td><td>%s</td><td></td></tr>",
      html_escape(ref$source),
      fmt4(ref$estimate), fmt4(ref$lo), fmt4(ref$hi),
      ifelse(isTRUE(ref$excludes_zero), "yes", "no")
    ))
  } else {
    lines <- c(lines, "<tr class='ref'><td colspan='6'>No All PGS IVW row for this comparison.</td></tr>")
  }
  if (nrow(rows) == 0L) {
    lines <- c(lines, "<tr><td colspan='6'>No model estimate for this ancestry in this bucket.</td></tr>")
  } else {
    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      lines <- c(lines, sprintf(
        "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s to %s</td><td>%s</td><td>%s</td></tr>",
        html_escape(r$model),
        html_escape(r$cohort_def),
        fmt4(r$estimate), fmt4(r$lo), fmt4(r$hi),
        ifelse(isTRUE(r$excludes_zero), "yes", ifelse(isFALSE(r$excludes_zero), "no", "NA")),
        html_escape(page_note(r))
      ))
    }
  }
  c(lines, "</tbody></table>")
}

fig_dir <- file.path(
  "results", paste0("review_sensitivity_", stamp), "q_random_effects_t2d_illustrations"
)
embed_png <- function(filename, alt) {
  path <- file.path(fig_dir, filename)
  if (!file.exists(path)) {
    stop(
      "Missing ", path, ". Run tests/q_random_effects_t2d_illustrations.R first.",
      call. = FALSE
    )
  }
  raw <- readBin(path, "raw", file.info(path)$size)
  b64 <- gsub("[\r\n]", "", jsonlite::base64_enc(raw))
  sprintf(
    "<img class='fig' alt='%s' src='data:image/png;base64,%s'>",
    html_escape(alt), b64
  )
}

buckets <- sort(unique(diabetes_eval$train_bucket))
target_ancestries_by_bucket <- lapply(buckets, function(b) {
  sort(setdiff(unique(diabetes_eval$ancestry_display[diabetes_eval$train_bucket == b]), "European"))
})
names(target_ancestries_by_bucket) <- buckets

css <- "
body { font-family: DejaVu Sans, Helvetica, Arial, sans-serif; margin: 24px; color: #1a1a1a; line-height: 1.45; max-width: 1100px; }
h1 { font-size: 1.6rem; margin-bottom: 0.4rem; }
h2 { font-size: 1.25rem; margin-top: 2rem; border-bottom: 1px solid #ddd; padding-bottom: 0.25rem; }
h3 { font-size: 1.05rem; margin-top: 1.4rem; }
.lead { background: #f4f7fb; border-left: 4px solid #332288; padding: 12px 14px; margin: 12px 0 20px; }
.box { background: #fafafa; border: 1px solid #ddd; padding: 12px 14px; margin: 12px 0; }
.q { background: #fff8e8; border-left: 4px solid #E69F00; padding: 10px 12px; margin: 10px 0 18px; font-weight: 600; }
table { border-collapse: collapse; width: 100%; font-size: 0.86rem; margin: 10px 0 18px; }
th, td { border: 1px solid #ccc; padding: 6px 8px; text-align: left; vertical-align: top; }
th { background: #f0f0f0; }
tr.ref { background: #eef6ee; font-weight: 600; }
code { font-family: DejaVu Sans Mono, Menlo, monospace; font-size: 0.9em; }
ul { margin-top: 0.4rem; }
.small { color: #444; font-size: 0.9rem; }
img.fig { max-width: 100%; height: auto; border: 1px solid #ddd; margin: 8px 0 4px; }
.cap { color: #444; font-size: 0.9rem; margin: 0 0 22px; }
"

html <- c(
  "<!DOCTYPE html>",
  "<html lang='en'><head><meta charset='utf-8'>",
  "<title>Type 2 diabetes — which model for PGS that share a sample set</title>",
  paste0("<style>", css, "</style>"),
  "</head><body>",
  "<h1>Type 2 diabetes: PGS that share a sample set</h1>",
  paste0(
    "<p class='small'>Catalog stamp ", html_escape(stamp),
    ". In the models, cohort means the sample set.</p>"
  ),

  "<div class='q'>",
  "The input I need: which of the four models should account for several PGS evaluated on the same sample set, ",
  "or whether you would rather drop rows than model the sharing.",
  "</div>",
  "<ul>",
  "<li><code>lm(AUC ~ ancestry + PGS)</code></li>",
  "<li><code>lm(AUC ~ ancestry + PGS, weights = 1/SE^2)</code></li>",
  "<li><code>lmer(AUC ~ ancestry + PGS + (1|cohort), weights = 1/SE^2)</code></li>",
  "<li><code>lmer(AUC ~ ancestry + (1|PGS) + (1|cohort), weights = 1/SE^2)</code></li>",
  "</ul>",

  "<div class='lead'>",
  "<p>An evaluation is one AUC: one PGS, one sample set, one ancestry group, and one standard error. ",
  "The sample-set term marks which evaluations share a sample set. It does not make the AUCs independent.</p>",
  "<p>The four models are fit on type 2 diabetes. They can separate an ancestry group from a sample set only when a comparison uses more than one target sample set. ",
  "Nine comparisons do. Seven are type 2 diabetes. The other two are breast cancer, European-only training: South Asian participants, and East Asian participants.</p>",
  "<p>Eighteen comparisons have several PGS and one target sample set. Ancestry group and sample set are the same column there. ",
  "The four models are not for that case. Five comparisons have one PGS.</p>",
  "</div>",

  "<h2>South Asian comparison</h2>",
  "<div class='box'>",
  "<p>Europeans and South Asian participants, European-only training. Thirteen PGS, eleven of them on one UK Biobank sample set.</p>",
  "<ul>",
  sprintf("<li>Inverse-variance pool of all 13: <strong>%s</strong> (%s to %s).</li>",
          fmt4(sa_all$delta_hat), fmt4(sa_all$lo), fmt4(sa_all$hi)),
  sprintf("<li>Pool of the 3 kept PGS, one per sample set: <strong>%s</strong> (%s to %s).</li>",
          fmt4(sa_one$delta_hat), fmt4(sa_one$lo), fmt4(sa_one$hi)),
  "</ul>",
  "</div>",

  "<h2>Two comparisons</h2>",
  "<p>South Asian participants, below, are one of the nine. ",
  "Hispanic or Latin American participants are one of the eighteen, shown so the boundary is visible. ",
  "The third figure lists the nine comparisons the question applies to.</p>",

  embed_png("02_south_asian_one_of_eight.png", "Type 2 diabetes, Europeans and South Asian participants, European-only training"),
  "<p class='cap'>Thirteen PGS. Eleven share PSS011749. Two other target sample sets are in the comparison. ",
  "The black line is All PGS IVW, the pool of all thirteen, about 0. One PGS per target sample set is −0.016 and excludes 0.</p>",

  embed_png("03_hispanic_one_of_seventeen.png", "Type 2 diabetes, Europeans and Hispanic or Latin American participants, European-only training"),
  "<p class='cap'>Eleven PGS, all on PSS011736. Ancestry group and sample set are the same column. ",
  "This comparison is not one of the four models.</p>",

  embed_png("05_eight_one_score_per_cohort.png", "Nine comparisons: every PGS, and one PGS per target sample set"),
  "<p class='cap'>Black is All PGS IVW, the same points as in the tables below. Purple is one PGS per target sample set. ",
  "The four models below are for these nine comparisons.</p>",

  "<h2>Evaluation table used here</h2>",
  "<div class='box'>",
  sprintf(
    "<p>%d evaluations. %d PGS. %d sample sets. %d study names. ",
    nrow(diabetes_eval),
    n_distinct(diabetes_eval$pgs_id),
    n_distinct(diabetes_eval$sampleset_id),
    n_distinct(diabetes_eval$study_name)
  ),
  "Corrected logit values from <code>auc_ci_to_logit()</code>. ",
  "AUC models use the delta-method SE (<code>SE_logit × AUC × (1 − AUC)</code>). ",
  "Fits are separate for each training set. European is the reference ancestry group.</p>",
  "<p class='small'>South Asian sample sets in the European-only training set:</p><ul>",
  paste(sprintf(
    "<li>%s — %s (%s), %d evaluation(s)</li>",
    html_escape(sa_cohorts$sampleset_id),
    html_escape(sa_cohorts$study_name),
    html_escape(sa_cohorts$cohort_txt),
    sa_cohorts$n
  ), collapse = "\n"),
  "</ul>",
  "<p class='small'>CSV: <code>diabetes_evaluation_rows.csv</code>.</p>",
  "</div>"
)

# Per-bucket comparison tables
bucket_title <- function(bucket) {
  if (grepl("European-only", bucket)) "European-only training" else "Multi-ancestry training"
}
for (bucket in buckets) {
  html <- c(html, sprintf("<h2>%s</h2>", html_escape(bucket_title(bucket))))
  for (anc in target_ancestries_by_bucket[[bucket]]) {
    html <- c(
      html,
      sprintf("<h3>%s</h3>", html_escape(anc)),
      cell_table_html(bucket, anc)
    )
  }
}

# Questions per model. Same four fits as the tables above.
html <- c(html, "<h2>Questions for Georg — one per model</h2>")
for (i in seq_len(nrow(model_blurbs))) {
  mb <- model_blurbs[i, ]
  sub <- page_rows %>% filter(model_key == mb$model_key)
  html <- c(
    html,
    sprintf("<h3>%s</h3>", html_escape(mb$title)),
    sprintf("<p>%s</p>", html_escape(mb$blurb)),
    "<p class='small'>Ancestry estimates:</p>",
    "<table><thead><tr><th>Bucket</th><th>Cohort</th><th>Ancestry</th><th>Estimate</th><th>95% CI</th><th>Excludes 0</th><th>Note</th></tr></thead><tbody>"
  )
  if (nrow(sub) == 0L) {
    html <- c(html, "<tr><td colspan='7'>No fit returned ancestry coefficients.</td></tr>")
  } else {
    sub <- sub %>% arrange(train_bucket, cohort_def, ancestry)
    for (j in seq_len(nrow(sub))) {
      r <- sub[j, ]
      html <- c(html, sprintf(
        "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s to %s</td><td>%s</td><td>%s</td></tr>",
        html_escape(bucket_title(r$train_bucket)), html_escape(r$cohort_def), html_escape(r$ancestry),
        fmt4(r$estimate), fmt4(r$lo), fmt4(r$hi),
        ifelse(isTRUE(r$excludes_zero), "yes", ifelse(isFALSE(r$excludes_zero), "no", "NA")),
        html_escape(page_note(r))
      ))
    }
  }
  html <- c(
    html,
    "</tbody></table>",
    sprintf("<div class='q'>%s</div>", html_escape(mb$question))
  )
}

html <- c(
  html,
  "<h2>If dropping rows is preferred to modelling the sharing</h2>",
  "<div class='box'>",
  "<p>Two ways to drop rows, counted on the 71 type 2 diabetes evaluations. The page does not pick one.</p>",
  "<ul>",
  sprintf(
    "<li>One PGS per sample set, the lowest PGS identifier: keeps %d evaluations, drops %d of %d.</li>",
    drop_counts$n_rows_kept[[1]], drop_counts$n_rows_dropped[[1]], drop_counts$n_rows_total[[1]]
  ),
  sprintf(
    "<li>The largest sample set when one PGS is evaluated more than once in an ancestry group: keeps %d evaluations, drops %d of %d.</li>",
    drop_counts$n_rows_kept[[2]], drop_counts$n_rows_dropped[[2]], drop_counts$n_rows_total[[2]]
  ),
  "</ul>",
  "</div>",
  "<div class='q'>If you would rather drop rows than model the sharing, which rule do you prefer: one PGS per sample set, or the largest sample set when one PGS is evaluated more than once?</div>",

  "<h2>Files</h2>",
  "<ul>",
  "<li><code>diabetes_evaluation_rows.csv</code></li>",
  "<li><code>diabetes_model_coefficients.csv</code> — the four models, plus the study-name fits and the metafor weight check</li>",
  "<li><code>diabetes_reduced_rule_counts.csv</code></li>",
  "<li><code>diabetes_cohort_models.html</code> (this page)</li>",
  "</ul>",
  "</body></html>"
)

html_path <- file.path(out_dir, "diabetes_cohort_models.html")
writeLines(html, html_path)
message("Wrote ", normalizePath(html_path))
message("Diabetes cohort models done: ", normalizePath(out_dir))
