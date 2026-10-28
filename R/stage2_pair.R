# Stage-2 paired delta AUC.
# se_scale = "logit_as_implemented" keeps the submitted standard error.
# se_scale = "delta_method" uses se_auc = se_logit * auc * (1 - auc).

auc_scale_se <- function(auc, se_logit, se_scale = c("logit_as_implemented", "delta_method")) {
  se_scale <- match.arg(se_scale)
  if (se_scale == "delta_method") se_logit * auc * (1 - auc) else se_logit
}

stage2_delta_interval <- function(logit_auc, logit_se,
                                  se_scale = c("logit_as_implemented", "delta_method")) {
  se_scale <- match.arg(se_scale)
  if (length(logit_auc) != 2L || length(logit_se) != 2L) {
    stop("logit_auc and logit_se must each have length 2 (target, then European).")
  }
  auc <- stats::plogis(logit_auc)
  se_used <- auc_scale_se(auc, logit_se, se_scale)
  delta <- auc[[1]] - auc[[2]]
  se_delta <- sqrt(sum(se_used^2))
  c(lo = delta - 1.96 * se_delta, hi = delta + 1.96 * se_delta)
}

pool_stage2_cells <- function(paired_df, se_scale = c("logit_as_implemented", "delta_method")) {
  se_scale <- match.arg(se_scale)
  out <- paired_df |>
    dplyr::mutate(
      tgt_se_use = auc_scale_se(tgt_auc, tgt_se, se_scale),
      eur_se_use = auc_scale_se(eur_auc, eur_se, se_scale),
      delta = tgt_auc - eur_auc,
      delta_se = sqrt(tgt_se_use^2 + eur_se_use^2),
      flag_not_pooled = !is_pooled
    ) |>
    dplyr::filter(is.finite(.data$delta), is.finite(.data$delta_se), .data$delta_se > 0) |>
    dplyr::group_by(.data$trait_label, .data$trained_bucket, .data$target_ancestry) |>
    dplyr::summarise(
      n_pairs = dplyr::n(),
      delta_hat = sum(.data$delta / (.data$delta_se^2)) / sum(1 / (.data$delta_se^2)),
      se_hat = sqrt(1 / sum(1 / (.data$delta_se^2))),
      all_not_pooled = all(.data$flag_not_pooled, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      lo = .data$delta_hat - 1.96 * .data$se_hat,
      hi = .data$delta_hat + 1.96 * .data$se_hat
    )
  assert_one_row_per_stage2_cell(out)
}

# Pair a Stage-1 table on trait, score, training set, and covariate class.
# se_col is the standard-error column on scat_df. lo_pooled, hi_pooled, and
# se_logit are copied to eur_* / tgt_* when they are present. The class
# column is the class on both sides. Call pick_one_class_per_score() before
# the Stage-2 average so one score contributes one delta.
pair_stage1_by_class <- function(scat_df, se_col = "se") {
  required <- c(
    "trait_label", "ancestry_display", "pgs_id", "auc_pooled",
    "k_eval", "trained_bucket", "covariate_class", se_col
  )
  missing <- setdiff(required, names(scat_df))
  if (length(missing)) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  classes <- as.character(scat_df$covariate_class)
  unknown <- setdiff(unique(classes), covariate_class_levels)
  if (length(unknown) || anyNA(classes)) {
    stop(
      "Unknown covariate class: ",
      paste(unique(classes[!classes %in% covariate_class_levels]), collapse = ", "),
      call. = FALSE
    )
  }

  side_frame <- function(df, prefix) {
    out <- tibble::tibble(
      trait_label = as.character(df$trait_label),
      pgs_id = as.character(df$pgs_id),
      trained_bucket = as.character(df$trained_bucket),
      covariate_class = as.character(df$covariate_class)
    )
    out[[paste0(prefix, "_auc")]] <- df$auc_pooled
    out[[paste0(prefix, "_se")]] <- df[[se_col]]
    out[[paste0(prefix, "_k")]] <- df$k_eval
    if ("lo_pooled" %in% names(df)) out[[paste0(prefix, "_lo")]] <- df$lo_pooled
    if ("hi_pooled" %in% names(df)) out[[paste0(prefix, "_hi")]] <- df$hi_pooled
    if ("se_logit" %in% names(df)) out[[paste0(prefix, "_se_logit")]] <- df$se_logit
    out
  }

  eur_df <- side_frame(scat_df[scat_df$ancestry_display == "European", , drop = FALSE], "eur")
  tgt_src <- scat_df[scat_df$ancestry_display != "European", , drop = FALSE]
  tgt_df <- side_frame(tgt_src, "tgt")
  tgt_df$target_ancestry <- as.character(tgt_src$ancestry_display)

  tgt_df %>%
    dplyr::inner_join(
      eur_df,
      by = c("trait_label", "pgs_id", "trained_bucket", "covariate_class")
    ) %>%
    dplyr::mutate(is_pooled = (.data$eur_k >= 2) & (.data$tgt_k >= 2)) %>%
    dplyr::filter(
      is.finite(.data$eur_auc), is.finite(.data$tgt_auc),
      is.finite(.data$eur_se), is.finite(.data$tgt_se),
      .data$eur_se > 0, .data$tgt_se > 0
    )
}

safe_weighted_paired_test <- function(eur_auc, tgt_auc, eur_se, tgt_se,
                                      se_scale = c("logit_as_implemented", "delta_method")) {
  se_scale <- match.arg(se_scale)
  eur_se <- auc_scale_se(eur_auc, eur_se, se_scale)
  tgt_se <- auc_scale_se(tgt_auc, tgt_se, se_scale)

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
