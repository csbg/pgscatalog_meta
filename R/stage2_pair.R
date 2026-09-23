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
