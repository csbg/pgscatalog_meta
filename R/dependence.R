# Covariance model for an ancestry difference when scores share people.
# The three calls are the metafor workflow for dependent estimates:
# vcalc builds the covariance, rma.mv fits the difference, robust gives
# the cluster-robust interval when the degrees of freedom allow it.
# Pustejovsky & Tipton (2022); metafor manual, "Some Recommended Practices".

ivw_estimate <- function(delta, se) {
  w <- 1 / se^2
  sum(w * delta) / sum(w)
}

difference_covariance <- function(se_eur, se_tgt, cluster_eur, cluster_tgt, pgs_id, rho) {
  if (any(!is.finite(c(se_eur, se_tgt, rho)))) {
    stop("Standard errors and rho must be finite.", call. = FALSE)
  }
  if (any(is.na(cluster_eur) | !nzchar(cluster_eur) | is.na(cluster_tgt) | !nzchar(cluster_tgt))) {
    stop("vcalc needs a non-missing cluster label on both sides.", call. = FALSE)
  }
  pgs_id <- as.character(pgs_id)
  if (anyDuplicated(pgs_id)) stop("Each score appears once in a comparison.", call. = FALSE)
  dat <- data.frame(
    vi_eur = se_eur^2,
    vi_tgt = se_tgt^2,
    cluster_eur = as.character(cluster_eur),
    cluster_tgt = as.character(cluster_tgt),
    pgs_id = pgs_id,
    stringsAsFactors = FALSE
  )
  V_eur <- metafor::vcalc(
    vi = vi_eur, cluster = cluster_eur, obs = pgs_id, data = dat, rho = rho
  )
  V_tgt <- metafor::vcalc(
    vi = vi_tgt, cluster = cluster_tgt, obs = pgs_id, data = dat, rho = rho
  )
  V_eur <- align_covariance(V_eur, pgs_id)
  V_tgt <- align_covariance(V_tgt, pgs_id)
  list(V = V_eur + V_tgt, V_eur = V_eur, V_tgt = V_tgt)
}

align_covariance <- function(V, pgs_id) {
  V <- as.matrix(V)
  if (!is.null(rownames(V)) && all(pgs_id %in% rownames(V))) {
    V <- V[pgs_id, pgs_id, drop = FALSE]
  }
  if (!isTRUE(all.equal(dim(V), c(length(pgs_id), length(pgs_id))))) {
    stop("vcalc returned a matrix of the wrong size.", call. = FALSE)
  }
  V
}

fit_difference_model <- function(delta, se_eur, se_tgt, cluster_eur, cluster_tgt, pgs_id,
                                 rho, sharing_group = NULL, sigma2 = NULL, min_robust_df = 4) {
  if (length(delta) == 1L) {
    se <- sqrt(se_eur^2 + se_tgt^2)
    return(list(
      estimate = delta, se = se, lo = delta - 1.96 * se, hi = delta + 1.96 * se,
      sig = (delta - 1.96 * se) > 0 | (delta + 1.96 * se) < 0,
      sigma2 = if (is.null(sigma2)) NA_real_ else sigma2,
      fit_note = "",
      V = matrix(se^2, 1, 1),
      V_eur = matrix(se_eur^2, 1, 1),
      V_tgt = matrix(se_tgt^2, 1, 1),
      robust_df = NA_real_, robust_reported = FALSE,
      robust_estimate = NA_real_, robust_se = NA_real_,
      robust_lo = NA_real_, robust_hi = NA_real_, robust_sig = NA
    ))
  }
  cov <- difference_covariance(se_eur, se_tgt, cluster_eur, cluster_tgt, pgs_id, rho)
  dat <- data.frame(yi = delta, pgs_id = as.character(pgs_id), stringsAsFactors = FALSE)
  warns <- character()
  fit <- withCallingHandlers(
    metafor::rma.mv(
      yi, cov$V,
      random = ~ 1 | pgs_id,
      data = dat,
      sigma2 = sigma2,
      test = "z"
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  est <- as.numeric(fit$beta)[1]
  se <- as.numeric(fit$se)[1]
  lo <- as.numeric(fit$ci.lb)[1]
  hi <- as.numeric(fit$ci.ub)[1]
  sig2 <- as.numeric(fit$sigma2)[1]
  out <- list(
    estimate = est,
    se = se,
    lo = lo,
    hi = hi,
    sig = is.finite(lo) && is.finite(hi) && (lo > 0 | hi < 0),
    sigma2 = sig2,
    fit_note = paste(unique(warns), collapse = " | "),
    V = cov$V,
    V_eur = cov$V_eur,
    V_tgt = cov$V_tgt,
    robust_df = NA_real_,
    robust_reported = FALSE,
    robust_estimate = NA_real_,
    robust_se = NA_real_,
    robust_lo = NA_real_,
    robust_hi = NA_real_,
    robust_sig = NA
  )
  n_groups <- if (is.null(sharing_group)) 0L else length(unique(sharing_group))
  if (is.null(sigma2) && n_groups >= 2L) {
    rob <- tryCatch(
      metafor::robust(fit, cluster = sharing_group, clubSandwich = TRUE),
      error = function(e) e
    )
    if (!inherits(rob, "error")) {
      ddf <- robust_df_value(rob)
      out$robust_df <- ddf
      if (is.finite(ddf) && ddf >= min_robust_df) {
        out$robust_reported <- TRUE
        out$robust_estimate <- as.numeric(rob$beta)[1]
        out$robust_se <- as.numeric(rob$se)[1]
        out$robust_lo <- as.numeric(rob$ci.lb)[1]
        out$robust_hi <- as.numeric(rob$ci.ub)[1]
        out$robust_sig <- out$robust_lo > 0 | out$robust_hi < 0
      }
    } else {
      out$fit_note <- paste(c(out$fit_note, conditionMessage(rob)), collapse = " | ")
    }
  }
  out
}

robust_df_value <- function(rob) {
  if (!is.null(rob$ddf)) return(as.numeric(rob$ddf)[1])
  if (!is.null(rob$dfs)) return(as.numeric(rob$dfs)[1])
  NA_real_
}

# Reported fits are only for comparisons whose scores share people.
# The rho = 0, sigma2 = 0 check is run for every comparison.
fit_shared_comparisons <- function(model_scores, comparisons,
                                   rho_main = 0.6,
                                   rho_grid = c(0, 0.3, 0.6, 0.9),
                                   min_robust_df = 4) {
  comp_key <- function(trait, bucket, ancestry) paste(trait, bucket, ancestry, sep = "\r")
  comp <- comparisons
  comp$goes <- comparisons$goes_through_covariance_model
  comp$key <- comp_key(comp$trait_label, comp$trained_bucket, comp$target_ancestry)
  model_scores$key <- comp_key(model_scores$trait_label, model_scores$trained_bucket, model_scores$target_ancestry)
  splits <- split(model_scores, model_scores$key)
  rows <- vector("list", nrow(comp))
  rho_rows <- list()
  check_gap <- numeric(nrow(comp))
  diag_gap <- numeric()
  for (i in seq_len(nrow(comp))) {
    df <- splits[[comp$key[[i]]]]
    if (is.null(df) || nrow(df) == 0L) stop("No scores for a comparison.", call. = FALSE)
    df <- df[order(df$pgs_id), , drop = FALSE]
    fe0 <- fit_difference_model(
      df$delta, df$eur_se, df$tgt_se, df$eur_cluster, df$tgt_cluster, df$pgs_id,
      rho = 0, sigma2 = 0
    )
    check_gap[[i]] <- abs(fe0$estimate - ivw_estimate(df$delta, df$delta_se))
    diag_gap <- c(diag_gap, max(abs(diag(fe0$V) - df$delta_se^2)))

    empty <- tibble::tibble(
      trait_label = comp$trait_label[[i]],
      trained_bucket = comp$trained_bucket[[i]],
      target_ancestry = comp$target_ancestry[[i]],
      delta_fe0 = fe0$estimate,
      delta_cov = NA_real_, se_cov = NA_real_, lo_cov = NA_real_, hi_cov = NA_real_,
      sig_cov = NA, sigma2_cov = NA_real_, sigma2_unstable = FALSE,
      fit_note = "",
      robust_df = NA_real_, robust_reported = FALSE,
      delta_robust = NA_real_, se_robust = NA_real_, lo_robust = NA_real_, hi_robust = NA_real_,
      sig_robust = NA,
      delta_cov_linked = NA_real_, se_cov_linked = NA_real_,
      lo_cov_linked = NA_real_, hi_cov_linked = NA_real_, sig_cov_linked = NA,
      delta_shared = NA_real_, se_shared = NA_real_,
      lo_shared = NA_real_, hi_shared = NA_real_, sig_shared = NA
    )
    if (!isTRUE(comp$goes[[i]])) {
      rows[[i]] <- empty
      next
    }
    grid <- lapply(rho_grid, function(rho) {
      fit_difference_model(
        df$delta, df$eur_se, df$tgt_se, df$eur_cluster, df$tgt_cluster, df$pgs_id,
        rho = rho, sharing_group = if (rho == rho_main) df$sharing_group else NULL,
        min_robust_df = min_robust_df
      )
    })
    names(grid) <- as.character(rho_grid)
    main <- grid[[as.character(rho_main)]]
    if (!is.finite(main$se)) {
      stop("Covariance model did not return a standard error for ", comp$trait_label[[i]],
           " / ", comp$target_ancestry[[i]], call. = FALSE)
    }
    linked <- fit_difference_model(
      df$delta, df$eur_se, df$tgt_se, df$eur_group, df$tgt_group, df$pgs_id,
      rho = rho_main
    )
    shared_only <- fit_difference_model(
      df$delta, df$eur_se, df$tgt_se, df$eur_cluster, df$tgt_cluster, df$pgs_id,
      rho = rho_main, sigma2 = 0
    )
    empty$delta_cov <- main$estimate
    empty$se_cov <- main$se
    empty$lo_cov <- main$lo
    empty$hi_cov <- main$hi
    empty$sig_cov <- main$sig
    empty$sigma2_cov <- main$sigma2
    empty$sigma2_unstable <- nrow(df) <= 3L
    empty$fit_note <- main$fit_note
    empty$robust_df <- main$robust_df
    empty$robust_reported <- main$robust_reported
    empty$delta_robust <- main$robust_estimate
    empty$se_robust <- main$robust_se
    empty$lo_robust <- main$robust_lo
    empty$hi_robust <- main$robust_hi
    empty$sig_robust <- main$robust_sig
    empty$delta_cov_linked <- linked$estimate
    empty$se_cov_linked <- linked$se
    empty$lo_cov_linked <- linked$lo
    empty$hi_cov_linked <- linked$hi
    empty$sig_cov_linked <- linked$sig
    empty$delta_shared <- shared_only$estimate
    empty$se_shared <- shared_only$se
    empty$lo_shared <- shared_only$lo
    empty$hi_shared <- shared_only$hi
    empty$sig_shared <- shared_only$sig
    rows[[i]] <- empty

    rho_rows[[length(rho_rows) + 1L]] <- tibble::tibble(
      trait_label = comp$trait_label[[i]],
      trained_bucket = comp$trained_bucket[[i]],
      target_ancestry = comp$target_ancestry[[i]],
      rho = rho_grid,
      delta_hat = vapply(grid, function(f) f$estimate, numeric(1)),
      se_hat = vapply(grid, function(f) f$se, numeric(1)),
      lo = vapply(grid, function(f) f$lo, numeric(1)),
      hi = vapply(grid, function(f) f$hi, numeric(1)),
      sig = vapply(grid, function(f) f$sig, logical(1)),
      sigma2 = vapply(grid, function(f) f$sigma2, numeric(1)),
      width = vapply(grid, function(f) f$hi - f$lo, numeric(1))
    )
    diag_gap <- c(diag_gap, max(abs(diag(main$V) - df$delta_se^2)))
  }
  list(
    by_comparison = dplyr::bind_rows(rows),
    rho = dplyr::bind_rows(rho_rows),
    max_abs_fe0 = max(check_gap),
    max_abs_diag = max(diag_gap)
  )
}

width_narrowed <- function(rho_tbl) {
  if (nrow(rho_tbl) == 0L) {
    return(rho_tbl[0, c("trait_label", "trained_bucket", "target_ancestry")])
  }
  rho_tbl |>
    dplyr::arrange(.data$trait_label, .data$trained_bucket, .data$target_ancestry, .data$rho) |>
    dplyr::group_by(.data$trait_label, .data$trained_bucket, .data$target_ancestry) |>
    dplyr::summarise(
      narrowed = any(diff(.data$width) < -1e-6),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$narrowed) |>
    dplyr::select("trait_label", "trained_bucket", "target_ancestry")
}
