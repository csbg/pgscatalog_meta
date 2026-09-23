# Checks on pooled summaries and on AUC intervals that would enter Stage 1.

auc_interval_issue <- function(lower, upper) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  finite <- is.finite(lower) & is.finite(upper)
  outside <- finite & (lower < 0 | lower > 1 | upper < 0 | upper > 1)
  nonpos <- finite & !outside & !((upper - lower) > 0)
  issue <- rep(NA_character_, length(lower))
  issue[outside] <- "bound_outside_unit_interval"
  issue[nonpos] <- "non_positive_width"
  issue
}

assert_i2_percent <- function(i2) {
  bad <- !is.na(i2) & (i2 < 0 | i2 > 100)
  if (any(bad)) stop("I2 is outside 0-100.")
  invisible(i2)
}

flag_high_i2 <- function(i2, k_eval, threshold = 80) {
  !is.na(i2) & k_eval >= 2 & i2 > threshold
}

assert_one_row_per_stage2_cell <- function(df,
                                           keys = c("trait_label", "trained_bucket", "target_ancestry")) {
  missing <- setdiff(keys, names(df))
  if (length(missing)) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  n_keys <- nrow(dplyr::distinct(df, dplyr::across(dplyr::all_of(keys))))
  if (n_keys != nrow(df)) {
    stop(
      "Stage-2 table has ", nrow(df), " rows for ", n_keys, " cells.",
      call. = FALSE
    )
  }
  invisible(df)
}
