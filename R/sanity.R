# Checks on pooled summaries.

assert_i2_percent <- function(i2) {
  bad <- !is.na(i2) & (i2 < 0 | i2 > 100)
  if (any(bad)) stop("I2 is outside 0-100.")
  invisible(i2)
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
