# One evaluation per score, sample set, and ancestry.
# Pass class_col to keep one evaluation inside each covariate class, so a
# score-only model and a covariate-adjusted model in the same sample set
# both survive. Prefer an effective sample size, then the larger one, then
# the larger number of people, then the score-only model. Do not break ties
# with AUC. Inside one class the last tie breaks on performance_id.

.score_only_covariate_tokens <- c(
  "", "0", "none", "no", "na", "n/a", "unadjusted", "no covariates", "-"
)

is_score_only_covariates <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[is.na(x)] <- ""
  x %in% .score_only_covariate_tokens
}

calc_neff <- function(cases, controls) {
  cases <- as.numeric(cases)
  controls <- as.numeric(controls)
  out <- rep(NA_real_, length(cases))
  ok <- !is.na(cases) & !is.na(controls) & cases > 0 & controls > 0
  out[ok] <- 4 / (1 / cases[ok] + 1 / controls[ok])
  out
}

# Keep one row per (pgs_id, sampleset_id, ancestry), and per covariate
# class when class_col is set. Callers may pass pss_id / ppm_id / ctrls /
# ancestry_eval aliases through the arguments.
keep_one_evaluation <- function(df,
                                pgs_col = "pgs_id",
                                sampleset_col = "sampleset_id",
                                ancestry_col = "ancestry_eval",
                                performance_col = "performance_id",
                                n_col = "n",
                                cases_col = "cases",
                                controls_col = "controls",
                                covariates_col = "covariates",
                                class_col = NULL,
                                return_dropped = FALSE) {
  required <- c(
    pgs_col, sampleset_col, ancestry_col, performance_col,
    n_col, cases_col, controls_col, covariates_col
  )
  if (!is.null(class_col)) required <- c(required, class_col)
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  work <- df
  work$.pgs <- as.character(work[[pgs_col]])
  work$.pss <- as.character(work[[sampleset_col]])
  work$.anc <- as.character(work[[ancestry_col]])
  work$.ppm <- as.character(work[[performance_col]])
  work$.n <- as.numeric(work[[n_col]])
  work$.cases <- as.numeric(work[[cases_col]])
  work$.controls <- as.numeric(work[[controls_col]])
  work$.neff <- calc_neff(work$.cases, work$.controls)
  work$.score_only <- is_score_only_covariates(work[[covariates_col]])
  if (is.null(class_col)) {
    work$.class <- ""
  } else {
    work$.class <- as.character(work[[class_col]])
    if (anyNA(work$.class) || any(!nzchar(work$.class))) {
      stop("covariate_class is missing.", call. = FALSE)
    }
  }

  ranked <- work %>%
    dplyr::group_by(.data$.pgs, .data$.pss, .data$.anc, .data$.class) %>%
    dplyr::arrange(
      dplyr::desc(!is.na(.data$.neff)),
      dplyr::desc(.data$.neff),
      dplyr::desc(.data$.n),
      dplyr::desc(.data$.score_only),
      .data$.ppm,
      .by_group = TRUE
    ) %>%
    dplyr::mutate(
      one_eval_rank = dplyr::row_number(),
      kept_one_eval = .data$one_eval_rank == 1L
    ) %>%
    dplyr::ungroup()

  drop_helpers <- c(
    ".pgs", ".pss", ".anc", ".ppm", ".n", ".cases", ".controls",
    ".neff", ".score_only", ".class"
  )
  ranked <- ranked %>% dplyr::select(-dplyr::all_of(drop_helpers))

  kept <- ranked %>%
    dplyr::filter(.data$kept_one_eval) %>%
    dplyr::select(-dplyr::any_of(c("one_eval_rank", "kept_one_eval")))

  if (isTRUE(return_dropped)) {
    dropped <- ranked %>%
      dplyr::filter(!.data$kept_one_eval) %>%
      dplyr::select(-dplyr::any_of("kept_one_eval"))
    return(list(kept = kept, dropped = dropped, ranked = ranked))
  }
  kept
}
