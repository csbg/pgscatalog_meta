# Covariate class of a Catalog evaluation, and the one class kept per score.
# Three classes. A delta uses the same class on the European side and the
# non-European side. When more than one class is present on both sides,
# the score-only delta is kept, then age/sex/PCs/technical, then clinical.

covariate_class_levels <- c(
  "None (PRS only)",
  "Age/sex/PCs/technical",
  "Clinical or other covariates"
)

covariate_colour_levels <- c(
  "Uses covariates",
  "Does not use covariates",
  "Both kinds of score"
)

demo_tech_regex <- paste0(
  "^(age|sex|gender|pc\\d*|pcs?|principal|genetic pc|genotyp|array|batch|cent(er|re)|birth|",
  "deprivation|townsend|assessment|recruitment|ancestry|region|study|cohort|site|platform|chip|country|",
  "year|season|ses\\b|socio|education|income|\\d+\\s*pcs?|pc\\s*1|genetic principal|first \\d+ pc|top \\d+ pc|",
  "ehr|electronic health|duration of|follow-?up|enrol)"
)

classify_one_covariate <- function(txt) {
  if (is_score_only_covariates(txt)) return(covariate_class_levels[[1]])
  x <- tolower(trimws(as.character(txt)))
  toks <- stringr::str_split(x, "[,;+/]|\\band\\b|\\bwith\\b")[[1]]
  toks <- trimws(toks)
  toks <- toks[nzchar(toks)]
  if (!length(toks)) return(covariate_class_levels[[1]])
  demo <- grepl(demo_tech_regex, toks, perl = TRUE) |
    grepl("(pc|principal component|genotyp|array|batch|age|sex)", toks, perl = TRUE)
  if (all(demo)) covariate_class_levels[[2]] else covariate_class_levels[[3]]
}

classify_covariates <- function(txt) {
  vapply(txt, classify_one_covariate, character(1), USE.NAMES = FALSE)
}

# One row per group. Class priority is score-only, then age/sex/PCs/technical,
# then clinical. Two rows of the same class in one group are an error.
pick_one_covariate_class <- function(df, group_cols, class_col = "covariate_class") {
  if (!nrow(df)) return(df)
  missing <- setdiff(c(group_cols, class_col), names(df))
  if (length(missing)) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  classes <- as.character(df[[class_col]])
  if (anyNA(classes) || any(!nzchar(classes))) {
    stop("covariate_class is missing.", call. = FALSE)
  }
  unknown <- setdiff(unique(classes), covariate_class_levels)
  if (length(unknown)) {
    stop("Unknown covariate class: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  dup <- df %>%
    dplyr::count(dplyr::across(dplyr::all_of(c(group_cols, class_col)))) %>%
    dplyr::filter(.data$n > 1L)
  if (nrow(dup)) {
    stop("More than one row for the same score and covariate class.", call. = FALSE)
  }
  priority <- stats::setNames(seq_along(covariate_class_levels), covariate_class_levels)
  df %>%
    dplyr::mutate(.class_priority = unname(priority[as.character(.data[[class_col]])])) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::arrange(.data$.class_priority, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-".class_priority")
}

pick_one_class_per_score <- function(paired) {
  pick_one_covariate_class(
    paired,
    c("trait_label", "trained_bucket", "target_ancestry", "pgs_id")
  )
}

# Colour of one comparison, from the class of each score in its average.
comparison_covariate_colour <- function(classes) {
  classes <- unique(as.character(classes))
  classes <- classes[!is.na(classes) & nzchar(classes)]
  if (!length(classes)) return(NA_character_)
  has_none <- any(classes == covariate_class_levels[[1]])
  has_cov <- any(classes != covariate_class_levels[[1]])
  if (has_none && has_cov) return(covariate_colour_levels[[3]])
  if (has_none) return(covariate_colour_levels[[2]])
  if (has_cov) return(covariate_colour_levels[[1]])
  NA_character_
}
