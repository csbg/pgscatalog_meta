# Confidence-interval parsers used by scripts 04 and 05.

coerce01 <- function(x) {
  xnum <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(xnum) & xnum > 1 & xnum <= 100, xnum / 100, xnum)
}

parse_est_ci <- function(x) {
  x <- as.character(x)
  est <- suppressWarnings(as.numeric(stringr::str_extract(x, "^[0-9]*\\.?[0-9]+")))
  lo  <- suppressWarnings(as.numeric(stringr::str_match(x, "\\[(\\d*\\.?\\d+),")[, 2]))
  hi  <- suppressWarnings(as.numeric(stringr::str_match(x, ",\\s*(\\d*\\.?\\d+)\\]")[, 2]))
  tibble::tibble(est = est, lo = lo, hi = hi)
}

parse_value_ci <- function(x) {
  parse_one <- function(x1) {
    s <- as.character(x1)
    if (is.na(s) || !nzchar(s)) return(c(NA_real_, NA_real_, NA_real_))
    s <- stringr::str_replace_all(s, "[\u2013\u2014\u2212]", "-")
    s <- stringr::str_squish(s)
    s <- stringr::str_replace_all(s, "(?i)\\b(AUROC|AUC|C-?index|C-stat(istic)?)\\s*[:=]?", "")
    s <- stringr::str_squish(s)

    num    <- "[-+]?\\d*[\\.,]?\\d+(?:[eE][-+]?\\d+)?"
    ci_any <- paste0("(?:\\[|\\()\\s*(", num, ")\\s*(?:,|;|\\-|to)\\s*(", num, ")\\s*(?:\\]|\\))")

    val1 <- stringr::str_match(s, paste0("\\b(", num, ")\\b"))[, 2]
    ci_m <- stringr::str_match(s, ci_any)
    lo1  <- ci_m[, 2]
    hi1  <- ci_m[, 3]
    if (is.na(lo1) || is.na(hi1)) {
      ci2 <- stringr::str_match(s, paste0("(?i)(?:ci|95%\\s*ci)\\s*(", num, ")\\s*(?:\\-|to)\\s*(", num, ")"))
      lo1 <- if (!is.na(ci2[, 1])) ci2[, 2] else lo1
      hi1 <- if (!is.na(ci2[, 1])) ci2[, 3] else hi1
    }

    norm_num <- function(z) suppressWarnings(as.numeric(gsub(",", ".", z, fixed = FALSE)))
    val_num <- norm_num(val1)
    lo_num <- norm_num(lo1)
    hi_num <- norm_num(hi1)

    clamp01 <- function(z) dplyr::case_when(
      is.na(z) ~ as.numeric(NA),
      z <= 1 ~ z,
      z > 1 & z <= 100 ~ z / 100,
      TRUE ~ z
    )
    c(
      estimate_value    = clamp01(val_num),
      estimate_ci_lower = clamp01(lo_num),
      estimate_ci_upper = clamp01(hi_num)
    )
  }
  m <- vapply(x, parse_one, FUN.VALUE = c(estimate_value = 0, estimate_ci_lower = 0, estimate_ci_upper = 0))
  tibble::as_tibble(t(m))
}
