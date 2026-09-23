# Ancestry helpers.
# to_display_cat() and make_ancestry_display() are the mappers the scripts
# already used. assign_evaluation_ancestry() is the evaluation-sample rule.

to_display_cat <- function(x) {
  x <- gsub("\\s*,\\s*", ",", x)
  sapply(strsplit(x, ","), function(v) {
    v <- unique(trimws(v))
    if (length(v) > 1) {
      if ("European" %in% v) "Multi-ancestry including European" else "Multi-ancestry excluding European"
    } else {
      vv <- v[1]
      dplyr::case_when(
        vv %in% c("European", "African", "East Asian", "South Asian",
                  "Hispanic or Latin American", "Middle Eastern or North African",
                  "Other/Mixed", "Not reported",
                  "Multi-ancestry including European", "Multi-ancestry excluding European") ~ vv,
        vv %in% c("African American or Afro-Caribbean", "African unspecified", "Sub-Saharan African") ~ "African",
        vv == "East Asian" ~ "East Asian",
        vv == "South Asian" ~ "South Asian",
        vv == "European" ~ "European",
        vv == "Hispanic or Latin American" ~ "Hispanic or Latin American",
        vv == "Greater Middle Eastern (Middle Eastern, North African or Persian)" ~ "Middle Eastern or North African",
        vv %in% c("Central Asian", "South East Asian", "Asian unspecified", "Oceanian", "Native American",
                  "Aboriginal Australian", "Other", "Other admixed ancestry") ~ "Other/Mixed",
        vv == "Not reported" ~ "Not reported",
        TRUE ~ "Other/Mixed"
      )
    }
  }, USE.NAMES = FALSE)
}

make_ancestry_display <- function(x) {
  x <- as.character(x)
  g <- function(p) grepl(p, x, ignore.case = TRUE, perl = TRUE)
  dplyr::case_when(
    is.na(x) | x == "" | g("Unknown|Not reported|\\bNR\\b|unspecified") ~ "Not reported",
    g("Multi-ancestry including European") ~ "Multi-ancestry including European",
    g("Multi-ancestry excluding European") ~ "Multi-ancestry excluding European",
    g("^European(,|$)|^EUR(,|$)|^European\\b") ~ "European",
    g("Sub-?Saharan|\\bSSA\\b") ~ "African",
    g("African\\s*American|Afro-?Caribbean|Black or African American|Black/African American") ~ "African",
    (g("^African") & !g("^African\\s*American")) ~ "African",
    g("Middle|Greater Middle|North African|\\bMENA\\b") ~ "Middle Eastern or North African",
    g("East\\s*Asian")  ~ "East Asian",
    g("South\\s*Asian") ~ "South Asian",
    g("Hispanic|Lat(in|inx|ino)|Latin American") ~ "Hispanic or Latin American",
    g("Other|Mixed|Admixed") ~ "Other/Mixed",
    TRUE ~ x
  )
}

# Exact Catalog broad_ancestry_category values that enter an evaluation.
# Commas inside parentheses are part of one label.
.evaluation_ancestry_map <- c(
  "European" = "European",
  "African American or Afro-Caribbean" = "African",
  "Sub-Saharan African" = "African",
  "East Asian" = "East Asian",
  "South Asian" = "South Asian",
  "Hispanic or Latin American" = "Hispanic or Latin American",
  "Greater Middle Eastern (Middle Eastern, North African or Persian)" = "Middle Eastern or North African"
)

split_ancestry_categories <- function(x) {
  if (length(x) != 1L) stop("split_ancestry_categories() takes one string.")
  if (is.na(x) || !nzchar(trimws(x))) return(character())
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  depth <- 0L
  buf <- character()
  parts <- character()
  for (ch in chars) {
    if (ch == "(") {
      depth <- depth + 1L
      buf <- c(buf, ch)
    } else if (ch == ")" && depth > 0L) {
      depth <- depth - 1L
      buf <- c(buf, ch)
    } else if (ch == "," && depth == 0L) {
      parts <- c(parts, paste(buf, collapse = ""))
      buf <- character()
    } else {
      buf <- c(buf, ch)
    }
  }
  parts <- c(parts, paste(buf, collapse = ""))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

assign_evaluation_ancestry <- function(x) {
  x <- as.character(x)
  vapply(x, function(one) {
    parts <- split_ancestry_categories(one)
    if (length(parts) != 1L || !(parts %in% names(.evaluation_ancestry_map))) {
      return(NA_character_)
    }
    unname(.evaluation_ancestry_map[[parts]])
  }, character(1), USE.NAMES = FALSE)
}
