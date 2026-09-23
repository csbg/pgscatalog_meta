# Reported-trait aliases. "Breast Carcinoma" is not mapped to "Breast cancer".

trait_dictionary <- c(
  "Type 2 Diabetes" = "Type 2 diabetes",
  "Breast Cancer" = "Breast cancer",
  "Prostate Cancer" = "Prostate cancer",
  "T1D" = "Type 1 diabetes",
  "Cancer of prostate" = "Prostate cancer",
  "Breast cancer [female]" = "Breast cancer",
  "Coronary artery disease" = "Coronary heart disease"
)

harmonize_trait_label <- function(efo_label, reported_trait) {
  na_empty <- function(x) {
    x <- as.character(x)
    ifelse(!is.na(x) & trimws(x) == "", NA_character_, x)
  }
  raw_label <- dplyr::coalesce(na_empty(efo_label), na_empty(reported_trait))
  temp_label <- trimws(raw_label)
  dplyr::recode(temp_label, !!!trait_dictionary)
}
