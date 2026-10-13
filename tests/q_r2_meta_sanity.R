#!/usr/bin/env Rscript

# ============================================================
# Sanity check for the R² sentences drafted for reviewer R.2.01.
#
# Those sentences say:
#   Yengo et al. 2022 has 14 height evaluations in the Catalog
#   (European R² 0.40; other ancestries about 0.09–0.21) and no
#   interval on any of them.
#   The only large set with a PGS-specific R², a confidence
#   interval, and both a European and a non-European evaluation
#   is 630 score–trait pairs, from Tanigawa et al. 2022
#   (PMID 35324888) and Tanigawa et al. 2023 (PMID 37890495).
#   The interval is on the unadjusted PGS R². The incremental R²
#   in those rows has no interval. Within a trait there is
#   usually only one score.
#
# "PGS-specific R²" is read here as the Catalog field
# "PGS R2 (no covariates)". A second table counts every other
# R² that already carries an interval, so that reading can be
# checked rather than assumed.
#
# Run from the repository root:
#   Rscript tests/q_r2_meta_sanity.R
#
# Tables land in
#   results/review_sensitivity_<stamp>/r2_meta_sanity/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

if (!file.exists("R/load.R")) {
  stop("R/load.R not found. Run this script from the repository root.", call. = FALSE)
}
source("R/load.R")

stamp <- pipeline_stamp
out_dir <- file.path("results", paste0("review_sensitivity_", stamp), "r2_meta_sanity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

num_pat <- "[-+]?(?:\\d+\\.?\\d*|\\d*\\.\\d+)(?:[eE][-+]?\\d+)?"

# label_regex identifies the metric. The optional bracket is the interval.
parse_metric_field <- function(text, label_regex) {
  pat <- sprintf(
    "%s\\s*=\\s*(%s)(?:\\s*\\[\\s*(%s)\\s*,\\s*(%s)\\s*\\])?",
    label_regex, num_pat, num_pat, num_pat
  )
  m <- str_match(text, regex(pat, ignore_case = TRUE))
  tibble(
    present = !is.na(m[, 2]),
    estimate = suppressWarnings(as.numeric(m[, 2])),
    ci_lower = suppressWarnings(as.numeric(m[, 3])),
    ci_upper = suppressWarnings(as.numeric(m[, 4]))
  )
}

read_catalog <- function(bulk_name) {
  read_csv(
    catalog_bulk_file(bulk_name, stamp),
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  )
}

message("Catalog stamp ", stamp)

ppm <- read_catalog("performance_metrics") %>%
  transmute(
    performance_id = `PGS Performance Metric (PPM) ID`,
    pgs_id = `Evaluated Score`,
    sampleset_id = `PGS Sample Set (PSS)`,
    pgp_id = `PGS Publication (PGP) ID`,
    reported_trait = `Reported Trait`,
    covariates = `Covariates Included in the Model`,
    other_info = `PGS Performance: Other Relevant Information`,
    pmid = `Publication (PMID)`,
    other_metric = `Other Metric(s)`,
    auc = `Area Under the Receiver-Operating Characteristic Curve (AUROC)`
  )

publications <- read_catalog("publications") %>%
  transmute(
    pgp_id = `PGS Publication/Study (PGP) ID`,
    first_author = `First Author`,
    title = Title,
    publication_year = substr(`Publication Date`, 1, 4),
    pmid_pub = `PubMed ID (PMID)`
  )

scores <- read_catalog("scores") %>%
  transmute(
    pgs_id = `Polygenic Score (PGS) ID`,
    pgs_name = `PGS Name`,
    n_variants = suppressWarnings(as.numeric(`Number of Variants`)),
    dev_method = `PGS Development Method`,
    gwas_ancestry = `Ancestry Distribution (%) - Source of Variant Associations (GWAS)`,
    score_pgp_id = `PGS Publication (PGP) ID`
  )

# One sample-set id can occupy more than one row. Keep the ancestry
# string only when every row agrees.
sample_sets <- read_catalog("evaluation_sample_sets") %>%
  transmute(
    sampleset_id = `PGS Sample Set (PSS)`,
    n_individuals = suppressWarnings(as.numeric(`Number of Individuals`)),
    n_cases = suppressWarnings(as.numeric(`Number of Cases`)),
    broad_ancestry_raw = `Broad Ancestry Category`,
    cohort = `Cohort(s)`,
    sample_note = `Additional Sample/Cohort Information`
  ) %>%
  group_by(sampleset_id) %>%
  summarise(
    n_pss_rows = n(),
    n_individuals = suppressWarnings(max(n_individuals, na.rm = TRUE)),
    any_cases = any(is.finite(n_cases)),
    broad_ancestry = {
      labs <- unique(broad_ancestry_raw[!is.na(broad_ancestry_raw) & nzchar(broad_ancestry_raw)])
      if (length(labs) == 1L) labs else paste(sort(labs), collapse = " | ")
    },
    ancestry_conflict = n_distinct(broad_ancestry_raw[!is.na(broad_ancestry_raw) & nzchar(broad_ancestry_raw)]) > 1,
    cohort = paste(unique(cohort[!is.na(cohort) & nzchar(cohort)]), collapse = " | "),
    sample_note = paste(unique(sample_note[!is.na(sample_note) & nzchar(sample_note)]), collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(
    n_individuals = ifelse(is.infinite(n_individuals), NA_real_, n_individuals),
    ancestry_display = assign_evaluation_ancestry(broad_ancestry),
    is_european = broad_ancestry == "European",
    is_non_european = !is.na(broad_ancestry) & nzchar(broad_ancestry) & broad_ancestry != "European"
  )

# The full-model field in this release is the superscript "R²".
# "PGS R2" and "Incremental R2" use an ASCII "R2", so the two
# patterns do not collide.
bare_r2 <- parse_metric_field(ppm$other_metric, "(?<![A-Za-z])R²")
pgs_r2 <- parse_metric_field(ppm$other_metric, "PGS R2 \\(no covariates\\)")
incr_r2 <- parse_metric_field(ppm$other_metric, "Incremental R2 \\(full-covars\\)")

metrics <- ppm %>%
  bind_cols(
    bare_r2 %>% rename_with(~ paste0("full_", .x)),
    pgs_r2 %>% rename_with(~ paste0("pgs_", .x)),
    incr_r2 %>% rename_with(~ paste0("incr_", .x))
  ) %>%
  mutate(
    pgs_has_ci = pgs_present & is.finite(pgs_ci_lower) & is.finite(pgs_ci_upper),
    full_has_ci = full_present & is.finite(full_ci_lower) & is.finite(full_ci_upper),
    incr_has_ci = incr_present & is.finite(incr_ci_lower) & is.finite(incr_ci_upper),
    any_r2_estimate = full_present | pgs_present | str_detect(
      coalesce(other_metric, ""),
      regex("R2|R²|variance explained", ignore_case = TRUE)
    )
  ) %>%
  left_join(sample_sets, by = "sampleset_id") %>%
  left_join(publications, by = "pgp_id") %>%
  left_join(scores, by = "pgs_id")

claims <- list()
add_claim <- function(id, statement, status, detail) {
  stopifnot(status %in% c("PASS", "PASS_WITH_SCOPE", "FAIL"))
  claims[[length(claims) + 1L]] <<- tibble(
    claim_id = id,
    status = status,
    statement = statement,
    detail = detail
  )
}

# ---------------------------------------------------------------------------
# Yengo et al. 2022 height
# ---------------------------------------------------------------------------

yengo <- metrics %>%
  filter(pgp_id == "PGP000382", reported_trait == "Height") %>%
  arrange(pgs_id, broad_ancestry, sampleset_id)

yengo_multi <- yengo %>% filter(pgs_id == "PGS002802")
yengo_multi_eur <- yengo_multi %>% filter(is_european)
yengo_multi_other <- yengo_multi %>% filter(is_non_european)
yengo_other_scores <- yengo %>% filter(pgs_id != "PGS002802", is_non_european)

yengo_n <- nrow(yengo)
yengo_ci_n <- sum(yengo$pgs_has_ci | yengo$full_has_ci | yengo$incr_has_ci)
yengo_multi_eur_r2 <- round(yengo_multi_eur$full_estimate, 2)
yengo_multi_other_r2 <- round(yengo_multi_other$full_estimate, 2)

yengo_ok <- yengo_n == 14L &&
  all(yengo$pmid == "36224396") &&
  all(yengo$first_author == "Yengo L") &&
  yengo_ci_n == 0L &&
  all(yengo$full_present) &&
  !any(yengo$ancestry_conflict) &&
  nrow(yengo_multi) == 7L &&
  all(yengo_multi_eur_r2 == 0.40) &&
  length(yengo_multi_other_r2) == 5L &&
  min(yengo_multi_other_r2) == 0.09 &&
  max(yengo_multi_other_r2) == 0.21

add_claim(
  "yengo_14_no_interval",
  "Yengo et al. 2022 has 14 height evaluations and no interval on any of them.",
  if (yengo_ok) "PASS" else "FAIL",
  sprintf(
    "PGP000382 / PMID 36224396: %d Height rows, %d with an interval. Multi-ancestry score PGS002802 European R² rounds to %s; its other ancestries round to %s–%s.",
    yengo_n,
    yengo_ci_n,
    paste(sprintf("%.2f", sort(unique(yengo_multi_eur_r2))), collapse = ", "),
    sprintf("%.2f", min(yengo_multi_other$full_estimate)),
    sprintf("%.2f", max(yengo_multi_other$full_estimate))
  )
)

outside_band <- yengo_other_scores %>%
  filter(round(full_estimate, 2) < 0.09 | round(full_estimate, 2) > 0.21)
lowest_yengo <- yengo %>% slice_min(full_estimate, n = 1, with_ties = FALSE)

add_claim(
  "yengo_range_scope",
  "European R² 0.40 and other ancestries about 0.09–0.21 describe the multi-ancestry score, not every one of the 14 rows.",
  if (nrow(outside_band) > 0L && yengo_ok) "PASS_WITH_SCOPE" else "FAIL",
  sprintf(
    "PGS002802 (%s, %s variants, %s) is the score evaluated in every ancestry group. These Catalog scores are not the 12,111-variant predictor in the reviewer's citation. The other %d rows are ancestry-specific scores. %d of those non-European rows fall outside 0.09–0.21 (lowest R² %.3f on %s in %s).",
    yengo_multi$pgs_name[1],
    format(yengo_multi$n_variants[1], big.mark = ","),
    yengo_multi$dev_method[1],
    sum(yengo$pgs_id != "PGS002802"),
    nrow(outside_band),
    lowest_yengo$full_estimate,
    lowest_yengo$pgs_id,
    lowest_yengo$broad_ancestry
  )
)

# ---------------------------------------------------------------------------
# Named field: PGS R2 (no covariates) with an interval
# ---------------------------------------------------------------------------

pgs_ci <- metrics %>% filter(pgs_has_ci)
pgs_no_ci <- metrics %>%
  filter(pgs_present, !pgs_has_ci) %>%
  mutate(
    nagelkerke_note = str_detect(
      coalesce(other_info, ""),
      regex("nagelkerke", ignore_case = TRUE)
    ),
    auroc_interval = str_detect(
      coalesce(other_metric, ""),
      regex("PGS AUROC \\(no covariates\\)\\s*=\\s*[-+0-9.eE]+\\s*\\[[^\\]]+\\]", ignore_case = TRUE)
    )
  )

pair_ancestry <- pgs_ci %>%
  group_by(pgs_id, reported_trait) %>%
  summarise(
    n_evaluations = n(),
    n_european = sum(is_european),
    n_non_european = sum(is_non_european),
    n_unlabelled = sum(!is_european & !is_non_european),
    ancestries = paste(sort(unique(broad_ancestry)), collapse = "; "),
    manuscript_ancestries = paste(sort(unique(na.omit(ancestry_display))), collapse = "; "),
    pmids = paste(sort(unique(pmid)), collapse = "; "),
    pgp_ids = paste(sort(unique(pgp_id)), collapse = "; "),
    first_authors = paste(sort(unique(first_author)), collapse = "; "),
    n_scores_check = n_distinct(pgs_id),
    incr_present_n = sum(incr_present),
    incr_ci_n = sum(incr_has_ci),
    full_ci_n = sum(full_has_ci),
    any_cases_n = sum(any_cases),
    auc_n = sum(!is.na(auc) & nzchar(auc)),
    .groups = "drop"
  )

both_ancestry <- pair_ancestry %>%
  filter(n_european > 0L, n_non_european > 0L)

trait_scores <- both_ancestry %>%
  group_by(reported_trait) %>%
  summarise(
    n_scores = n_distinct(pgs_id),
    pmids = paste(sort(unique(unlist(str_split(pmids, "; ")))), collapse = "; "),
    .groups = "drop"
  )

within_paper <- pgs_ci %>%
  distinct(pmid, reported_trait, pgs_id) %>%
  count(pmid, reported_trait, name = "n_scores")

n_pairs <- nrow(both_ancestry)
pair_pmids <- sort(unique(unlist(str_split(both_ancestry$pmids, "; "))))
expected_pmids <- c("35324888", "37890495")

tanigawa_structure_ok <- n_pairs == 630L &&
  identical(pair_pmids, expected_pmids) &&
  all(both_ancestry$n_evaluations == 5L) &&
  all(both_ancestry$n_european == 2L) &&
  all(both_ancestry$incr_present_n == both_ancestry$n_evaluations) &&
  all(both_ancestry$incr_ci_n == 0L) &&
  all(both_ancestry$full_ci_n == both_ancestry$n_evaluations) &&
  all(both_ancestry$any_cases_n == 0L) &&
  all(both_ancestry$auc_n == 0L) &&
  nrow(pgs_ci) == 3150L &&
  setequal(pgs_ci$pmid, expected_pmids)

add_claim(
  "tanigawa_630_pairs",
  "630 score–trait pairs have a PGS R2 (no covariates) interval and both a European and a non-European evaluation, and they come from Tanigawa 2022 and Tanigawa 2023.",
  if (tanigawa_structure_ok) "PASS" else "FAIL",
  sprintf(
    "%d evaluations, %d score–trait pairs, PMIDs %s. Each pair has 5 evaluations (2 European). Incremental R² is present on every row and has an interval on %d rows. The full-model R² on the same rows does have an interval (%d/%d).",
    nrow(pgs_ci),
    n_pairs,
    paste(pair_pmids, collapse = ", "),
    sum(pgs_ci$incr_has_ci),
    sum(pgs_ci$full_has_ci),
    nrow(pgs_ci)
  )
)

one_score_within_paper <- all(within_paper$n_scores == 1L)
n_traits <- nrow(trait_scores)
n_traits_one_score <- sum(trait_scores$n_scores == 1L)
n_traits_two_scores <- sum(trait_scores$n_scores == 2L)

add_claim(
  "tanigawa_r2_without_interval",
  "The remaining named PGS R² values have no interval, so they cannot be meta-analyzed either.",
  if (nrow(pgs_no_ci) == 1189L && all(pgs_no_ci$pmid == "35324888") && all(pgs_no_ci$any_cases) && all(pgs_no_ci$nagelkerke_note) && all(pgs_no_ci$auroc_interval)) {
    "PASS"
  } else {
    "FAIL"
  },
  sprintf(
    "%d rows, all Tanigawa 2022 (PMID 35324888), all with case counts, all noting Nagelkerke's R². The interval in those rows is on the AUROC, not on the R².",
    nrow(pgs_no_ci)
  )
)

add_claim(
  "one_score_per_trait",
  "Within a trait there is usually only one score.",
  if (one_score_within_paper && n_traits_one_score == 586L && n_traits_two_scores == 22L) "PASS" else "FAIL",
  sprintf(
    "Within each paper every trait has one score. Across the two papers, %d of %d traits have one score and %d have two (the same trait released in both papers).",
    n_traits_one_score, n_traits, n_traits_two_scores
  )
)

# The two European rows are UK Biobank partitions, not two studies.
eur_notes <- pgs_ci %>%
  filter(is_european) %>%
  group_by(pgs_id, reported_trait) %>%
  summarise(
    n_holdout = sum(str_detect(sample_note, regex("holdout|testing cohort", ignore_case = TRUE))),
    n_blank_note = sum(!nzchar(sample_note)),
    cohorts = paste(sort(unique(cohort)), collapse = "; "),
    .groups = "drop"
  )

add_claim(
  "not_a_multi_score_meta_analysis",
  "These 630 pairs do not support a within-trait, multi-score meta-analysis.",
  if (one_score_within_paper && all(eur_notes$cohorts == "UKB") && all(eur_notes$n_holdout + eur_notes$n_blank_note == 2L)) {
    "PASS"
  } else {
    "FAIL"
  },
  sprintf(
    "Each paper contributes one score per trait. All %d European pairs are UK Biobank. %d of 630 pairs label one European sample as a held-out testing cohort; the other European sample is a second UK Biobank partition. Non-European samples are one UK Biobank set per ancestry label.",
    nrow(eur_notes),
    sum(eur_notes$n_holdout == 1L)
  )
)

manuscript_non_eur <- pgs_ci %>%
  filter(is_non_european) %>%
  count(broad_ancestry, ancestry_display, name = "n")

# ---------------------------------------------------------------------------
# Other R² intervals, so "the only large set" can be checked
# ---------------------------------------------------------------------------

# An interval on some R² or "variance explained" phrase, excluding the
# named PGS R2 (no covariates) field already counted above.
other_ci_pat <- regex(
  "(?i)((?:incremental\\s+variance\\s+explained[^,=]{0,80}|R²|R2))\\s*=\\s*[-+0-9.eE]+\\s*\\[[^\\]]+\\]",
  ignore_case = TRUE
)

other_ci <- metrics %>%
  filter(!pgs_has_ci) %>%
  mutate(has_other_r2_ci = str_detect(coalesce(other_metric, ""), other_ci_pat)) %>%
  filter(has_other_r2_ci)

other_pairs <- other_ci %>%
  group_by(pmid, first_author, publication_year, pgp_id, title, pgs_id, reported_trait) %>%
  summarise(
    n_evaluations = n(),
    n_european = sum(is_european),
    n_non_european = sum(is_non_european),
    ancestries = paste(sort(unique(broad_ancestry)), collapse = "; "),
    covariates = paste(sort(unique(covariates)), collapse = " || "),
    other_info = paste(sort(unique(na.omit(other_info[nzchar(other_info)]))), collapse = " || "),
    .groups = "drop"
  )

other_both <- other_pairs %>%
  filter(n_european > 0L, n_non_european > 0L)

other_by_paper <- other_both %>%
  group_by(pmid, first_author, publication_year, pgp_id, title) %>%
  summarise(
    n_score_trait_pairs = n(),
    n_traits = n_distinct(reported_trait),
    n_scores = n_distinct(pgs_id),
    max_scores_in_one_trait = max(as.numeric(table(reported_trait))),
    ancestries = paste(sort(unique(unlist(str_split(ancestries, "; ")))), collapse = "; "),
    covariate_values = paste(sort(unique(covariates)), collapse = " || "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_score_trait_pairs))

large_other <- other_by_paper %>% filter(n_score_trait_pairs >= 20L)

add_claim(
  "only_large_set_named_field",
  "Under the named field 'PGS R2 (no covariates)', no other publication has a score–trait pair with an interval and both ancestries.",
  if (nrow(pgs_ci) == 3150L && setequal(unique(pgs_ci$pmid), expected_pmids)) "PASS" else "FAIL",
  sprintf(
    "Publications with that named field and an interval: %s.",
    paste(sort(unique(paste(pgs_ci$first_author, pgs_ci$pmid))), collapse = "; ")
  )
)

monti <- other_both %>% filter(pmid == "38908374")
monti_quant <- monti %>% filter(str_detect(other_info, fixed("pearson correlation")))

large_other_detail <- if (nrow(large_other) == 0L) {
  "No other paper has 20 or more score–trait pairs with an R² interval and both ancestries."
} else {
  paste(
    "That sentence is true for the named field only.",
    paste(sprintf(
      "%s %s (PMID %s) has %d score–trait pairs (%d traits, up to %d scores in one trait).",
      large_other$first_author,
      large_other$publication_year,
      large_other$pmid,
      large_other$n_score_trait_pairs,
      large_other$n_traits,
      large_other$max_scores_in_one_trait
    ), collapse = " ")
  )
}

add_claim(
  "only_large_set_any_r2",
  "The two Tanigawa papers are the only large cross-ancestry set with an R² interval.",
  if (nrow(large_other) == 0L) "PASS" else "PASS_WITH_SCOPE",
  large_other_detail
)

# ---------------------------------------------------------------------------
# Write the audit tables and print the verdict
# ---------------------------------------------------------------------------

yengo_out <- yengo %>%
  transmute(
    performance_id, pgs_id, pgs_name, n_variants, dev_method, gwas_ancestry,
    reported_trait, sampleset_id, broad_ancestry, ancestry_display,
    n_individuals, cohort, covariates,
    r2 = full_estimate, r2_ci_lower = full_ci_lower, r2_ci_upper = full_ci_upper,
    other_metric, pmid, pgp_id
  )

pairs_out <- both_ancestry %>%
  left_join(
    pgs_ci %>% distinct(pgs_id, pgs_name, dev_method),
    by = "pgs_id"
  ) %>%
  arrange(pmids, reported_trait, pgs_id)

write_csv(yengo_out, file.path(out_dir, "yengo_height_evaluations.csv"))
write_csv(pairs_out, file.path(out_dir, "pgs_r2_ci_score_trait_pairs.csv"))
write_csv(manuscript_non_eur, file.path(out_dir, "pgs_r2_ci_non_european_labels.csv"))
write_csv(other_by_paper, file.path(out_dir, "other_r2_ci_cross_ancestry_papers.csv"))
write_csv(trait_scores, file.path(out_dir, "pgs_r2_ci_scores_per_trait.csv"))

claim_tbl <- bind_rows(claims)
write_csv(claim_tbl, file.path(out_dir, "claim_checks.csv"))

cat("\n")
cat("R² meta-analysis sanity check, Catalog stamp ", stamp, "\n", sep = "")
cat(strrep("=", 72), "\n", sep = "")
for (i in seq_len(nrow(claim_tbl))) {
  cat("\n[", claim_tbl$status[i], "] ", claim_tbl$claim_id[i], "\n", sep = "")
  cat(claim_tbl$statement[i], "\n", sep = "")
  cat(claim_tbl$detail[i], "\n", sep = "")
}

cat("\n", strrep("-", 72), "\n", sep = "")
cat("Yengo height rows (multi-ancestry score first)\n")
yengo_print <- yengo_out %>%
  mutate(role = if_else(pgs_id == "PGS002802", "multi-ancestry score", "ancestry-specific score")) %>%
  arrange(desc(role == "multi-ancestry score"), broad_ancestry, pgs_id)
yengo_show <- yengo_print %>%
  transmute(
    role, pgs_id, broad_ancestry,
    n = format(n_individuals, big.mark = ","),
    r2 = sprintf("%.3f", r2),
    variants = format(n_variants, big.mark = ",")
  )
print(as.data.frame(yengo_show), row.names = FALSE)

cat("\nNon-European labels on the 3,150 PGS-R²-with-interval rows\n")
cat("ancestry_display is the manuscript rule in R/ancestry.R; NA means that rule drops the label.\n")
print(as.data.frame(manuscript_non_eur), row.names = FALSE)

cat("\nOther papers with an R² interval and both ancestries\n")
if (nrow(other_by_paper) == 0L) {
  cat("None.\n")
} else {
  print(as.data.frame(other_by_paper %>%
    transmute(
      year = publication_year,
      author = first_author,
      pmid,
      pairs = n_score_trait_pairs,
      traits = n_traits,
      max_scores_one_trait = max_scores_in_one_trait
    )), row.names = FALSE)
}

cat("\nMonti et al. 2024 quantitative rows (performance note contains 'pearson correlation'): ",
    nrow(monti_quant), " score–trait pairs with both ancestries.\n", sep = "")
cat("Monti covariates field, unique values: ",
    paste(sort(unique(monti$covariates)), collapse = " | "), "\n", sep = "")

cat("\nAnswers to the two supervisor questions\n")
cat("1. Can we meta-analyze Catalog R²?\n")
cat("   Yengo height: no. None of the 14 rows has an interval, so there is no standard error.\n")
cat(sprintf(
  "   Tanigawa also reports PGS R² with no interval on %d rows (PMID %s). Those rows have case counts, and the performance note says the R² is Nagelkerke's.\n",
  nrow(pgs_no_ci),
  paste(sort(unique(pgs_no_ci$pmid)), collapse = ", ")
))
cat("   Where the unadjusted PGS R² does have an interval, each paper still has one score per trait,\n")
cat("   and each non-European ancestry is one UK Biobank sample. That is a single estimate per cell,\n")
cat("   not a meta-analysis across scores. The incremental R², which removes the covariates, has no interval.\n")
cat("2. How many score–trait pairs have a PGS-specific R², an interval, and both ancestries?\n")
cat("   Named field 'PGS R2 (no covariates)': ", n_pairs, " pairs, both Tanigawa UK Biobank papers.\n", sep = "")
cat("   Any other R² that already has an interval: see the paper table above. The large extra set is Monti et al. 2024.\n")

cat("\nWrote tables to ", out_dir, "\n", sep = "")

n_fail <- sum(claim_tbl$status == "FAIL")
if (n_fail > 0L) {
  stop(sprintf("%d claim(s) failed. See %s", n_fail, file.path(out_dir, "claim_checks.csv")), call. = FALSE)
}
