test_that("study names and unnamed sample sets keep separate cluster labels", {
  expect_equal(normalize_study(c("UKB", "G&H", "AllofUs", "MEC|WHI")), c(
    "UK Biobank", "Genes & Health", "All of Us", "MEC|WHI"
  ))
  tokens <- sample_set_source_token(
    c("UKB", "UKB", NA, "NA", "FinnGen"),
    c("PSS1", "PSS2", "PSS3", "PSS4", "PSS5")
  )
  expect_equal(tokens[1:2], c("UK Biobank", "UK Biobank"))
  expect_equal(tokens[3:4], c("unnamed:PSS3", "unnamed:PSS4"))
  expect_false(any(is.na(tokens)))
  expect_equal(source_cluster_label(c("UK Biobank", "FinnGen", "UK Biobank")), "FinnGen + UK Biobank")
})

test_that("a shared study links scores only in the connected-group cluster", {
  pairs <- tibble::tibble(
    trait_label = "Example",
    trained_bucket = "Train: European-only",
    target_ancestry = "South Asian",
    pgs_id = c("PGS1", "PGS2", "PGS3"),
    covariate_class = "score_only",
    delta = c(-0.02, -0.04, -0.01),
    delta_se = c(0.03, 0.04, 0.05),
    eur_se = c(0.02, 0.03, 0.04),
    tgt_se = c(0.022, 0.026, 0.03)
  )
  keys <- tibble::tibble(
    pgs_id = c("PGS1", "PGS2", "PGS3", "PGS1", "PGS2", "PGS3"),
    trait_label = "Example",
    ancestry_display = c("European", "European", "European", "South Asian", "South Asian", "South Asian"),
    covariate_class = "score_only",
    sampleset_id = c("PSS_A", "PSS_B", "PSS_C", "PSS_T1", "PSS_T2", "PSS_T3")
  )
  sample_sets <- tibble::tibble(
    sampleset_id = c("PSS_A", "PSS_B", "PSS_C", "PSS_T1", "PSS_T2", "PSS_T3"),
    cohort_txt = c("UKB", "UKB", "FinnGen", "G&H", "PAGE", NA)
  )
  built <- build_source_sharing(pairs, keys, sample_sets)
  eur <- built$scores[built$scores$side == "European", ]
  expect_equal(eur$cluster_label[eur$pgs_id == "PGS1"], "UK Biobank")
  expect_equal(eur$cluster_label[eur$pgs_id == "PGS2"], "UK Biobank")
  expect_equal(eur$cluster_label[eur$pgs_id == "PGS3"], "FinnGen")
  tgt <- built$scores[built$scores$side == "Target", ]
  expect_equal(tgt$cluster_label[tgt$pgs_id == "PGS3"], "unnamed:PSS_T3")
  expect_equal(unnamed_token_collisions(built$scores), 0L)
  expect_true(built$comparisons$eur_shared_source)
  expect_false(built$comparisons$tgt_shared_source)
  expect_equal(built$comparisons$share_side, "European")
  expect_true(built$comparisons$source_adds_link)
  expect_equal(built$comparisons$n_sharing_groups, 2L)
  expect_equal(built$comparisons$n_partial_pairs, 0L)

  # Partial overlap: PGS2 also evaluated in FinnGen, so its label differs from PGS1
  # but the two still share UK Biobank.
  keys2 <- dplyr::bind_rows(keys, keys[2, ])
  keys2$sampleset_id[nrow(keys2)] <- "PSS_D"
  sample_sets2 <- dplyr::bind_rows(
    sample_sets,
    tibble::tibble(sampleset_id = "PSS_D", cohort_txt = "FinnGen")
  )
  sample_sets2$cohort_txt[sample_sets2$sampleset_id == "PSS_C"] <- "HUNT"
  built2 <- build_source_sharing(pairs, keys2, sample_sets2)
  eur2 <- built2$model_scores
  expect_equal(eur2$eur_cluster[eur2$pgs_id == "PGS1"], "UK Biobank")
  expect_equal(eur2$eur_cluster[eur2$pgs_id == "PGS2"], "FinnGen + UK Biobank")
  expect_equal(built2$comparisons$n_partial_pairs, 1L)
  expect_equal(eur2$eur_group[eur2$pgs_id == "PGS1"], eur2$eur_group[eur2$pgs_id == "PGS2"])
})

test_that("vcalc blocks add, and rho = 0 with sigma2 = 0 matches the inverse-variance average", {
  skip_if_not_installed("metafor")
  se_eur <- c(0.02, 0.03, 0.04)
  se_tgt <- c(0.02, 0.02, 0.05)
  delta <- c(-0.02, -0.05, 0.01)
  pgs <- c("A", "B", "C")
  # A and B share the European source. Targets are all different. C is independent.
  cov <- difference_covariance(
    se_eur, se_tgt,
    cluster_eur = c("UK Biobank", "UK Biobank", "FinnGen"),
    cluster_tgt = c("G1", "G2", "G3"),
    pgs_id = pgs,
    rho = 0.6
  )
  expect_equal(diag(cov$V_eur), se_eur^2)
  expect_equal(cov$V_eur[1, 2], 0.6 * 0.02 * 0.03)
  expect_equal(cov$V_eur[1, 3], 0)
  expect_equal(diag(cov$V), se_eur^2 + se_tgt^2)
  expect_equal(cov$V[1, 2], 0.6 * 0.02 * 0.03)

  delta_se <- sqrt(se_eur^2 + se_tgt^2)
  fe0 <- fit_difference_model(
    delta, se_eur, se_tgt,
    cluster_eur = c("UK Biobank", "UK Biobank", "FinnGen"),
    cluster_tgt = c("G1", "G2", "G3"),
    pgs_id = pgs,
    rho = 0,
    sigma2 = 0
  )
  expect_equal(fe0$estimate, ivw_estimate(delta, delta_se))
  expect_equal(sum(fe0$weight), 1)
  expect_equal(sum(fe0$weight * delta), fe0$estimate)
  expect_true(is.finite(fe0$reml_ll))

  # The linked European group puts rho on A and B even if the exact labels differed.
  linked <- difference_covariance(
    se_eur[1:2], se_tgt[1:2],
    cluster_eur = c("group", "group"),
    cluster_tgt = c("G1", "G2"),
    pgs_id = c("A", "B"),
    rho = 0.6
  )
  exact <- difference_covariance(
    se_eur[1:2], se_tgt[1:2],
    cluster_eur = c("UK Biobank", "FinnGen + UK Biobank"),
    cluster_tgt = c("G1", "G2"),
    pgs_id = c("A", "B"),
    rho = 0.6
  )
  expect_equal(exact$V[1, 2], 0)
  expect_equal(linked$V[1, 2], 0.6 * se_eur[1] * se_eur[2])
})
