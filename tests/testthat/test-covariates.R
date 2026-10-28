test_that("covariate class separates score-only, demographic, and clinical text", {
  expect_equal(classify_covariates(c("", "-", "None", NA_character_)), rep("None (PRS only)", 4))
  expect_equal(classify_covariates("age, sex, PC1"), "Age/sex/PCs/technical")
  expect_equal(classify_covariates("age, BMI"), "Clinical or other covariates")
  expect_false(is_score_only_covariates("Age-specific absolute risk adjusted by PGS relative risk"))
})

test_that("pairing keeps one class on both sides and one class per score", {
  scat <- tibble::tibble(
    trait_label = "Trait",
    trained_bucket = "Train: European-only",
    pgs_id = "PGS1",
    ancestry_display = c("European", "European", "African", "African"),
    covariate_class = c(
      "None (PRS only)", "Clinical or other covariates",
      "None (PRS only)", "Clinical or other covariates"
    ),
    auc_pooled = c(0.70, 0.80, 0.60, 0.81),
    se = c(0.02, 0.02, 0.02, 0.02),
    k_eval = 1
  )
  paired <- pair_stage1_by_class(scat, se_col = "se")
  expect_equal(nrow(paired), 2L)
  expect_equal(sort(paired$covariate_class), sort(covariate_class_levels[c(1, 3)]))
  picked <- pick_one_class_per_score(paired)
  expect_equal(picked$covariate_class, "None (PRS only)")
  expect_equal(picked$eur_auc, 0.70)
  expect_equal(picked$tgt_auc, 0.60)
  expect_equal(
    comparison_covariate_colour(paired$covariate_class),
    "Both kinds of score"
  )
  expect_equal(
    comparison_covariate_colour(picked$covariate_class),
    "Does not use covariates"
  )
  expect_equal(
    comparison_covariate_colour(c("Age/sex/PCs/technical", "Clinical or other covariates")),
    "Uses covariates"
  )
})

test_that("a class present on one side only is not a delta", {
  scat <- tibble::tibble(
    trait_label = "Trait",
    trained_bucket = "Train: European-only",
    pgs_id = "PGS1",
    ancestry_display = c("European", "European", "African"),
    covariate_class = c(
      "None (PRS only)", "Clinical or other covariates", "Clinical or other covariates"
    ),
    auc_pooled = c(0.70, 0.80, 0.81),
    se = c(0.02, 0.02, 0.02),
    k_eval = 1
  )
  paired <- pair_stage1_by_class(scat, se_col = "se")
  expect_equal(nrow(paired), 1L)
  expect_equal(paired$covariate_class, "Clinical or other covariates")
  expect_equal(pick_one_class_per_score(paired)$eur_auc, 0.80)
})
