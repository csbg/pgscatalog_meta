test_that("blank and dash covariates are score-only; absolute-risk text is not", {
  expect_true(is_score_only_covariates("-"))
  expect_true(is_score_only_covariates(""))
  expect_true(is_score_only_covariates(NA_character_))
  expect_true(is_score_only_covariates("None"))
  expect_false(
    is_score_only_covariates("Age-specific absolute risk adjusted by PGS relative risk")
  )
})

breast_tie_df <- function() {
  tibble::tibble(
    pgs_id = "PGS005345",
    sampleset_id = c("PSS012140", "PSS012140", "PSS012137", "PSS012137"),
    ancestry_eval = "European",
    performance_id = c("PPM023137", "PPM023162", "PPM023147", "PPM023172"),
    n = c(212360, 212360, 2654, 2654),
    cases = c(3083, 3083, 194, 194),
    controls = c(209277, 209277, 2460, 2460),
    covariates = c(
      "-",
      "Age-specific absolute risk adjusted by PGS relative risk",
      "-",
      "Age-specific absolute risk adjusted by PGS relative risk"
    ),
    auc = c(0.641, 0.645, 0.658, 0.647)
  )
}

test_that("without a class column, breast-carcinoma ties keep the score-only row", {
  # Same Neff and n within a sample set; AUC would pick opposite models.
  kept <- keep_one_evaluation(breast_tie_df())
  expect_equal(nrow(kept), 2L)
  expect_equal(sort(kept$performance_id), c("PPM023137", "PPM023147"))
  expect_equal(sort(kept$sampleset_id), c("PSS012137", "PSS012140"))
})

test_that("breast-carcinoma ties keep both classes in each sample set", {
  df <- breast_tie_df() %>%
    dplyr::mutate(covariate_class = classify_covariates(covariates))
  kept <- keep_one_evaluation(df, class_col = "covariate_class")
  expect_equal(nrow(kept), 4L)
  expect_equal(
    sort(kept$performance_id),
    c("PPM023137", "PPM023147", "PPM023162", "PPM023172")
  )
  expect_equal(sort(unique(kept$sampleset_id)), c("PSS012137", "PSS012140"))
  expect_equal(length(unique(kept$covariate_class)), 2L)
  expect_true("None (PRS only)" %in% kept$covariate_class)
})

test_that("ties inside one covariate class keep the first performance id", {
  df <- tibble::tibble(
    pgs_id = "PGS000001",
    sampleset_id = "PSS000001",
    ancestry_eval = "European",
    performance_id = c("PPM000002", "PPM000001"),
    n = c(200, 200),
    cases = c(50, 50),
    controls = c(150, 150),
    covariates = c("age, sex", "age and sex"),
    covariate_class = c("Age/sex/PCs/technical", "Age/sex/PCs/technical")
  )
  kept <- keep_one_evaluation(df, class_col = "covariate_class")
  expect_equal(nrow(kept), 1L)
  expect_equal(kept$performance_id, "PPM000001")
})

test_that("different sample-set identifiers are both kept", {
  df <- tibble::tibble(
    pgs_id = "PGS005345",
    sampleset_id = c("PSS012140", "PSS012137"),
    ancestry_eval = "European",
    performance_id = c("PPM023137", "PPM023147"),
    n = c(212360, 2654),
    cases = c(3083, 194),
    controls = c(209277, 2460),
    covariates = c("-", "-")
  )

  kept <- keep_one_evaluation(df)
  expect_equal(nrow(kept), 2L)
  expect_equal(sort(kept$sampleset_id), c("PSS012137", "PSS012140"))
})

test_that("effective sample size outranks score-only when both differ", {
  df <- tibble::tibble(
    pgs_id = "PGS000001",
    sampleset_id = "PSS000001",
    ancestry_eval = "European",
    performance_id = c("PPM000002", "PPM000001"),
    n = c(100, 200),
    cases = c(10, 50),
    controls = c(90, 150),
    covariates = c("-", "age, sex")
  )

  kept <- keep_one_evaluation(df)
  expect_equal(kept$performance_id, "PPM000001")
})
