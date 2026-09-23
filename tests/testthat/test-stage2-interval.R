test_that("PGS000004 intervals match the submitted and delta-method bounds", {
  logit_auc <- c(0.291508, 0.405465)
  logit_se <- c(0.01630624, 0.02126195)

  submitted <- stage2_delta_interval(logit_auc, logit_se, "logit_as_implemented")
  expect_lt(max(abs(submitted - c(-0.080153, 0.024883))), 1e-5)

  delta_method <- stage2_delta_interval(logit_auc, logit_se, "delta_method")
  expect_lt(max(abs(delta_method - c(-0.040332, -0.014937))), 1e-5)

  auc <- stats::plogis(logit_auc)
  one <- tibble::tibble(
    trait_label = "example",
    trained_bucket = "Train: European-only",
    target_ancestry = "African",
    tgt_auc = auc[[1]],
    eur_auc = auc[[2]],
    tgt_se = logit_se[[1]],
    eur_se = logit_se[[2]],
    is_pooled = TRUE
  )
  pooled <- pool_stage2_cells(one, se_scale = "logit_as_implemented")
  expect_equal(nrow(pooled), 1L)
  expect_lt(max(abs(c(pooled$lo, pooled$hi) - c(-0.080153, 0.024883))), 1e-5)

  pooled_delta <- pool_stage2_cells(one, se_scale = "delta_method")
  expect_lt(max(abs(c(pooled_delta$lo, pooled_delta$hi) - c(-0.040332, -0.014937))), 1e-5)
})

test_that("the weighted paired test uses the selected standard-error scale", {
  eur_auc <- c(0.80, 0.70, 0.60)
  tgt_auc <- c(0.60, 0.65, 0.55)
  eur_se <- c(0.02, 0.05, 0.01)
  tgt_se <- c(0.02, 0.05, 0.01)
  submitted <- safe_weighted_paired_test(
    eur_auc, tgt_auc, eur_se, tgt_se, se_scale = "logit_as_implemented"
  )
  delta_method <- safe_weighted_paired_test(
    eur_auc, tgt_auc, eur_se, tgt_se, se_scale = "delta_method"
  )
  expect_equal(submitted$n_pairs, 3L)
  expect_gt(abs(submitted$estimate - delta_method$estimate), 1e-8)
})
