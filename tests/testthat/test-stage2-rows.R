test_that("Stage-2 pooling returns one row per trait, training bucket, and target ancestry", {
  paired <- tibble::tibble(
    trait_label = c("Trait A", "Trait A", "Trait B"),
    trained_bucket = c("Train: European-only", "Train: European-only", "Train: European-only"),
    target_ancestry = c("African", "African", "East Asian"),
    tgt_auc = c(0.60, 0.62, 0.55),
    eur_auc = c(0.70, 0.68, 0.66),
    tgt_se = c(0.02, 0.03, 0.02),
    eur_se = c(0.02, 0.02, 0.02),
    is_pooled = c(TRUE, FALSE, TRUE)
  )
  out <- pool_stage2_cells(paired, se_scale = "logit_as_implemented")
  expect_equal(nrow(out), 2L)
  expect_equal(out$n_pairs, c(2L, 1L))
  assert_one_row_per_stage2_cell(out)
})
