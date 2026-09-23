test_that("bad AUC intervals are the zero-width and out-of-range cases", {
  parsed <- parse_value_ci(c("0.81 [0.81,0.81]", "0.721 [0.704,737.0]", "0.72 [0.70, 0.74]"))
  issue <- auc_interval_issue(parsed$estimate_ci_lower, parsed$estimate_ci_upper)
  expect_equal(issue[[1]], "non_positive_width")
  expect_equal(issue[[2]], "bound_outside_unit_interval")
  expect_true(is.na(issue[[3]]))
})

test_that("the named catalog stamp points at the frozen download", {
  expect_equal(pipeline_stamp, "20260923")
  expect_true(file.exists(catalog_bulk_file("performance_metrics")))
  expect_true(file.exists(catalog_bulk_file("evaluation_sample_sets")))
})
