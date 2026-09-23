test_that("AUC strings parse to an estimate and interval", {
  parsed <- parse_value_ci(c("0.81 [0.81, 0.81]", "0.72 [0.70, 0.74]"))
  expect_equal(parsed$estimate_value, c(0.81, 0.72))
  expect_equal(parsed$estimate_ci_lower, c(0.81, 0.70))
  expect_equal(parsed$estimate_ci_upper, c(0.81, 0.74))
})
