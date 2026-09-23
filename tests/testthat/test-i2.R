test_that("I2 matches the hand calculation", {
  # weights 4 and 4, pooled eta 0.5, Q = 2, df = 1
  # I2 = (2 - 1) / 2 * 100 = 50
  out <- ivw_pool_logit(c(0, 1), c(0.5, 0.5))
  expect_equal(out$Q, 2)
  expect_equal(out$I2, 50)

  # Q is below df, so I2 is floored at 0
  near <- ivw_pool_logit(c(0, 0.01), c(1, 1))
  expect_equal(near$I2, 0)

  # one study has no heterogeneity
  one <- ivw_pool_logit(0.2, 0.1)
  expect_equal(one$k_eval, 1L)
  expect_true(is.na(one$I2))
})

test_that("I2 above 80 with two or more evaluations is the high-heterogeneity flag", {
  # weights 1 and 1, pooled eta 2, Q = 8, df = 1
  # I2 = (8 - 1) / 8 * 100 = 87.5
  high <- ivw_pool_logit(c(0, 4), c(1, 1))
  expect_equal(high$Q, 8)
  expect_equal(high$I2, 87.5)
  expect_equal(high$k_eval, 2L)
  expect_true(flag_high_i2(high$I2, high$k_eval))

  expect_false(flag_high_i2(80, 2))
  expect_false(flag_high_i2(90, 1))
  expect_false(flag_high_i2(NA_real_, 3))
})
