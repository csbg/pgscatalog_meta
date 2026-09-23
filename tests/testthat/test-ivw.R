test_that("Stage-1 IVW matches a two-study calculation and reports I2 as a percent", {
  out <- ivw_pool_logit(c(0, 1), c(0.5, 0.5))
  expect_equal(out$eta, 0.5)
  expect_equal(out$se, sqrt(1 / 8))
  expect_equal(out$I2, 50)
  expect_true(out$I2 >= 0 && out$I2 <= 100)
  expect_gt(out$I2, 1)

  near <- ivw_pool_logit(c(0, 0.01), c(1, 1))
  expect_equal(near$I2, 0)
  expect_true(near$I2 >= 0 && near$I2 <= 100)
})
