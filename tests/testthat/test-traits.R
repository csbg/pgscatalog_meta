test_that("reported trait labels keep the existing dictionary", {
  expect_equal(harmonize_trait_label(NA_character_, "Breast Carcinoma"), "Breast Carcinoma")
  expect_equal(harmonize_trait_label(NA_character_, "Breast Cancer"), "Breast cancer")
  expect_equal(
    harmonize_trait_label("Breast Carcinoma", "Breast Cancer"),
    "Breast Carcinoma"
  )
})
