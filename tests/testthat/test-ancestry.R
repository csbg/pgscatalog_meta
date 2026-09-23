test_that("evaluation ancestry keeps one allowed category and drops the rest", {
  middle_eastern <- "Greater Middle Eastern (Middle Eastern, North African or Persian)"
  expect_equal(
    assign_evaluation_ancestry(middle_eastern),
    "Middle Eastern or North African"
  )
  expect_equal(
    assign_evaluation_ancestry("African American or Afro-Caribbean"),
    "African"
  )
  expect_equal(
    assign_evaluation_ancestry("Sub-Saharan African"),
    "African"
  )
  expect_true(is.na(assign_evaluation_ancestry("African unspecified")))
  expect_true(is.na(assign_evaluation_ancestry(
    "East Asian, European, Hispanic or Latin American"
  )))
  south_east <- assign_evaluation_ancestry("South East Asian")
  expect_true(is.na(south_east))
  expect_false(identical(south_east, "East Asian"))
})

test_that("the display mapper used by the scripts is unchanged", {
  expect_equal(to_display_cat("African unspecified"), "African")
  expect_equal(to_display_cat("European"), "European")
})
