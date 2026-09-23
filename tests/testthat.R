library(testthat)

if (!file.exists("tests/testthat")) {
  stop("Run tests from the repository root.", call. = FALSE)
}

test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
