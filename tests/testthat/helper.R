suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
})

root <- normalizePath(testthat::test_path("..", ".."))
source(file.path(root, "R", "load.R"))
