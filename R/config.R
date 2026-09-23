# Named input for this revision. Scripts read this stamp instead of the newest file.

pipeline_stamp <- "20260923"

pipeline_repo_root <- function() {
  if (exists("pipeline_r_dir", inherits = TRUE)) dirname(pipeline_r_dir) else "."
}

# Paper Catalog outputs: Stage-1 tables and Stage-2 / Figure 1C files.
catalog_results_dir <- function(stamp = pipeline_stamp) {
  file.path(pipeline_repo_root(), "results", paste0("pgs_catalog_", stamp))
}

# Per-trait exploratory PDFs and run logs (not manuscript Stage 2).
catalog_diagnostics_dir <- function(stamp = pipeline_stamp) {
  file.path(pipeline_repo_root(), "diagnostics", stamp)
}

catalog_bulk_file <- function(bulk_name, stamp = pipeline_stamp) {
  root <- pipeline_repo_root()
  fp <- file.path(
    root, "data", "catalog_bulk", stamp,
    paste0("bulk_", bulk_name, "_", stamp, ".csv")
  )
  if (!file.exists(fp)) {
    stop("Catalog file not found: ", fp, call. = FALSE)
  }
  normalizePath(fp, mustWork = TRUE)
}
