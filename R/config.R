# Named input for this revision. Scripts read this stamp instead of the newest file.

pipeline_stamp <- "20260923"

catalog_bulk_file <- function(bulk_name, stamp = pipeline_stamp) {
  root <- if (exists("pipeline_r_dir", inherits = TRUE)) dirname(pipeline_r_dir) else "."
  fp <- file.path(
    root, "data", "catalog_bulk", stamp,
    paste0("bulk_", bulk_name, "_", stamp, ".csv")
  )
  if (!file.exists(fp)) {
    stop("Catalog file not found: ", fp, call. = FALSE)
  }
  normalizePath(fp, mustWork = TRUE)
}
