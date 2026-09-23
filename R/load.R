# Source the shared pipeline functions. Run analysis scripts from the repository root.

pipeline_r_dir <- (function() {
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile)))
  }
  if (file.exists(file.path("R", "ci_parse.R"))) return(normalizePath("R"))
  stop("Cannot locate the R/ directory. Run scripts from the repository root.")
})()

source(file.path(pipeline_r_dir, "sanity.R"))
source(file.path(pipeline_r_dir, "ci_parse.R"))
source(file.path(pipeline_r_dir, "ancestry.R"))
source(file.path(pipeline_r_dir, "traits.R"))
source(file.path(pipeline_r_dir, "stage1_ivw.R"))
source(file.path(pipeline_r_dir, "stage2_pair.R"))
