# pgscatalog_meta

Code for systematic evaluation of polygenic score (PGS) portability across ancestries using data from the PGS Catalog.

## Overview

This repository contains the analysis pipeline used to:

- Download and process PGS Catalog performance data

- Construct evaluation-level datasets across ~3,900 PGS evaluations

- Perform inverse-variance weighted meta-analysis of AUC

- Quantify portability as ΔAUC (target ancestry − European)

- Generate publication figures, including faceted ΔAUC forest plots

## Data source

All data are derived from the PGS Catalog using the `quincunx` R package.

## Repository structure

analysis/
01_Installation_API_and_perfomance_pull.R
02_loading_PGS_matrix.R
03_bulk_metada_pull.R
06_PGS_systemic_portability_unique_pss.R
07_2_pgs_auc_ci_audit.R
12_meta_roadmap_two_stages.R
12_i_square_filter.R
12_meta_results_overview_pythonStyle.R

All scripts are stored in the `analysis/` folder and are ordered chronologically according to the development of the analysis.

## Main analysis pipeline

To reproduce the publication figure from raw data, run the scripts in the following order:

1. **01_Installation_API_and_perfomance_pull.R**  
   Downloads performance metrics from the PGS Catalog API and stores normalized tables in `data/pgs_cache`.

2. **02_loading_PGS_matrix.R**  
   Builds and validates relationships between performance metrics, sample sets, and ancestry.

3. **03_bulk_metada_pull.R**  
   Retrieves score metadata and constructs training ancestry categories (`train_bucket`).

4. **06_PGS_systemic_portability_unique_pss.R**  
   Builds the evaluation-level dataset with unique sample sets.

5. **07_2_pgs_auc_ci_audit.R**  
   Harmonizes AUC values and confidence intervals and produces the final evaluation table.

6. **12_meta_roadmap_two_stages.R**  
   Performs Stage-1 meta-analysis by pooling AUC across cohorts for each PGS and ancestry using inverse-variance weighting.

7. **12_i_square_filter.R**  
   Filters pooled estimates based on heterogeneity (I² threshold).

8. **12_meta_results_overview_pythonStyle.R**  
   Generates publication figures, including the faceted ΔAUC forest plot.

## Key output

The main publication figure is a **faceted ΔAUC forest plot**, which shows:

- Differences in predictive performance across ancestries
- ΔAUC relative to European evaluation cohorts
- Aggregated results across PGS using inverse-variance weighting

## Notes on reproducibility

- Scripts are designed to be run sequentially.
- Intermediate results are written to disk and reused in later steps.
- Running from a fresh environment requires starting from step 1.
- Some scripts assume the presence of previously generated files in `data/` and `results/`.

## Requirements

R (≥ 4.0 recommended)

Main packages used:

- quincunx  
- dplyr  
- tidyr  
- readr  
- ggplot2  
- stringr  
- purrr  
- forcats  

## License

See `LICENSE` file.
