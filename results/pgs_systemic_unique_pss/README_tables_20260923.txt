PGS systemic portability — diagnostics (20260923)
eval_df: 17387 rows | 6773 PGS | 826 EFO

AUC-only rows (eval_auc): 5064
PGS with ≥3 display ancestries (pgs_keep): 294

Top 10 PGS by display-ancestry coverage:
# A tibble: 10 × 4
   pgs_id    n_ancestries_display n_eval_rows train_bucket          
   <chr>                    <int>       <int> <chr>                 
 1 PGS004859                    6          14 Train: Multi incl. EUR
 2 PGS000036                    6          12 Train: Multi incl. EUR
 3 PGS000013                    5          22 Train: Multi incl. EUR
 4 PGS001781                    5          15 Train: European-only  
 5 PGS005091                    5          14 Train: Multi incl. EUR
 6 PGS005092                    5          14 Train: Multi incl. EUR
 7 PGS002308                    5          12 Train: Multi incl. EUR
 8 PGS012555                    5          11 Train: Multi incl. EUR
 9 PGS000014                    5          10 Train: European-only  
10 PGS001357                    5          10 Train: European-only  

Evaluation counts per display ancestry (eval_work):
# A tibble: 6 × 2
  ancestry_display                    n
  <fct>                           <int>
1 European                          653
2 South Asian                       314
3 East Asian                        292
4 Middle Eastern or North African    80
5 Hispanic or Latin American         79
6 African                            46

Target pairs coverage (non-empty pairs only):
# A tibble: 6 × 2
  pair                                          n
  <chr>                                     <int>
1 Hispanic or Latin American vs European       50
2 Hispanic or Latin American vs East Asian     36
3 African vs European                          32
4 Hispanic or Latin American vs South Asian    31
5 African vs South Asian                       14
6 African vs East Asian                         9

Selected for barplots (n=10): PGS004859, PGS000036, PGS000013, PGS001781, PGS005091, PGS005092, PGS002308, PGS012555, PGS000014, PGS001357
