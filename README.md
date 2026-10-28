# CovSF-NS — stability-based ensemble covariate selection for non-stationary flood frequency analysis

R code that reproduces every number, table and figure in the manuscript. Tested with R 4.3.3.

## Packages
`readxl, gamlss, gamlss.dist, gamboostLSS, mboost, stabs, modifiedmk, lmom, parallel`;
for the ensemble-learner (ELC) baseline also `glmnet, mgcv, RRF, xgboost, data.table, stringi`.

## Layout
Run every script from this folder, e.g. `Rscript R/01_screening.R`. Inputs are read from `data/` and outputs written to `out/` (override with the `COVSF_DATA`, `COVSF_OUT` and `COVSF_R` environment variables).

## Data (`data/`)
| file | content |
|---|---|
| `narmada_flow.xlsx` | annual maximum discharge, 1980–2019, one column per gauge |
| `annual_rainfall_narmada.xlsx`, `DMR_narmada.xlsx` | catchment annual rainfall, annual max daily rainfall |
| `Narmada_cov_multi.xlsx` | 15 climate indices (one value per year) |
| `MRI_narmada.xlsx` | original reservoir index (year-specific runoff); used **only** to demonstrate leakage and in the ELC baseline |
| `dams_narmada.csv` *(recommended; start from `dams_narmada_TEMPLATE.csv`)* | `dam, A_km2, C_bcm, start` — catchment area, gross storage, first year of operation. If absent, the public-source table in `00_setup.R` is used. Upstream membership is set in `UPSTREAM` in `00_setup.R`. |

## Main workflow (manuscript version: stabilised four-learner ensemble)
Run in this order: `01_screening.R` → `03_distribution.R` → `04_selection.R` (all-subsets table) → `04b_cluster_stability.R` (needed for the clustering file and as an alternative selector) → `11_ensemble_stability.R` (**CovSF-NS**: GLM elastic net, GAM null-space, guided RF, XGBoost consensus repeated on 100 half-samples) → copy `out/ens_final_models.csv` merged with `out/sel_unconstrained_best.csv` to `out/sel_final_models.csv` → `14_original_style_outputs.R`, `07_benchmarks_returnlevels.R`, `08_extra_checks.R`, `09_rl_sensitivity.R` → `05_montecarlo.R`, `12_mc_ensemble.R` → `10_nested_loo.R`, `13_nested_ensemble.R`.
The ensemble learners and the stability layer are in `ens_functions.R`.

## All scripts (outputs go to `out/`, or set `COVSF_OUT`)
| script | stage | main outputs |
|---|---|---|
| `00_setup.R` | data, leak-free reservoir index, candidate frames (sourced by all) | – |
| `01_screening.R` | Stage 0 trend/change-point tests with BH + field significance; Stage 1 leakage check, lag screening, redundancy clustering | `stage0_*.csv`, `stage1_*.csv/rds` |
| `02_engine.R` | likelihood engine (WEI3, GA, LOGNO, RG, IG, EXP, GG, GEV-GML), LOO scores, bootstrap GOF (sourced) | – |
| `03_distribution.R` | regional distribution: summed AICc, L-moments, KS/AD | `dist_*.csv`, `lmoment_ratios.csv` |
| `04_selection.R` | exhaustive all-subsets search (k ≤ 3), AICc and LR tests | `sel_level3_exhaustive.csv`, `sel_unconstrained_best.csv` |
| `04b_cluster_stability.R` | Stage 2a regional + 2b station cluster stability selection; final models | `sel_level1_regional.csv`, `sel_level2_*.csv`, `sel_final_models.csv` |
| `05_montecarlo.R` | Monte Carlo recovery experiment | `mc_raw_*.csv`, `mc_summary_*.csv` |
| `06_v1_sensitivity.R` | learner-subset and bootstrap sensitivity of the ELC baseline | `v1_*.csv` |
| `07_benchmarks_returnlevels.R` | benchmarks vs competing NS models; return levels with bootstrap CIs; GOF | `benchmarks.csv`, `returnlevels.csv`, `effective_rl100.csv`, `final_*.csv` |
| `08_extra_checks.R` | LR test for covariates in the shape parameter; GEV shapes | `extra_sigma_LR_gevshape.csv` |
| `09_rl_sensitivity.R` | return-level sensitivity to the distribution family | `rl_sensitivity.csv` |
| `10_nested_loo.R` | nested leave-one-out (selection repeated in every fold) | `nested_loo.csv` |
| `ens_functions.R` | the four learners (continuous target), original intersection/union consensus, stability layer (sourced) | – |
| `11_ensemble_stability.R` | stabilised ensemble selection, single-learner comparison, final models | `ens_*.csv` |
| `12_mc_ensemble.R` | Monte Carlo for one-shot vs stabilised ensemble (same series as 05) | `mc_raw_ensemble.csv` |
| `13_nested_ensemble.R` | nested LOO for the ensemble | `nested_loo_with_ensemble.csv` |
| `14_original_style_outputs.R` | quantities for Figs 6–8 and Table 2 (LOO PR curves, LOO MAE/NSE, AIC) | `orig_style_*.csv` |
| `stab_functions.R` | fast likelihood boosting + complementary-pairs cluster stability selection (sourced) | – |
| `v1_embedded_selection.R` | ELC baseline (GLM-elastic net, GAM, RF, XGBoost consensus) | – |
| `check_engine.R`, `check_booster.R` (validation, run from this folder) | validation: engine vs `gamlss()`; fast booster vs `gamboostLSS` | console |

Approximate run times on 2 cores: 01–04b < 5 min; 07 ≈ 15 min; 05 ≈ 30 min; 06 ≈ 30 min; 10 ≈ 40 min.

## Key settings
`FAMILY = "WEI3"` (regional distribution), stability cutoff `0.75`, `PFER = 1`, `B = 50` complementary pairs,
max 3 covariates, families: Rainfall {AR, DMR}, Trend/Regulation {Year, RI}, plus one representative per climate-index cluster.
Rainfall covariates enter as log(AR), log(DMR) (power law under the log link).
