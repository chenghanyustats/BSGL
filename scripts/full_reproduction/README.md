# Full Reproduction Scripts

These are the longer source-repository workflows used to reproduce retained
simulation, real-data, interaction, and figure results. They are not installed
with the R package and should be run from the **repository root**.

Before a long run, execute:

```r
source("scripts/01_repository_preflight.R")
```

The quick checks `scripts/02_summarize_main_results.R` and
`scripts/04_figure_inventory.R` do not refit the models.

## Data locations

- `data/data_rep10/`: ten-replicate simulation objects.
- `data/data_indep/`: independent simulation data sets.
- `data/data_cleaned_small_expanded.RData`: prepared MODIS reviewer subset.
- `data/geo_naturalearth_10m/`: retained geographic layers used by map scripts.

These larger inputs are intentionally excluded from the built R package by
`.Rbuildignore`; use a GitHub source checkout for full reproduction.

## Recommended execution order

### Core simulation

1. `sim_10.R` — full BSGL/GGP-GAM/GSVC/MGWR simulation comparison.
2. `rerun_mgwr.R` — optional MGWR completion/recovery for the same retained grid.
3. `sim_gpsvc.R` — GSVC-only runs on the retained `data/data_indep` paper grid.

### Four-method n=1000, m=5 figure workflow

1. `plot_four_method_gray_n1000_p5.R` — fits methods and writes the retained
   `gray_run_fit.rds` under `figures/reproduction_outputs/simulation/`.
2. `plot_beta5_ownscale_onebar.R` — consumes that fit object.

### Abrupt-boundary workflow

1. `run_abrupt_boundary_mainsettings.R` — generates the fit/result object.
2. `plot_abrupt_boundary_mainsettings_figures.R` — consumes the generated fit
   and writes reproduction figures.

### Simulated interaction workflow

1. `interaction_pilot_main_effects.R` — produces
   `results/interaction/main_effects_rerun/main_interaction_pilot_fit.rds`.
2. Then run any of:
   - `interaction_scp_metrics.R`
   - `plot_interaction_beta_recovery.R`
   - `plot_interaction_x1_contribution.R`

### MODIS no-interaction workflow

1. `real_nointeraction_scp_compare.R` — produces the BSGL/GSVC fit and SCP files.
2. Then run as needed:
   - `organize_real_10variable_results.R`
   - `plot_real_10variable_scp_maps_mean_lc.R`
   - `plot_test_evi_predictive_quantiles.R`
3. `real_prediction_four_method_compare.R` writes the retained four-method
   manuscript prediction summary; it does not refit all four methods.

### MODIS interaction workflow

1. `real_interaction_evi_validation.R` — produces the interaction fit object.
2. `plot_real_interaction_scp_maps_1x2.R` — consumes the retained/generated fit.

### Main Figure 2

`regenerate_main_fig2_sim_scp_map.R` is self-contained apart from the retained
`data/data_rep10/n5000_p10.RData` input and writes to `figures/main/`.

## Computational note

Several workflows launch long MCMC or MGWR fits. The repository preflight checks
paths, syntax, retained dependencies, quick summaries, and vignette rendering
without automatically launching those long fits.

## Reproducibility validation

Before launching the computationally intensive full reproduction, the repository can be validated in three stages. Run the following commands from the repository root:

```r
source("scripts/01_repository_preflight.R")
source("scripts/05_runtime_smoke_tests.R")
source("scripts/06_remaining_runtime_smoke_tests.R")
```

The repository preflight checks required code and data files, retained result objects, the paper-figure inventory, script parsing, and, optionally, vignette rendering. The two runtime smoke-test scripts then exercise the major producer and downstream/post-processing execution paths using shortened computational settings and temporary output directories.

`05_runtime_smoke_tests.R` and `06_remaining_runtime_smoke_tests.R` should be run sequentially in the same R session because the latter reuses temporary fit objects created by the former.

These smoke tests validate execution paths and dependencies. They do **not** constitute a full rerun of all manuscript analyses at the original simulation, MCMC, and other computational settings.
