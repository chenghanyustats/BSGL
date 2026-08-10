# BSGL: Bayesian Spatial Group Lasso

This repository contains a reviewer-facing R package, prepared analysis data, result tables, and generated figures for the BSGL spatially varying coefficient paper submitted to JABES.

The reviewer workflows install the package, load `BSGL`, and refit models locally for the prepared simulation and MODIS EVI examples.

## Quick Start

Run from the repository root:

```sh
Rscript run_demo.R
Rscript run_real_demo.R
```

The simulation demo fits BSGL, GSVC, GGP-GAM, and MGWR on:

```text
data/simulation_demo/n1000_p5
```

It writes:

```text
results/demo_method_comparison.csv
```

Use the default simulation-demo MCMC settings for the expected comparison
table. Shorter iteration settings are useful only for checking that the script
runs end to end.

The MODIS EVI demo fits BSGL, GSVC, and GGP-GAM by default on the prepared real-data subset and writes:

```text
results/real_demo_method_comparison.csv
```

Both demo outputs contain the columns `Method`, `Test_MSE`, and `Test_MAE`. The real-data demo includes an MGWR row in the output table, but skips the local MGWR fit by default to keep the smoke test fast. Set `REAL_DEMO_RUN_MGWR=true` to attempt the optional local MGWR fit. The optional MGWR run uses a warm-start bandwidth initialization and skips automatically when the training set exceeds `REAL_DEMO_MGWR_MAX_TRAIN` (default: 1000).

The retained CSV files under `results/` provide the manuscript-reference values reported in the paper. The simulation demo is a single local run on the prepared \(n = 1000, m = 5\) data set; the manuscript-reference simulation results are the retained 10-replicate summaries in:

```text
results/simulation/rep10_metrics_by_replicate.csv
results/simulation/rep10_scp_by_replicate.csv
results/simulation/rep10_metrics_summary.csv
results/simulation/rep10_scp_summary.csv
```

Because the demo scripts refit the models locally, small numerical differences across local runs may occur due to runtime settings, stochastic fitting components, or seed/split differences.

The retained four-method MODIS EVI prediction table, including MGWR, is stored under:

```text
results/real_data
```

## Package Install

The package can also be installed manually from the repository root:

```sh
R CMD INSTALL .
```

After installation, the core functions can be used directly from R:

```r
library(BSGL)

fit <- fit_bsgl(y, X, coords)
pred <- pred_bsgl(fit, X_new, coords_new)
```

## Vignettes

Two package vignettes are included under `vignettes/`:

```text
vignettes/simulation_demo.Rmd
vignettes/modis_real_data_demo.Rmd
```

They provide compact worked examples for the prepared simulation data and the
cleaned MODIS EVI subset. The code chunks are written in runnable form and are
kept unevaluated during package installation so that installation does not
automatically run MCMC.

## Results and Figures

Paper summaries are available under:

```text
results/simulation
results/real_data
results/interaction
results/abrupt_boundary
```

Generated figures are available under:

```text
figures/main
figures/supplement
figures/reproduction_outputs
```

For compact checks of tables and figures, run:

```sh
Rscript scripts/02_summarize_main_results.R
Rscript scripts/04_figure_inventory.R
```

## Method Overview

Let `z_i = (u_i, v_i)` denote spatial coordinates. The model is:

```text
y(z_i) = sum_j x_j(z_i) beta_j(z_i) + epsilon_i,
epsilon_i ~ N(0, sigma^2).
```

Each coefficient surface is represented with a tensor-product B-spline basis. BSGL applies a Bayesian group lasso prior to the basis coefficients, shrinking negligible coefficient surfaces toward zero while retaining spatial structure for informative predictors.

Posterior inference provides coefficient-surface samples, pointwise credible intervals, credible-effect maps, and spatial coverage probabilities (SCPs).

## Main Functions

The primary BSGL functions are:

```r
bsgl_fit <- fit_bsgl(
  y = train_data$y,
  X = train_data$X,
  coords = train_data$coords,
  L = 25,
  niter = 1000,
  nburn = 100,
  a_lam = 15,
  b_lam = 1
)

bsgl_pred <- pred_bsgl(
  bsgl_fit,
  test_data$X,
  test_data$coords
)
```

The GSVC comparison functions are:

```r
gsvc_fit <- fit_gs(
  y = train_data$y,
  X = train_data$X,
  coords = train_data$coords,
  L = 25,
  niter = 1000,
  nburn = 100
)

gsvc_pred <- pred_gs(
  gsvc_fit,
  test_data$X,
  test_data$coords
)
```

## Repository Layout

```text
R/                 Core BSGL, GSVC, GGP-GAM, and MGWR helper functions
data/              Prepared simulation and analysis data
results/           Paper result tables
figures/           Code-generated paper figures
scripts/           Summary and reproduction scripts
vignettes/         Package vignettes
run_demo.R         Package-based simulation demo script
run_real_demo.R    Package-based MODIS EVI demo script
README.md          Repository overview
```

The paper uses `m` for the number of covariates. Some R code retains the conventional internal variable name `p`.

## Dependencies

The package dependencies are listed in `DESCRIPTION`. Main dependencies include `mvtnorm`, `MCMCpack`, `splines`, `GIGrvg`, `mgcv`, `GWmodel`, `sp`, `ggplot2`, `cowplot`, and `viridis`.

## Citation

Please cite the associated JABES paper:

Q. Zhan, C.-H. Yu, Y. Chen, Z. Dong, and R. Guhaniyogi, *Mapping Drivers of Greenness: Spatial Variable Selection for MODIS Vegetation Indices*.
