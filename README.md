# BSGL: Bayesian Spatial Group Lasso

This repository contains the reviewer-facing R implementation, prepared analysis data, retained result tables, and generated figures for the BSGL spatially varying coefficient paper submitted to the *Journal of Agricultural, Biological, and Environmental Statistics* (JABES).

The demonstration workflows are designed to make the proposed BSGL method and the comparison methods easy to run from a standard R or RStudio session. The demo scripts install the BSGL package and required R-package dependencies when needed, obtain only the prepared data required for the selected example, fit the available methods independently, and write a common comparison table.

## Quick Start in R

No terminal commands are required for the demonstration workflows.

### Simulation demo

If `run_demo.R` is available in the current working directory, run:

```r
source("run_demo.R")
```

Alternatively, the current GitHub version can be sourced directly from R:

```r
source(
  "https://raw.githubusercontent.com/Qishi7/BSGL/main/run_demo.R"
)
```

If BSGL is not already installed, the script installs `remotes` when necessary and then installs BSGL directly from GitHub together with its required R-package dependencies.

By default, the simulation demo attempts all four methods on the same prepared simulation data set:

- BSGL (proposed method)
- GSVC
- GGP-GAM
- MGWR

The simulation demo uses the prepared \(n = 1000, m = 5\) example and writes:

```text
results/demo_method_comparison.csv
```

To skip MGWR in the simulation demo:

```r
Sys.setenv(DEMO_RUN_MGWR = "false")
source("run_demo.R")
```

To return to the default behavior in a later R session or run:

```r
Sys.unsetenv("DEMO_RUN_MGWR")
```

The default simulation-demo MCMC settings are intended for the reviewer comparison. Environment variables such as `DEMO_NITER` and `DEMO_NBURN` can be changed for shorter software checks, but abbreviated runs are not intended to reproduce the retained manuscript-reference results.

### MODIS EVI real-data demo

If `run_real_demo.R` is available in the current working directory, run:

```r
source("run_real_demo.R")
```

Alternatively:

```r
source(
  "https://raw.githubusercontent.com/Qishi7/BSGL/main/run_real_demo.R"
)
```

The real-data demo fits BSGL, GSVC, and GGP-GAM by default on the prepared MODIS EVI subset.

MGWR is **not run by default** for the real-data demo because `GWmodel::gwr.multiscale()` can be substantially more computationally intensive for this example. The MGWR row is still retained in the comparison table with status `Not run`.

To attempt the MGWR fit:

```r
Sys.setenv(REAL_DEMO_RUN_MGWR = "true")
source("run_real_demo.R")
```

The real-data workflow also uses `REAL_DEMO_MGWR_MAX_TRAIN` as a safeguard against unintentionally launching a very large MGWR fit. If MGWR is requested but the training sample exceeds this threshold, its status is reported as `Not run` with an explanatory message.

The real-data demo writes:

```text
results/real_demo_method_comparison.csv
```

To return to the default real-data behavior:

```r
Sys.unsetenv("REAL_DEMO_RUN_MGWR")
```

## Interpretation of Demo Output

Both demo scripts produce a common comparison table with the columns:

```text
Method, Status, Test_MSE, Test_MAE, Message
```

`Status` has three possible values:

- `Success`: the method completed and produced finite predictions and comparison metrics.
- `Failed`: the method was attempted but generated an error or invalid predictions. The original error information is recorded in `Message`.
- `Not run`: the method was deliberately skipped by the current configuration.

A failure of one method does **not** stop the remaining methods from being attempted. This makes it possible to distinguish a method-specific software or numerical problem from a failure of the complete demonstration workflow.

For `Failed` or `Not run` methods, unavailable MSE and MAE values are left blank in the written CSV file rather than displayed as the text `NA`.

At completion, the scripts also report whether all attempted methods succeeded or whether any method-specific failures occurred.

## Optional Manual Package Installation

The demonstration scripts install BSGL automatically if necessary. Users who prefer to install the package first can do so entirely within R:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("Qishi7/BSGL")
library(BSGL)
```

After installation, the core functions can be called directly:

```r
fit <- fit_bsgl(y, X, coords)
pred <- pred_bsgl(fit, X_new, coords_new)
```

To force a fresh installation of the current GitHub version before running either demo:

```r
Sys.setenv(BSGL_DEMO_INSTALL = "true")
source("run_demo.R")
```

or

```r
Sys.setenv(BSGL_DEMO_INSTALL = "true")
source("run_real_demo.R")
```

Afterward, the setting can be cleared with:

```r
Sys.unsetenv("BSGL_DEMO_INSTALL")
```

## Prepared Demo Data

The standalone demo scripts do not require the full repository to be cloned before execution.

Each script first checks whether its prepared demo data are available within the installed BSGL package. If not, it downloads only the data file or files required for that demonstration from this GitHub repository into a temporary directory.

The simulation demo uses the prepared training and test data corresponding to:

```text
data/simulation_demo/n1000_p5/
```

The real-data demo uses the prepared MODIS EVI data file:

```text
data/data_cleaned_small_expanded.RData
```

## Manuscript-Reference Results

The demo scripts refit the models locally and are intended to provide a transparent, runnable comparison of the proposed and competing methods.

The retained manuscript-reference simulation results are stored in:

```text
results/simulation/rep10_metrics_by_replicate.csv
results/simulation/rep10_scp_by_replicate.csv
results/simulation/rep10_metrics_summary.csv
results/simulation/rep10_scp_summary.csv
```

The retained four-method MODIS EVI prediction results, including MGWR, are stored under:

```text
results/real_data/
```

Because the demo scripts refit stochastic models locally, small numerical differences can occur across software environments or runs. The retained files under `results/` provide the reference values used for the manuscript and supplementary material.

## Vignettes

Two package vignettes are included under:

```text
vignettes/simulation_demo.Rmd
vignettes/modis_real_data_demo.Rmd
```

They provide compact worked examples for the prepared simulation data and the cleaned MODIS EVI subset. The code chunks are written in runnable form and are kept unevaluated during package installation so that package installation does not automatically launch MCMC fitting.

## Results and Figures

Paper summaries are available under:

```text
results/simulation/
results/real_data/
results/interaction/
results/abrupt_boundary/
```

Generated figures are available under:

```text
figures/main/
figures/supplement/
figures/reproduction_outputs/
```

If the full repository has been downloaded, compact checks of the retained tables and figure inventory can be run from R with:

```r
source("scripts/02_summarize_main_results.R")
source("scripts/04_figure_inventory.R")
```

## Method Overview

Let \(z_i = (u_i, v_i)\) denote spatial coordinates. The spatially varying coefficient model is

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

The package also includes wrappers used by the reviewer demos for GGP-GAM and MGWR:

```r
gam_fit <- fit_gam(X, y, coords)
gam_pred <- pred_gam(gam_fit, X_new, coords_new)

mgwr_fit <- fit_mgwr(X, y, coords)
mgwr_pred <- pred_mgwr(mgwr_fit, X_new, coords_new)
```

## Repository Layout

```text
R/                    Core BSGL, GSVC, GGP-GAM, and MGWR functions
data/                 Prepared simulation and analysis data
results/              Retained result tables
figures/              Code-generated paper figures
scripts/              Summary and reproduction scripts
vignettes/            Package vignettes
run_demo.R            Standalone simulation comparison demo
run_real_demo.R       Standalone MODIS EVI comparison demo
README.md             Repository overview and reviewer instructions
```

The paper uses `m` for the number of covariates. Some R code retains the conventional internal variable name `p`.

## Dependencies

BSGL requires R >= 4.1.0. Package dependencies are declared in `DESCRIPTION` and are installed through the standard R package installation process.

Main dependencies currently include `GIGrvg`, `GWmodel`, `MCMCpack`, `mgcv`, `mvtnorm`, `sf`, `sp`, `spdep`, `ggplot2`, `gridExtra`, `cowplot`, `rnaturalearth`, and `viridis`, among others listed in `DESCRIPTION`.

Some spatial R packages may require platform-specific system libraries. If installation of an R dependency fails, the underlying installation error should be resolved before rerunning the affected method.

## Citation

Please cite the associated JABES paper:

Q. Zhan, C.-H. Yu, Y. Chen, Z. Dong, and R. Guhaniyogi, *Mapping Drivers of Greenness: Spatial Variable Selection for MODIS Vegetation Indices*.
