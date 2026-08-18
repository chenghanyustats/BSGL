# Prepared Simulation Data

This directory contains the prepared simulation subset used by the reviewer-facing
reviewer result script.

## Data Specifications

- `n = 1000`
- `m = 5`
- `sigma2 = 0.1`
- independent predictors
- seed `123`

## Files

- `train_data.RData`: training response, covariates, and coordinates
- `test_data.RData`: held-out response, covariates, and coordinates
- `grid_data.RData`: grid points and true coefficient surfaces
- `meta_info.RData`: simulation metadata

The main result script loads these files through:

```r
source("run_demo.R")
```
