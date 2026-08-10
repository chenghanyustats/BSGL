# Figure Files

This directory contains generated figure outputs retained for the reviewer package.

- `main/`: figures used in the main paper.
- `supplement/`: figures used in the supplementary material.
- `reproduction_outputs/`: additional generated outputs retained for traceability.

The source script for each retained figure is listed in `figure_manifest.csv`.

Run the figure inventory check from the repository root:

```sh
Rscript scripts/04_figure_inventory.R
```
