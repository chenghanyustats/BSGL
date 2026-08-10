# Full Reproduction Scripts

These scripts are the cleaned versions of the longer analysis and figure-generation scripts used during the paper preparation. Run them from the repository root directory.

The quick entry points `../02_summarize_main_results.R` and `../04_figure_inventory.R` are intended for fast checks. The package-based result script is `../../run_demo.R`. The scripts in this folder are longer analysis scripts for the retained paper results and figures.

Important scripts:

- `sim_10.R`: full simulation comparison over the paper grid.
- `rerun_mgwr.R`: MGWR simulation runner for the retained comparison grid.
- `sim_gpsvc.R`: GSVC simulation runner for the retained comparison grid.
- `plot_four_method_gray_n1000_p5.R`: four-method beta-surface figure generation for n = 1000, m = 5.
- `real_prediction_four_method_compare.R`: MODIS four-method prediction comparison.
- `real_nointeraction_scp_compare.R`: MODIS ten-variable BSGL vs GSVC comparison.
- `real_interaction_evi_validation.R`: MODIS interaction validation.
- `plot_real_10variable_scp_maps_mean_lc.R`: ten-variable real-data SCP map.
- `run_abrupt_boundary_mainsettings.R`: abrupt-boundary simulation.

Simulation data are stored in `data/data_rep10` and `data/data_indep`. The cleaned MODIS subset is stored in `data/data_cleaned_small_expanded.RData`.
