# Check retained generated paper figures.

source("scripts/00_project_setup.R")

required_figures <- c(
  "figures/main/fig1a_beta3.png",
  "figures/main/fig1b_beta5.png",
  "figures/main/fig2_simulation_scp_map.png",
  "figures/main/fig3_evi_comparison.png",
  "figures/main/fig4_modis_scp_maps.png",
  "figures/supplement/figS1_trace_plots.pdf",
  "figures/supplement/figS2_parameter_dist.png",
  "figures/supplement/figS3_abrupt_beta2.png",
  "figures/supplement/figS4_abrupt_beta5.png",
  "figures/supplement/figS5_prediction_errors.png",
  "figures/supplement/figS6_beta12_surface.png",
  "figures/supplement/figS7_beta12_scp.png",
  "figures/supplement/figS8_x1_contribution.png",
  "figures/supplement/figS9_real_interaction_scp.png"
)

status <- data.frame(
  figure = required_figures,
  exists = file.exists(project_path(required_figures)),
  size_kb = ifelse(
    file.exists(project_path(required_figures)),
    round(file.info(project_path(required_figures))$size / 1024, 1),
    NA_real_
  )
)

print(status)

if (!all(status$exists)) {
  stop("Some retained generated figures are missing.")
}
message("All retained generated paper figures are present.")
