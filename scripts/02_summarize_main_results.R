# Summarize the already consolidated paper results.

source("scripts/00_project_setup.R")

sim_metrics <- read.csv(project_path("results", "simulation", "rep10_metrics_summary.csv"))
sim_scp <- read.csv(project_path("results", "simulation", "rep10_scp_summary.csv"))
real_pred <- read.csv(project_path("results", "real_data", "real_prediction_four_methods.csv"))
real_scp <- read.csv(project_path("results", "real_data", "real_scp_10variable_wide.csv"))

cat("\nSimulation prediction summary: n, m, Method, MSPE\n")
print(sim_metrics[, intersect(c("n", "m", "Method", "MSPE_mean", "MSPE_sd", "Coverage_mean"), names(sim_metrics))])

cat("\nSimulation SCP summary: first rows\n")
print(head(sim_scp, 20))

cat("\nMODIS EVI prediction summary\n")
print(real_pred)

cat("\nMODIS EVI ten-variable SCP summary\n")
print(real_scp)
