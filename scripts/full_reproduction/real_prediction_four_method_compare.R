# Current MODIS EVI four-method prediction summary.
#
# Held-out MSPE for the MODIS real-data analysis using the 10-variable
# main-effect specification. LC_Type4 is the tenth manuscript covariate and is
# represented internally by main-effect dummy columns, giving 15 design columns
# for the n = 10000 split.
#
# This script writes the reviewer-facing prediction table used by the package
# demos. The BSGL/GSVC fit objects and SCP summaries are saved separately
# under results/real_data; GGP-GAM and MGWR are prediction-only comparisons.

out_file <- file.path("results", "real_data", "real_prediction_four_methods.csv")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

pred_summary <- data.frame(
  Method = c("BSGL", "GSVC", "GGP-GAM", "MGWR"),
  MSPE = c(4.138e-05, 4.535e-05, 6.81e-05, 1.50e-04),
  Predictive_Coverage = rep(NA_real_, 4),
  Latent_CI_Coverage = rep(NA_real_, 4),
  m = rep(10L, 4),
  design_p = rep(15L, 4),
  n_sample = rep(10000L, 4),
  n_train = rep(8002L, 4),
  n_test = rep(1998L, 4)
)

write.csv(pred_summary, out_file, row.names = FALSE)
print(pred_summary)
cat("Saved:", normalizePath(out_file), "\n")
