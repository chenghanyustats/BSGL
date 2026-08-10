library(ggplot2)

in_dir <- "results/real_data/nointeraction_scp_compare_n10000_rerun"
out_dir <- "results/real_data/clean_real_data_10variable_rerun"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

scp_raw <- read.csv(file.path(in_dir, "real_nointeraction_scp_by_column.csv"))
pred <- read.csv(file.path(in_dir, "real_nointeraction_prediction_summary.csv"))
diagnostics_file <- file.path(out_dir, "real_prediction_coverage_diagnostics.csv")
if (file.exists(diagnostics_file)) {
  pred_diag <- read.csv(diagnostics_file)
  pred <- pred_diag[, c("Method", "MSPE", "Predictive_Coverage", "Latent_CI_Coverage")]
} else {
  names(pred)[names(pred) == "Coverage"] <- "Predictive_Coverage"
}

continuous_order <- c(
  "red_reflectance",
  "NIR_reflectance",
  "blue_reflectance",
  "MIR_reflectance",
  "GPP",
  "LE",
  "view_zenith_angle",
  "sun_zenith_angle",
  "relative_azimuth_angle"
)

labels <- c(
  red_reflectance = "Red Reflectance",
  NIR_reflectance = "NIR Reflectance",
  blue_reflectance = "Blue Reflectance",
  MIR_reflectance = "MIR Reflectance",
  GPP = "GPP",
  LE = "LE",
  view_zenith_angle = "View Zenith Angle",
  sun_zenith_angle = "Sun Zenith Angle",
  relative_azimuth_angle = "Relative Azimuth Angle",
  LC_Type4 = "LC Type4"
)

collapse_one_method <- function(df) {
  cont <- df[df$variable %in% continuous_order, c("Method", "variable", "label", "SCP")]
  lc <- df[grepl("^LC_Type4_", df$variable), ]
  lc_row <- data.frame(
    Method = unique(df$Method),
    variable = "LC_Type4",
    label = "LC Type4",
    SCP = mean(lc$SCP, na.rm = TRUE)
  )
  out <- rbind(cont, lc_row)
  out$variable <- factor(out$variable, levels = c(continuous_order, "LC_Type4"))
  out <- out[order(out$variable), ]
  out$variable <- as.character(out$variable)
  out$label <- unname(labels[out$variable])
  out
}

scp_10 <- do.call(
  rbind,
  lapply(split(scp_raw, scp_raw$Method), collapse_one_method)
)

write.csv(pred, file.path(out_dir, "real_prediction_bsgl_gaussian_svc.csv"), row.names = FALSE)
write.csv(scp_10, file.path(out_dir, "real_scp_10variable_bsgl_gaussian_svc.csv"), row.names = FALSE)

wide <- reshape(
  scp_10[, c("variable", "label", "Method", "SCP")],
  idvar = c("variable", "label"),
  timevar = "Method",
  direction = "wide"
)
names(wide) <- sub("^SCP\\.", "", names(wide))
wide$variable <- factor(wide$variable, levels = c(continuous_order, "LC_Type4"))
wide <- wide[order(wide$variable), ]
wide$variable <- as.character(wide$variable)
write.csv(wide, file.path(out_dir, "real_scp_10variable_wide.csv"), row.names = FALSE)

plot_df <- scp_10
plot_df$label <- factor(plot_df$label, levels = rev(unname(labels[c(continuous_order, "LC_Type4")])))

p <- ggplot(plot_df, aes(x = label, y = SCP, fill = Method)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  coord_flip() +
  scale_fill_manual(values = c("BSGL" = "gray25", "Gaussian SVC" = "gray70")) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(y = "SCP")

ggsave(file.path(out_dir, "real_scp_10variable_bsgl_gaussian_svc.png"), p, width = 7.5, height = 5.5, dpi = 300)

latex_lines <- c(
  "\\begin{table}[!t]",
  "\\centering",
  "\\caption{Real-data comparison using the MODIS EVI data.}",
  "\\label{tab:real_data_gaussian_svc_comparison}",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Variable & BSGL SCP & Gaussian SVC SCP \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(wide))) {
  latex_lines <- c(
    latex_lines,
    sprintf(
      "%s & %.3f & %.3f \\\\",
      wide$label[i],
      wide$BSGL[i],
      wide$`Gaussian SVC`[i]
    )
  )
}

latex_lines <- c(
  latex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(latex_lines, file.path(out_dir, "real_scp_10variable_table.tex"))

cat("Saved clean real-data results to:", out_dir, "\n")
print(wide)
cat("\nPrediction summary\n")
print(pred)
