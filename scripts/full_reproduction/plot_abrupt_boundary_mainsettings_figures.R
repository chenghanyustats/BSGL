source(file.path("R", "beta_functions.R"))
source(file.path("R", "helpers_gplasso.R"))

library(ggplot2)
library(viridis)

fit_path <- "results/abrupt_boundary/main_simulation_settings_rerun/abrupt_boundary_mainsettings_fit.rds"
out_dir <- "figures/reproduction_outputs/abrupt_boundary"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(fit_path)
grid_data <- obj$grid_data
meta <- obj$meta

est_betas <- betas_bsgl(obj$bsgl_fit, grid_data)
true_betas <- grid_data$true_betas
coords <- as.data.frame(grid_data$grid_points)
colnames(coords) <- c("x", "y")

plot_one <- function(j) {
  df_all <- rbind(
    data.frame(
      x = coords$x,
      y = coords$y,
      value = true_betas[, j],
      Type = "True"
    ),
    data.frame(
      x = coords$x,
      y = coords$y,
      value = est_betas[, j],
      Type = "BSGL"
    ),
    data.frame(
      x = coords$x,
      y = coords$y,
      value = est_betas[, j] - true_betas[, j],
      Type = "Difference"
    )
  )

  df_all$Type <- factor(df_all$Type, levels = c("True", "BSGL", "Difference"))

  p <- ggplot(df_all, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    facet_wrap(~ Type, nrow = 1) +
    coord_fixed(expand = FALSE) +
    scale_fill_viridis_c() +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme_void(base_size = 16) +
    theme(
      strip.text = element_text(face = "bold", size = 17),
      legend.text = element_text(size = 12),
      plot.margin = margin(2, 2, 2, 2)
    )

  out_file <- file.path(
    out_dir,
    sprintf("beta%d_BSGL_abrupt_mainsettings_n%d_p%d.png", j, meta$n, meta$p)
  )
  ggsave(out_file, p, width = 9.5, height = 2.8, dpi = 300, bg = "white")
  out_file
}

files <- vapply(c(2, 5), plot_one, character(1))
print(files)
