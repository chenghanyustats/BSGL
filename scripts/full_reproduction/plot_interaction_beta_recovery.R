rm(list = ls())

source(file.path("R", "helpers_gplasso.R"))
source(file.path("R", "helpers_gpprior.R"))
source(file.path("R", "helpers_gam.R"))
source(file.path("R", "helpers_mgwr.R"))

library(ggplot2)

fit_path <- file.path(
  "results/interaction/main_effects_rerun",
  "main_interaction_pilot_fit.rds"
)
out_dir <- "results/interaction/main_effects_rerun"

obj <- readRDS(fit_path)
dat <- obj$data
grid_data <- dat$grid
grid <- as.data.frame(grid_data$grid_points)
colnames(grid) <- c("x", "y")

var_names <- colnames(dat$train$X)

bsgl_beta <- betas_bsgl(obj$bsgl_fit, grid_data)
gs_beta <- betas_gs(obj$gs_fit, grid_data)
gam_beta <- get_gam_betas(obj$gam_fit, as.matrix(grid_data$grid_points))
mgwr_beta <- get_mgwr_betas(obj$mgwr_fit, as.matrix(grid_data$grid_points))

colnames(bsgl_beta) <- var_names
colnames(gs_beta) <- var_names
colnames(gam_beta) <- var_names
colnames(mgwr_beta) <- var_names

plot_one_beta <- function(j, file_stub, title_text) {
  plot_df <- rbind(
    data.frame(grid, method = "True", beta = grid_data$true_betas[, j]),
    data.frame(grid, method = "BSGL", beta = bsgl_beta[, j]),
    data.frame(grid, method = "Gaussian SVC", beta = gs_beta[, j]),
    data.frame(grid, method = "GGP-GAM", beta = gam_beta[, j]),
    data.frame(grid, method = "MGWR", beta = mgwr_beta[, j])
  )

  method_levels <- c("True", "BSGL", "Gaussian SVC", "GGP-GAM", "MGWR")
  plot_df$method <- factor(plot_df$method, levels = method_levels)

  p <- ggplot(plot_df, aes(x, y, fill = beta)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_viridis_c(option = "D") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 12),
      plot.title = element_text(size = 13, face = "bold")
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = NULL,
      fill = "beta"
    )

  out_file <- file.path(out_dir, paste0(file_stub, ".png"))
  ggsave(out_file, p, width = 13, height = 3.2, dpi = 240)
  cat("Saved:", normalizePath(out_file), "\n")
}

plot_one_beta(
  j = 6,
  file_stub = "beta6_x1x2_recovery",
  title_text = "Beta 6: centered scaled X1:X2"
)

plot_one_beta(
  j = 7,
  file_stub = "beta7_x1x5_recovery",
  title_text = "Beta 7: centered scaled X1:X5"
)
