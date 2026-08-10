source(file.path("R", "beta_functions.R"))
source(file.path("R", "helpers_ggpgam.R"))
source(file.path("R", "helpers_gplasso.R"))
source(file.path("R", "helpers_gpprior.R"))
source(file.path("R", "helpers_mgwr.R"))
source(file.path("R", "calc_beta.R"))

library(ggplot2)
library(cowplot)
library(grid)

fit_file <- "figures/simulation/gray_four_method_n1000_p5_rerun/gray_run_fit.rds"
out_dir <- "figures/simulation/gray_four_method_n1000_p5_rerun"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(fit_file)
grid_data <- obj$grid_data
meta_info <- obj$meta_info
grid_points <- as.data.frame(grid_data$grid_points)
colnames(grid_points) <- c("x", "y")

j <- 5
beta_surfaces <- list(
  "True" = calc_true_beta(as.matrix(grid_points), meta_info$p)[, j],
  "BSGL" = betas_bsgl(obj$bsgl_fit, grid_data)[, j],
  "GGP-GAM" = get_gam_betas(obj$gam_model, as.matrix(grid_points))[, j],
  "GSVC" = betas_gs(obj$gs_fit, grid_data)[, j],
  "MGWR" = get_mgwr_betas(obj$mgwr_model, as.matrix(grid_points))[, j]
)

plot_data <- lapply(beta_surfaces, function(z) {
  data.frame(x = grid_points$x, y = grid_points$y, beta = z)
})

panel_lims <- lapply(beta_surfaces, function(z) {
  r <- range(z, na.rm = TRUE)
  if (!all(is.finite(r))) r <- c(-1, 1)
  if (diff(r) < 1e-8) {
    center <- mean(r)
    r <- center + c(-1, 1) * 1e-6
  }
  r
})

make_panel <- function(df, ttl, lims, keep_legend = FALSE) {
  p <- ggplot(df, aes(x = x, y = y, fill = beta)) +
    geom_raster() +
    coord_fixed() +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.margin = margin(1, 1, 1, 1),
      panel.border = element_blank(),
      legend.title = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0, "pt")
    ) +
    ggtitle(ttl) +
    scale_fill_gradient(
      low = "white",
      high = "black",
      limits = lims,
      name = NULL,
      guide = guide_colorbar(
        barheight = unit(4.2, "cm"),
        barwidth = unit(0.22, "cm"),
        ticks = TRUE
      )
    )

  if (!keep_legend) {
    p + theme(legend.position = "none")
  } else {
    p + theme(
      legend.position = "right",
      legend.text = element_text(size = 10, face = "bold")
    )
  }
}

method_order <- c("True", "BSGL", "GGP-GAM", "GSVC", "MGWR")
panels <- lapply(method_order, function(m) {
  make_panel(
    plot_data[[m]],
    m,
    panel_lims[[m]],
    keep_legend = (m == "MGWR")
  )
})
names(panels) <- method_order

legend_plot <- panels[["MGWR"]]
lgd <- cowplot::get_legend(legend_plot)
panels[["MGWR"]] <- panels[["MGWR"]] + theme(legend.position = "none")

main_row <- cowplot::plot_grid(
  plotlist = panels,
  nrow = 1,
  rel_widths = rep(1, length(panels)),
  align = "h",
  axis = "tb"
)

final_plot <- cowplot::plot_grid(
  main_row,
  lgd,
  nrow = 1,
  rel_widths = c(1, 0.04)
)

out_file <- file.path(out_dir, "beta5_ownscale_onebar_n1000_p5.png")
ggsave(out_file, final_plot, width = 13, height = 2.6, dpi = 300, bg = "white")

range_df <- data.frame(
  Method = method_order,
  min = vapply(method_order, function(m) panel_lims[[m]][1], numeric(1)),
  max = vapply(method_order, function(m) panel_lims[[m]][2], numeric(1))
)
write.csv(range_df, file.path(out_dir, "beta5_ownscale_ranges.csv"), row.names = FALSE)
print(range_df)
cat("Saved:", out_file, "\n")
