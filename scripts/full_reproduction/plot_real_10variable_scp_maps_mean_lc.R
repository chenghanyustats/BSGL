source(file.path("R", "helpers_gplasso.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(gridExtra)
  library(grid)
  library(cowplot)
})

out_dir <- "results/real_data/clean_real_data_10variable_rerun"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit_obj <- readRDS("results/real_data/nointeraction_scp_compare_n10000_rerun/real_nointeraction_scp_compare_fit.rds")
fit <- fit_obj$bsgl_fit
var_names <- fit_obj$variables

continuous_vars <- c(
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
lc_vars <- var_names[grepl("^LC_Type4_", var_names)]

display_vars <- c(continuous_vars, "LC_Type4")
display_labels <- c(
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

load(file.path("data", "data_cleaned_small_expanded.RData"))
real_data <- data_cleaned_small

set.seed(20260531)
sample_size <- 10000
idx <- if (nrow(real_data) > sample_size) sample(nrow(real_data), sample_size) else seq_len(nrow(real_data))
sample_data <- real_data[idx, ]

set.seed(456)
test_indices <- sample(seq_len(nrow(sample_data)), size = floor(0.2 * nrow(sample_data)))
train_indices <- setdiff(seq_len(nrow(sample_data)), test_indices)
boundary_indices <- unique(c(
  which(sample_data$scaled_x == min(sample_data$scaled_x)),
  which(sample_data$scaled_x == max(sample_data$scaled_x)),
  which(sample_data$scaled_y == min(sample_data$scaled_y)),
  which(sample_data$scaled_y == max(sample_data$scaled_y))
))
train_indices <- sort(unique(c(train_indices, boundary_indices)))

train_coords <- as.matrix(sample_data[train_indices, c("scaled_x", "scaled_y")])
colnames(train_coords) <- c("coords_x", "coords_y")

grid_n <- 40
grid_coords <- expand.grid(
  x = seq(min(train_coords[, 1]), max(train_coords[, 1]), length.out = grid_n),
  y = seq(min(train_coords[, 2]), max(train_coords[, 2]), length.out = grid_n)
)
grid_points <- as.matrix(grid_coords)

to_lonlat_grid <- function(grid_points, sample_data) {
  grid_x_original <- grid_points[, 1] * diff(range(sample_data$x)) + min(sample_data$x)
  grid_y_original <- grid_points[, 2] * diff(range(sample_data$y)) + min(sample_data$y)
  modis_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
  wgs84_crs <- "+proj=longlat +datum=WGS84 +no_defs"
  coords_modis <- st_as_sf(
    data.frame(x = grid_x_original, y = grid_y_original),
    coords = c("x", "y"),
    crs = modis_crs
  )
  st_coordinates(st_transform(coords_modis, crs = wgs84_crs))
}

read_ne_zip <- function(file_name) {
  zip_path <- normalizePath(file.path("data/geo_naturalearth_10m", file_name))
  shp_name <- sub("\\.zip$", ".shp", file_name)
  st_read(file.path("/vsizip", zip_path, shp_name), quiet = TRUE)
}

make_geo_layers <- function() {
  sf_use_s2(FALSE)
  ocean <- read_ne_zip("ne_10m_ocean.zip")
  lakes <- read_ne_zip("ne_10m_lakes.zip")
  rivers <- read_ne_zip("ne_10m_rivers_lake_centerlines.zip")
  coastline <- read_ne_zip("ne_10m_coastline.zip")
  study_area <- st_as_sfc(st_bbox(
    c(xmin = -125, xmax = -103, ymin = 29, ymax = 41),
    crs = st_crs(4326)
  ))
  layers <- list(
    ocean = st_intersection(ocean, study_area),
    lakes = st_intersection(lakes, study_area),
    rivers = st_intersection(rivers, study_area),
    coast = st_intersection(coastline, study_area)
  )
  sf_use_s2(TRUE)
  layers
}

surface_excludes <- function(fit, grid_points, ci_level = 0.95) {
  alpha <- 1 - ci_level
  lower_prob <- alpha / 2
  upper_prob <- 1 - alpha / 2
  psi_grid <- make_basis(grid_points, fit$L, saved_knots = fit$saved_knots)
  n_grid <- nrow(grid_points)
  n_vars <- fit$m
  n_samples <- dim(fit$beta)[3]
  out <- matrix(FALSE, nrow = n_grid, ncol = n_vars)
  colnames(out) <- var_names

  for (j in seq_len(n_vars)) {
    samples <- matrix(0, nrow = n_grid, ncol = n_samples)
    for (s in seq_len(n_samples)) {
      samples[, s] <- psi_grid %*% fit$beta[, j, s]
    }
    lower <- apply(samples, 1, quantile, lower_prob)
    upper <- apply(samples, 1, quantile, upper_prob)
    out[, j] <- (lower > 0) | (upper < 0)
  }
  out
}

excludes <- surface_excludes(fit, grid_points)
point_scp <- as.data.frame(excludes)
point_scp$LC_Type4 <- rowMeans(point_scp[, lc_vars, drop = FALSE])
scp_rate <- c(
  colMeans(point_scp[, continuous_vars, drop = FALSE]),
  LC_Type4 = mean(point_scp$LC_Type4)
)

lonlat <- to_lonlat_grid(grid_points, sample_data)
layers <- make_geo_layers()

plot_one <- function(var) {
  plot_df <- data.frame(
    lon = lonlat[, 1],
    lat = lonlat[, 2],
    scp_value = as.numeric(point_scp[[var]])
  )
  plot_df <- plot_df[
    plot_df$lon >= -125 & plot_df$lon <= -103 &
      plot_df$lat >= 29 & plot_df$lat <= 41,
  ]

  ggplot() +
    geom_sf(data = layers$ocean, fill = "#C6DBEF", color = NA, alpha = 0.6) +
    geom_point(
      data = plot_df,
      aes(x = lon, y = lat, color = scp_value),
      alpha = 0.95,
      size = 1.2,
      shape = 15
    ) +
    geom_sf(data = layers$lakes, fill = "#4292C6", color = "#2171B5", linewidth = 0.2) +
    geom_sf(data = layers$rivers, color = "#2171B5", linewidth = 0.4, alpha = 0.7) +
    geom_sf(data = layers$coast, color = "gray30", linewidth = 0.6, fill = NA) +
    scale_color_gradient(
      low = "#DEEBF7",
      high = "#2171B5",
      limits = c(0, 1),
      breaks = c(0, 1),
      labels = c("CI includes 0", "CI excludes 0"),
      name = "",
      guide = guide_colorbar(
        barwidth = unit(6.0, "cm"),
        barheight = unit(0.35, "cm"),
        title.position = "top"
      )
    ) +
    coord_sf(xlim = c(-125, -103), ylim = c(29, 41), expand = FALSE) +
    theme_bw(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "#E8F4F8"),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(1, 0, 1, 0)),
      plot.margin = margin(1, 1, 1, 1),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    labs(title = sprintf("%s (SCP = %.3f)", display_labels[var], scp_rate[var]))
}

plots <- lapply(display_vars, plot_one)

legend_plot <- plot_one("NIR_reflectance") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 11, face = "bold"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )
legend <- cowplot::get_legend(legend_plot)

combined <- grid.arrange(
  arrangeGrob(grobs = plots, nrow = 5, ncol = 2),
  legend,
  ncol = 1,
  heights = c(1, 0.055)
)

out_file <- file.path(out_dir, "real_scp_10variable_maps_bsgl_mean_lc.png")
ggsave(out_file, combined, width = 12, height = 20, dpi = 300, bg = "white")

write.csv(
  data.frame(variable = display_vars, label = unname(display_labels[display_vars]), SCP = as.numeric(scp_rate[display_vars])),
  file.path(out_dir, "real_scp_10variable_maps_bsgl_mean_lc.csv"),
  row.names = FALSE
)

cat("Saved:", normalizePath(out_file), "\n")
print(data.frame(variable = display_vars, SCP = as.numeric(scp_rate[display_vars])))
