source(file.path("R", "helpers_gplasso.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(gridExtra)
  library(grid)
})

in_dir <- "results/interaction/evi_interaction_validation_n10000_lcdummy_nir_viewz_nir_gpp"
fit_file <- file.path(in_dir, "evi_interaction_validation_fit.rds")
out_file <- file.path(in_dir, "interaction_scp_maps_1x2.png")
out_geo_file <- file.path(in_dir, "interaction_scp_maps_1x2_geo.png")

obj <- readRDS(fit_file)
grid_data <- obj$grid_data
surface_ci <- obj$surface_ci
variables <- obj$variables

interaction_vars <- c(
  "NIR_reflectance_x_view_zenith_angle",
  "NIR_reflectance_x_GPP"
)
interaction_labels <- c(
  NIR_reflectance_x_view_zenith_angle = "NIR x View Zenith",
  NIR_reflectance_x_GPP = "NIR x GPP"
)

grid_df <- data.frame(
  sx = grid_data$grid_points[, 1],
  sy = grid_data$grid_points[, 2]
)

load(file.path("data", "data_cleaned_small_expanded.RData"))
real_data <- data_cleaned_small

set.seed(20260531)
sample_size <- 10000
if (nrow(real_data) > sample_size) {
  sample_indices <- sample(nrow(real_data), sample_size)
  sample_data <- real_data[sample_indices, ]
} else {
  sample_data <- real_data
}

x_range_original <- range(sample_data$x)
y_range_original <- range(sample_data$y)

grid_x_original <- grid_df$sx * diff(x_range_original) + x_range_original[1]
grid_y_original <- grid_df$sy * diff(y_range_original) + y_range_original[1]

modis_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
wgs84_crs <- "+proj=longlat +datum=WGS84 +no_defs"

coords_modis <- st_as_sf(
  data.frame(x = grid_x_original, y = grid_y_original),
  coords = c("x", "y"),
  crs = modis_crs
)
coords_latlon <- st_transform(coords_modis, crs = wgs84_crs)
lonlat_grid <- st_coordinates(coords_latlon)

grid_df$lon <- lonlat_grid[, 1]
grid_df$lat <- lonlat_grid[, 2]

read_ne_zip <- function(file_name) {
  zip_path <- normalizePath(file.path("data/geo_naturalearth_10m", file_name))
  shp_name <- sub("\\.zip$", ".shp", file_name)
  st_read(file.path("/vsizip", zip_path, shp_name), quiet = TRUE)
}

sf_use_s2(FALSE)
ocean <- read_ne_zip("ne_10m_ocean.zip")
lakes <- read_ne_zip("ne_10m_lakes.zip")
rivers <- read_ne_zip("ne_10m_rivers_lake_centerlines.zip")
coastline <- read_ne_zip("ne_10m_coastline.zip")

study_area <- st_as_sfc(st_bbox(
  c(xmin = -125, xmax = -103, ymin = 29, ymax = 41),
  crs = st_crs(4326)
))

ocean_crop <- st_intersection(ocean, study_area)
lakes_crop <- st_intersection(lakes, study_area)
rivers_crop <- st_intersection(rivers, study_area)
coast_crop <- st_intersection(coastline, study_area)
sf_use_s2(TRUE)

plot_one_scp <- function(var) {
  j <- match(var, variables)
  if (is.na(j)) stop("Variable not found in fit: ", var)

  excludes0 <- surface_ci$excludes0[, j]
  scp <- mean(excludes0)
  plot_df <- data.frame(
    lon = grid_df$lon,
    lat = grid_df$lat,
    scp_binary = factor(as.numeric(excludes0), levels = c(0, 1))
  )
  plot_df <- plot_df[
    plot_df$lon >= -125 & plot_df$lon <= -103 &
      plot_df$lat >= 29 & plot_df$lat <= 41,
  ]

  ggplot() +
    geom_sf(data = ocean_crop, fill = "#C6DBEF", color = NA, alpha = 0.6) +
    geom_point(
      data = plot_df,
      aes(x = lon, y = lat, color = scp_binary),
      alpha = 0.95,
      size = 1.2,
      shape = 15
    ) +
    geom_sf(data = lakes_crop, fill = "#4292C6", color = "#2171B5", linewidth = 0.2) +
    geom_sf(data = rivers_crop, color = "#2171B5", linewidth = 0.4, alpha = 0.7) +
    geom_sf(data = coast_crop, color = "gray30", linewidth = 0.6, fill = NA) +
    scale_color_manual(
      values = c("0" = "#DEEBF7", "1" = "#2171B5"),
      labels = c("CI includes 0", "CI excludes 0"),
      name = "",
      guide = guide_legend(
        keywidth = unit(0.7, "cm"),
        keyheight = unit(0.7, "cm"),
        override.aes = list(size = 4.5, shape = 15)
      )
    ) +
    coord_sf(xlim = c(-125, -103), ylim = c(29, 41), expand = FALSE) +
    theme_bw(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "#E8F4F8"),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        face = "bold",
        margin = margin(0, 0, 1, 0)
      ),
      plot.margin = margin(0, 1, 0, 1),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(sprintf("%s (SCP = %.3f)", interaction_labels[var], scp))
}

plot_one_scp_scaled <- function(var) {
  j <- match(var, variables)
  if (is.na(j)) stop("Variable not found in fit: ", var)

  excludes0 <- surface_ci$excludes0[, j]
  scp <- mean(excludes0)
  plot_df <- data.frame(
    x = grid_df$sx,
    y = grid_df$sy,
    scp_binary = factor(as.numeric(excludes0), levels = c(0, 1))
  )

  ggplot(plot_df, aes(x, y)) +
    geom_point(
      aes(color = scp_binary),
      alpha = 0.95,
      size = 2.0,
      shape = 15
    ) +
    scale_color_manual(
      values = c("0" = "#DEEBF7", "1" = "#2171B5"),
      labels = c("CI includes 0", "CI excludes 0"),
      name = "",
      guide = guide_legend(
        keywidth = unit(0.7, "cm"),
        keyheight = unit(0.7, "cm"),
        override.aes = list(size = 4.5, shape = 15)
      )
    ) +
    coord_fixed(expand = FALSE) +
    theme_bw(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "#E8F4F8"),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        face = "bold",
        margin = margin(2, 0, 2, 0)
      ),
      plot.margin = margin(2, 2, 2, 2),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(sprintf("%s (SCP = %.3f)", interaction_labels[var], scp))
}

plots <- lapply(interaction_vars, plot_one_scp)

legend_df <- data.frame(
  category = factor(
    c("CI includes 0", "CI excludes 0"),
    levels = c("CI includes 0", "CI excludes 0")
  ),
  x = 1:2,
  y = 1
)

legend_plot <- ggplot(legend_df, aes(x, y, color = category)) +
  geom_point(size = 5, shape = 15) +
  scale_color_manual(
    values = c("CI includes 0" = "#DEEBF7", "CI excludes 0" = "#2171B5"),
    name = ""
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-4, 0, 0, 0)
  )

legend <- cowplot::get_legend(legend_plot)

combined <- grid.arrange(
  arrangeGrob(grobs = plots, ncol = 2),
  legend,
  ncol = 1,
  heights = c(1, 0.075)
)

ggsave(out_file, combined, width = 8.8, height = 3.7, dpi = 300, bg = "white")
ggsave(out_geo_file, combined, width = 8.8, height = 3.7, dpi = 300, bg = "white")

cat("Saved:", normalizePath(out_file), "\n")
cat("Saved:", normalizePath(out_geo_file), "\n")
