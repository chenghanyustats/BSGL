source(file.path("R", "helpers_gplasso.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(gridExtra)
  library(grid)
  library(cowplot)
})

out_dir <- "figures/reproduction_outputs/real_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

scale01_train_test <- function(train_x, test_x) {
  train_scaled <- train_x
  test_scaled <- test_x
  for (j in seq_len(ncol(train_x))) {
    lo <- min(train_x[, j], na.rm = TRUE)
    hi <- max(train_x[, j], na.rm = TRUE)
    denom <- hi - lo
    if (!is.finite(denom) || denom == 0) {
      train_scaled[, j] <- 0
      test_scaled[, j] <- 0
    } else {
      train_scaled[, j] <- (train_x[, j] - lo) / denom
      test_scaled[, j] <- (test_x[, j] - lo) / denom
    }
  }
  colnames(train_scaled) <- colnames(train_x)
  colnames(test_scaled) <- colnames(test_x)
  list(train = train_scaled, test = test_scaled)
}

make_lc_dummy_design <- function(train_lc, test_lc) {
  train_levels <- sort(unique(train_lc))
  ref_level <- train_levels[1]
  dummy_levels <- train_levels[train_levels != ref_level]

  train_dummy <- sapply(dummy_levels, function(level) as.numeric(train_lc == level))
  test_dummy <- sapply(dummy_levels, function(level) as.numeric(test_lc == level))
  if (length(dummy_levels) == 1) {
    train_dummy <- matrix(train_dummy, ncol = 1)
    test_dummy <- matrix(test_dummy, ncol = 1)
  }

  colnames(train_dummy) <- paste0("LC_Type4_", dummy_levels)
  colnames(test_dummy) <- paste0("LC_Type4_", dummy_levels)
  list(train = train_dummy, test = test_dummy)
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

to_lonlat <- function(sample_data, coords_scaled) {
  x_original <- coords_scaled[, 1] * diff(range(sample_data$x)) + min(sample_data$x)
  y_original <- coords_scaled[, 2] * diff(range(sample_data$y)) + min(sample_data$y)

  modis_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
  wgs84_crs <- "+proj=longlat +datum=WGS84 +no_defs"
  coords_modis <- st_as_sf(
    data.frame(x = x_original, y = y_original),
    coords = c("x", "y"),
    crs = modis_crs
  )
  st_coordinates(st_transform(coords_modis, crs = wgs84_crs))
}

load(file.path("data", "data_cleaned_small_expanded.RData"))
real_data <- data_cleaned_small

set.seed(20260531)
idx <- if (nrow(real_data) > 10000) sample(nrow(real_data), 10000) else seq_len(nrow(real_data))
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
test_indices <- setdiff(test_indices, boundary_indices)
train_indices <- sort(unique(c(train_indices, boundary_indices)))

test_y <- log(sample_data$EVI[test_indices] + 1)

scaled_cont <- scale01_train_test(
  as.matrix(sample_data[train_indices, continuous_vars]),
  as.matrix(sample_data[test_indices, continuous_vars])
)
lc <- make_lc_dummy_design(
  train_lc = sample_data$LC_Type4[train_indices],
  test_lc = sample_data$LC_Type4[test_indices]
)

test_X <- cbind(scaled_cont$test, lc$test)
test_coords <- as.matrix(sample_data[test_indices, c("scaled_x", "scaled_y")])
colnames(test_coords) <- c("coords_x", "coords_y")

fit <- readRDS("results/real_data/nointeraction_scp_compare_n10000_rerun/real_nointeraction_scp_compare_fit.rds")
set.seed(1)
pred <- pred_bsgl(fit$bsgl_fit, test_X, test_coords, return_samples = TRUE)

q05 <- apply(pred$y_samples, 1, quantile, probs = 0.05)
q95 <- apply(pred$y_samples, 1, quantile, probs = 0.95)
latent_q05 <- apply(pred$mu_samples, 1, quantile, probs = 0.05)
latent_q95 <- apply(pred$mu_samples, 1, quantile, probs = 0.95)
lonlat <- to_lonlat(sample_data, test_coords)

plot_data <- data.frame(
  lon = lonlat[, 1],
  lat = lonlat[, 2],
  observed = test_y,
  q05 = q05,
  mean = pred$mean,
  q95 = q95,
  latent_q05 = latent_q05,
  latent_q95 = latent_q95
)

value_range <- range(c(plot_data$q05, plot_data$mean, plot_data$q95), na.rm = TRUE)
layers <- make_geo_layers()

plot_evi_map <- function(value_col, title_text, legend_title = "log(EVI+1)", show_legend = FALSE) {
  ggplot() +
    geom_sf(data = layers$ocean, fill = "#C6DBEF", color = NA, alpha = 0.6) +
    geom_point(
      data = plot_data,
      aes(x = lon, y = lat, color = .data[[value_col]]),
      size = 0.9,
      alpha = 0.78
    ) +
    geom_sf(data = layers$lakes, fill = "#4292C6", color = "#2171B5", linewidth = 0.2) +
    geom_sf(data = layers$rivers, color = "#2171B5", linewidth = 0.4, alpha = 0.7) +
    geom_sf(data = layers$coast, color = "gray30", linewidth = 0.6, fill = NA) +
    scale_color_viridis_c(
      name = legend_title,
      limits = value_range,
      guide = guide_colorbar(
        barheight = unit(2.4, "cm"),
        barwidth = unit(0.28, "cm"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    coord_sf(xlim = c(-125, -103), ylim = c(29, 41), expand = FALSE) +
    theme_bw(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "#E8F4F8"),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(0, 0, 1, 0)),
      plot.margin = margin(0, 1, 0, 1),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(title = title_text)
}

p_q05 <- plot_evi_map("q05", "5% Predictive Quantile")
p_mean <- plot_evi_map("mean", "Posterior Mean")
p_q95 <- plot_evi_map("q95", "95% Predictive Quantile")
shared_legend <- get_legend(plot_evi_map("mean", "Posterior Mean", show_legend = TRUE))

combined_maps <- grid.arrange(
  arrangeGrob(p_q05, p_mean, p_q95, ncol = 3),
  shared_legend,
  ncol = 2,
  widths = c(1, 0.07)
)
ggsave(
  file.path(out_dir, "test_evi_predictive_q05_mean_q95_maps.png"),
  combined_maps,
  width = 13.8,
  height = 3.45,
  dpi = 300,
  bg = "white"
)

density_data <- rbind(
  data.frame(type = "5% predictive quantile", value = plot_data$q05),
  data.frame(type = "Posterior mean", value = plot_data$mean),
  data.frame(type = "95% predictive quantile", value = plot_data$q95),
  data.frame(type = "Observed test EVI", value = plot_data$observed)
)
density_data$type <- factor(
  density_data$type,
  levels = c("Observed test EVI", "5% predictive quantile", "Posterior mean", "95% predictive quantile")
)

p_density <- ggplot(density_data, aes(x = value, color = type, fill = type)) +
  geom_density(alpha = 0.12, linewidth = 1.1) +
  scale_color_manual(values = c(
    "Observed test EVI" = "gray25",
    "5% predictive quantile" = "#2C7BB6",
    "Posterior mean" = "#1A9850",
    "95% predictive quantile" = "#D7191C"
  )) +
  scale_fill_manual(values = c(
    "Observed test EVI" = "gray25",
    "5% predictive quantile" = "#2C7BB6",
    "Posterior mean" = "#1A9850",
    "95% predictive quantile" = "#D7191C"
  )) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Test EVI Predictive Distribution",
    x = "log(EVI+1)",
    y = "Density"
  )

ggsave(
  file.path(out_dir, "test_evi_predictive_q05_mean_q95_density.png"),
  p_density,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

summary_df <- data.frame(
  Quantity = c(
    "Observed test EVI",
    "5% predictive quantile",
    "Posterior mean",
    "95% predictive quantile",
    "5% latent quantile",
    "95% latent quantile",
    "Predictive interval width q95-q05",
    "Latent interval width q95-q05"
  ),
  Mean = c(
    mean(plot_data$observed),
    mean(plot_data$q05),
    mean(plot_data$mean),
    mean(plot_data$q95),
    mean(plot_data$latent_q05),
    mean(plot_data$latent_q95),
    mean(plot_data$q95 - plot_data$q05),
    mean(plot_data$latent_q95 - plot_data$latent_q05)
  ),
  SD = c(
    sd(plot_data$observed),
    sd(plot_data$q05),
    sd(plot_data$mean),
    sd(plot_data$q95),
    sd(plot_data$latent_q05),
    sd(plot_data$latent_q95),
    sd(plot_data$q95 - plot_data$q05),
    sd(plot_data$latent_q95 - plot_data$latent_q05)
  ),
  Min = c(
    min(plot_data$observed),
    min(plot_data$q05),
    min(plot_data$mean),
    min(plot_data$q95),
    min(plot_data$latent_q05),
    min(plot_data$latent_q95),
    min(plot_data$q95 - plot_data$q05),
    min(plot_data$latent_q95 - plot_data$latent_q05)
  ),
  Max = c(
    max(plot_data$observed),
    max(plot_data$q05),
    max(plot_data$mean),
    max(plot_data$q95),
    max(plot_data$latent_q05),
    max(plot_data$latent_q95),
    max(plot_data$q95 - plot_data$q05),
    max(plot_data$latent_q95 - plot_data$latent_q05)
  )
)

write.csv(
  summary_df,
  file.path(out_dir, "test_evi_predictive_q05_mean_q95_summary.csv"),
  row.names = FALSE
)

cat("Saved maps:", normalizePath(file.path(out_dir, "test_evi_predictive_q05_mean_q95_maps.png")), "\n")
cat("Saved density:", normalizePath(file.path(out_dir, "test_evi_predictive_q05_mean_q95_density.png")), "\n")
cat("Saved summary:", normalizePath(file.path(out_dir, "test_evi_predictive_q05_mean_q95_summary.csv")), "\n")
print(summary_df)
