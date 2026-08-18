# ===========================================
# 地图绘制函数集合
# ===========================================

library(ggplot2)
library(sf)
library(viridis)
library(gridExtra)
library(ggspatial)

# ===========================================
# 绘制观测EVI地图
# ===========================================

plot_observed_evi <- function(data, ocean_crop, coast_crop, lakes_crop, rivers_crop, 
                              save_path = "figures/real_data/observed_evi.png") {
  
  p <- ggplot() +
    geom_sf(data = ocean_crop, fill = "#C6DBEF", color = NA, alpha = 0.6) +
    geom_sf(data = coast_crop, color = "gray30", linewidth = 0.6) +
    geom_sf(data = lakes_crop, fill = "#4292C6", color = "#2171B5", linewidth = 0.2) +
    geom_sf(data = rivers_crop, color = "#2171B5", linewidth = 0.4, alpha = 0.7) +
    
    geom_point(data = data, 
               aes(x = lon, y = lat, color = log_EVI),
               alpha = 0.8, size = 1.3) +
    
    scale_color_viridis_c(
      name = "",
      option = "viridis",
      breaks = seq(0, 0.5, 0.1),
      labels = sprintf("%.1f", seq(0, 0.5, 0.1)),
      guide = guide_colorbar(
        barwidth = 0.8,
        barheight = 10,
        frame.colour = NA,
        ticks.colour = "black",
        ticks.linewidth = 0.5
      )
    ) +
    
    annotate("text", x = -124, y = 35, label = "Pacific\nOcean", 
             color = "#08519C", size = 4.5, fontface = "italic") +
    
    scale_x_continuous(
      breaks = seq(-124, -104, 4),
      labels = function(x) paste0(abs(x), "°W")
    ) +
    scale_y_continuous(
      breaks = seq(30, 40, 2),
      labels = function(y) paste0(y, "°N")
    ) +
    
    annotation_north_arrow(
      location = "tl", which_north = "true",
      pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
      style = north_arrow_nautical()
    ) +
    annotation_scale(
      location = "bl", width_hint = 0.25,
      pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in")
    ) +
    
    coord_sf(xlim = c(-125, -103), ylim = c(29, 41)) +
    
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "#E8F4F8"),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 9)
    ) +
    labs(title = "Observed EVI", x = "Longitude", y = "Latitude")
  
  ggsave(save_path, p, width = 12, height = 9, dpi = 300, bg = "white")
  return(p)
}

# ===========================================
# 绘制观测vs预测对比图
# ===========================================

plot_observed_vs_predicted <- function(test_observed_plot, test_predicted_plot,
                                       ocean_crop, coast_crop, lakes_crop, rivers_crop,
                                       save_path = "figures/real_data/observed_vs_predicted.png") {
  
  # 观测图
  p_observed <- ggplot() +
    geom_sf(data = ocean_crop, fill = "#C6DBEF", color = NA, alpha = 0.6) +
    geom_sf(data = coast_crop, color = "gray30", linewidth = 0.6) +
    geom_sf(data = lakes_crop, fill = "#4292C6", color = "#2171B5", linewidth = 0.2) +
    geom_sf(data = rivers_crop, color = "#2171B5", linewidth = 0.4, alpha = 0.7) +
    
    geom_point(data = test_observed_plot, 
               aes(x = lon, y = lat, color = observed),
               alpha = 0.8, size = 1.5) +
    
    scale_color_viridis_c(
      name = "",
      option = "viridis",
      breaks = seq(0, 0.5, 0.1),
      labels = sprintf("%.1f", seq(0, 0.5, 0.1)),
      guide = guide_colorbar(
        barwidth = 1, barheight = 12,
        frame.colour = NA,
        ticks.colour = "black",
        ticks.linewidth = 0.5
      )
    ) +
    
    annotate("text", x = -124, y = 35, label = "Pacific\nOcean", 
             color = "#08519C", size = 5, fontface = "bold.italic") +
    
    scale_x_continuous(
      breaks = seq(-124, -104, 4),
      labels = function(x) paste0(abs(x), "°W")
    ) +
    scale_y_continuous(
      breaks = seq(30, 40, 2),
      labels = function(y) paste0(y, "°N")
    ) +
    
    annotation_north_arrow(
      location = "tl", which_north = "true",
      pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
      style = north_arrow_nautical()
    ) +
    annotation_scale(
      location = "bl", width_hint = 0.25,
      pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in")
    ) +
    
    coord_sf(xlim = c(-125, -103), ylim = c(29, 41)) +
    
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "#E8F4F8"),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", 
                                margin = margin(2, 0, 2, 0)),
      plot.margin = margin(2, 5, 2, 5),
      legend.position = "right",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 10, face = "bold"),
      legend.margin = margin(0, 0, 0, 0),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 11, face = "bold")
    ) +
    labs(title = "Observed EVI", x = "Longitude", y = "Latitude")
  
  # 预测图
  p_predicted <- ggplot() +
    geom_sf(data = ocean_crop, fill = "#C6DBEF", color = NA, alpha = 0.6) +
    geom_sf(data = coast_crop, color = "gray30", linewidth = 0.6) +
    geom_sf(data = lakes_crop, fill = "#4292C6", color = "#2171B5", linewidth = 0.2) +
    geom_sf(data = rivers_crop, color = "#2171B5", linewidth = 0.4, alpha = 0.7) +
    
    geom_point(data = test_predicted_plot, 
               aes(x = lon, y = lat, color = predicted),
               alpha = 0.8, size = 1.5) +
    
    scale_color_viridis_c(
      name = "",
      option = "viridis",
      breaks = seq(0, 0.5, 0.1),
      labels = sprintf("%.1f", seq(0, 0.5, 0.1)),
      guide = guide_colorbar(
        barwidth = 1, barheight = 12,
        frame.colour = NA,
        ticks.colour = "black",
        ticks.linewidth = 0.5
      )
    ) +
    
    annotate("text", x = -124, y = 35, label = "Pacific\nOcean", 
             color = "#08519C", size = 5, fontface = "bold.italic") +
    
    scale_x_continuous(
      breaks = seq(-124, -104, 4),
      labels = function(x) paste0(abs(x), "°W")
    ) +
    scale_y_continuous(
      breaks = seq(30, 40, 2),
      labels = function(y) paste0(y, "°N")
    ) +
    
    annotation_north_arrow(
      location = "tl", which_north = "true",
      pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
      style = north_arrow_nautical()
    ) +
    annotation_scale(
      location = "bl", width_hint = 0.25,
      pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in")
    ) +
    
    coord_sf(xlim = c(-125, -103), ylim = c(29, 41)) +
    
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "#E8F4F8"),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", 
                                margin = margin(2, 0, 2, 0)),
      plot.margin = margin(2, 5, 2, 5),
      legend.position = "right",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 10, face = "bold"),
      legend.margin = margin(0, 0, 0, 0),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 11, face = "bold")
    ) +
    labs(title = "Predicted EVI (95% CI)", x = "Longitude", y = "Latitude")
  
  # 组合保存
  combined <- grid.arrange(
    p_observed, p_predicted, 
    ncol = 2,
    padding = unit(0.2, "line")
  )
  
  ggsave(save_path, combined, 
         width = 24, height = 8.5, dpi = 300, bg = "white", limitsize = FALSE)
  
  return(list(observed = p_observed, predicted = p_predicted, combined = combined))
}

# ===========================================
# 绘制系数地图
# ===========================================

plot_coefficient_maps <- function(lasso_gam, train_coords, sample_real_data,
                                  ocean_crop, coast_crop, pip_rate,
                                  n_grid = 80,
                                  save_path = "figures/real_data/coefficient_maps.png") {
  
  sf_use_s2(FALSE)
  
  # 变量名
  variable_names <- c(
    "β1" = "Red Reflectance", "β2" = "NIR Reflectance", 
    "β3" = "Blue Reflectance", "β4" = "MIR Reflectance",
    "β5" = "GPP", "β6" = "LE",
    "β7" = "View Zenith Angle", "β8" = "Sun Zenith Angle",
    "β9" = "Relative Azimuth Angle", "β10" = "LC Type4"
  )
  
  # 选择变量
  top3_indices <- order(pip_rate, decreasing = TRUE)[1:3]
  top3_vars <- names(pip_rate)[top3_indices]
  bottom3_indices <- order(pip_rate, decreasing = FALSE)[1:3]
  bottom3_vars <- names(pip_rate)[bottom3_indices]
  
  # 创建密集grid
  x_range <- range(train_coords[,1])
  y_range <- range(train_coords[,2])
  
  grid_x <- seq(x_range[1], x_range[2], length.out = n_grid)
  grid_y <- seq(y_range[1], y_range[2], length.out = n_grid)
  grid_coords_dense <- expand.grid(x = grid_x, y = grid_y)
  
  grid_data_dense <- list(grid_points = as.matrix(grid_coords_dense))
  
  # 计算beta
  grid_betas_dense <- betas_bsgl(lasso_gam, grid_data_dense)
  
  # 转换坐标
  x_range_original <- range(sample_real_data$x)
  y_range_original <- range(sample_real_data$y)
  
  grid_x_original <- grid_coords_dense[,1] * (x_range_original[2] - x_range_original[1]) + x_range_original[1]
  grid_y_original <- grid_coords_dense[,2] * (y_range_original[2] - y_range_original[1]) + y_range_original[1]
  
  modis_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
  wgs84_crs <- "+proj=longlat +datum=WGS84 +no_defs"
  
  temp_df <- data.frame(x = grid_x_original, y = grid_y_original)
  coords_modis <- st_as_sf(temp_df, coords = c("x", "y"), crs = modis_crs)
  coords_latlon <- st_transform(coords_modis, crs = wgs84_crs)
  lonlat_grid <- st_coordinates(coords_latlon)
  
  # 绘图函数
  plot_single_beta <- function(beta_values, var_name, pip_value, lon, lat) {
    
    plot_data <- data.frame(lon = lon, lat = lat, beta = beta_values)
    plot_data <- plot_data[plot_data$lon >= -125 & plot_data$lon <= -103 &
                             plot_data$lat >= 29 & plot_data$lat <= 41, ]
    
    title_text <- paste0(variable_names[var_name], " (PIP = ", sprintf("%.3f", pip_value), ")")
    
    p <- ggplot() +
      geom_sf(data = ocean_crop, fill = "#C6DBEF", color = NA, alpha = 0.4) +
      
      geom_point(data = plot_data, 
                 aes(x = lon, y = lat, color = beta),
                 alpha = 0.95, size = 1.5, shape = 15) +
      
      geom_sf(data = coast_crop, color = "gray30", linewidth = 0.5, fill = NA) +
      
      scale_color_viridis_c(
        name = "",
        option = "viridis",
        guide = guide_colorbar(
          barwidth = 0.8, barheight = 8,
          frame.colour = NA,
          ticks.colour = "black",
          ticks.linewidth = 0.4
        )
      ) +
      
      coord_sf(xlim = c(-125, -103), ylim = c(29, 41), expand = FALSE) +
      
      theme_bw() +
      theme(
        panel.background = element_rect(fill = "#E8F4F8"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold",
                                  margin = margin(0, 0, 0, 0)),
        plot.margin = margin(0, 5, 0, 5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0)
      ) +
      labs(title = title_text)
    
    return(p)
  }
  
  # 创建6个maps
  selected_vars <- c(top3_vars, bottom3_vars)
  selected_indices <- as.numeric(gsub("β", "", selected_vars))
  
  plot_list <- list()
  
  for (i in 1:6) {
    var_idx <- selected_indices[i]
    var_name <- selected_vars[i]
    
    p <- plot_single_beta(
      beta_values = grid_betas_dense[, var_idx],
      var_name = var_name,
      pip_value = pip_rate[var_name],
      lon = lonlat_grid[, 1],
      lat = lonlat_grid[, 2]
    )
    
    plot_list[[i]] <- p
  }
  
  # 组合保存
  combined_beta <- grid.arrange(
    grobs = plot_list,
    nrow = 2,
    ncol = 3,
    padding = unit(0, "line")
  )
  
  ggsave(save_path, combined_beta,
         width = 20, height = 9.5, dpi = 300, bg = "white")
  
  sf_use_s2(TRUE)
  
  return(combined_beta)
}
