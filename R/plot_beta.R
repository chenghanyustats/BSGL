# plot_beta <- function(models_to_plot, grid_data, meta_info, j, save_path = NULL, grayscale = FALSE) {
#   
#   # 加载必要的包
#   library(ggplot2)
#   library(viridis)
#   library(cowplot)
#   library(dplyr)
#   
#   # 提取数据
#   grid_points <- grid_data$grid_points
#   
#   # 准备模型名称映射
#   model_names_map <- list("Lasso" = "BSGL", "GAM" = "GGP-GAM", "MGWR" = "MGWR")
#   
#   # 创建绘图数据列表
#   plot_data_list <- list()
#   
#   # 添加True beta作为第一个
#   true_beta <- calc_true_beta(grid_points, meta_info$p)
#   df_true <- data.frame(
#     x = grid_points$x,
#     y = grid_points$y,
#     beta = true_beta[, j]
#   )
#   plot_data_list[["True"]] <- df_true
#   
#   # 遍历模型并获取 beta 值
#   for (model_name in names(models_to_plot)) {
#     model_data <- models_to_plot[[model_name]]
#     if (model_name == "Lasso") {
#       pred_beta <- get_lasso_betas(model_data, grid_data)
#     } else if (model_name == "GAM") {
#       pred_beta <- get_gam_betas(model_data, grid_points)
#     } else if (model_name == "MGWR") {
#       pred_beta <- get_mgwr_betas(model_data, grid_points)
#     } else {
#       next
#     }
#     
#     # 准备该模型的数据框
#     df <- data.frame(
#       x = grid_points$x,
#       y = grid_points$y,
#       beta = pred_beta[, j]
#     )
#     plot_data_list[[model_names_map[[model_name]]]] <- df
#   }
#   
#   # 文件名和路径
#   dir.create("plots", showWarnings = FALSE)
#   num_methods <- length(models_to_plot) + 1  # +1 for True
#   filename_base <- paste0("beta", j, "_", num_methods, "_comparison_n", 
#                           meta_info$n, "_p", meta_info$p, 
#                           if(grayscale) "_gray" else "", ".png")
#   if (is.null(save_path)) {
#     save_path <- file.path("plots", filename_base)
#   }
#   
#   # 创建图（顺序：True, BSGL, GGP-GAM, MGWR）
#   plot_list <- list()
#   method_order <- c("True", "BSGL", "GGP-GAM", "MGWR")
#   
#   for (i in 1:length(method_order)) {
#     method <- method_order[i]
#     if (!method %in% names(plot_data_list)) next
#     
#     df <- plot_data_list[[method]]
#     
#     # 每个图用自己的range绘制
#     p <- ggplot(df, aes(x = x, y = y, fill = beta)) +
#       geom_raster() +
#       coord_fixed() +
#       xlab("") + ylab("") +
#       ggtitle(method) +
#       theme_bw() + 
#       theme(
#         legend.position = "right",
#         legend.key.width = unit(0.3, "cm"),
#         legend.key.height = unit(1.2, "cm"),
#         legend.text = element_text(size = 12, face = "bold"),
#         legend.margin = margin(l = 1, r = 1),
#         legend.box.margin = margin(l = 1, r = 1),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(), 
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
#         plot.margin = margin(t = 1, r = 1, b = 1, l = 1)
#       )
#     
#     # 根据参数选择颜色方案
#     if (grayscale) {
#       p <- p + scale_fill_gradient(
#         low = "white", 
#         high = "black",
#         name = NULL,
#         labels = function(x) sprintf("%.1f", x)
#       )
#     } else {
#       p <- p + scale_fill_viridis_c(
#         name = NULL,
#         labels = function(x) sprintf("%.1f", x)
#       )
#     }
#     
#     plot_list[[method]] <- p
#   }
#   
#   # 组合图形
#   if (length(plot_list) > 0) {
#     combined_plot <- plot_grid(
#       plotlist = plot_list,
#       nrow = 1, ncol = length(plot_list),
#       align = "hv",
#       rel_widths = rep(1, length(plot_list))
#     )
#     
#     # 保存图片
#     ggsave(save_path, combined_plot, 
#            width = 16, height = 4, dpi = 300, bg = "white")
#     
#     cat("Saved:", save_path, "\n")
#   }
# }


plot_beta <- function(models_to_plot, grid_data, meta_info, j,
                      save_path = NULL, grayscale = FALSE, show_title = TRUE) {
  
  # packages
  library(ggplot2)
  library(viridis)
  library(cowplot)
  library(dplyr)
  
  grid_points <- grid_data$grid_points
  
  model_names_map <- list("Lasso" = "BSGL", "GAM" = "GGP-GAM", "MGWR" = "MGWR")
  
  plot_data_list <- list()
  
  # True beta
  true_beta <- calc_true_beta(grid_points, meta_info$p)
  plot_data_list[["True"]] <- data.frame(
    x = grid_points$x,
    y = grid_points$y,
    beta = true_beta[, j]
  )
  
  # models
  for (model_name in names(models_to_plot)) {
    model_data <- models_to_plot[[model_name]]
    
    if (model_name == "Lasso") {
      pred_beta <- get_lasso_betas(model_data, grid_data)
    } else if (model_name == "GAM") {
      pred_beta <- get_gam_betas(model_data, grid_points)
    } else if (model_name == "MGWR") {
      pred_beta <- get_mgwr_betas(model_data, grid_points)
    } else {
      next
    }
    
    plot_data_list[[model_names_map[[model_name]]]] <- data.frame(
      x = grid_points$x,
      y = grid_points$y,
      beta = pred_beta[, j]
    )
  }
  
  # file path
  dir.create("plots", showWarnings = FALSE)
  num_methods <- length(models_to_plot) + 1
  filename_base <- paste0(
    "beta", j, "_", num_methods, "_comparison_n",
    meta_info$n, "_p", meta_info$p,
    if (grayscale) "_gray" else "",
    if (show_title) "" else "_notitle",
    ".png"
  )
  if (is.null(save_path)) save_path <- file.path("plots", filename_base)
  
  # plots: True, BSGL, GGP-GAM, MGWR
  plot_list <- list()
  method_order <- c("True", "BSGL", "GGP-GAM", "MGWR")
  
  for (method in method_order) {
    if (!method %in% names(plot_data_list)) next
    df <- plot_data_list[[method]]
    
    p <- ggplot(df, aes(x = x, y = y, fill = beta)) +
      geom_raster() +
      coord_fixed() +
      xlab("") + ylab("") +
      theme_bw() +
      theme(
        legend.position = "right",              # 每个图都显示 legend
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(1.2, "cm"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.margin = margin(l = 1, r = 1),
        legend.box.margin = margin(l = 1, r = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1)
      )
    
    # show/hide panel title
    if (show_title) {
      p <- p + ggtitle(method)
    }
    
    # per-panel scale (NO limits=...)
    if (grayscale) {
      p <- p + scale_fill_gradient(
        low = "white",
        high = "black",
        name = NULL,
        labels = function(x) sprintf("%.1f", x)
      )
    } else {
      p <- p + scale_fill_viridis_c(
        name = NULL,
        labels = function(x) sprintf("%.1f", x)
      )
    }
    
    plot_list[[method]] <- p
  }
  
  # combine
  if (length(plot_list) > 0) {
    combined_plot <- plot_grid(
      plotlist = plot_list,
      nrow = 1,
      ncol = length(plot_list),
      align = "hv",
      rel_widths = rep(1, length(plot_list))
    )
    
    ggsave(save_path, combined_plot, width = 16, height = 4, dpi = 300, bg = "white")
    cat("Saved:", save_path, "\n")
  }
}



# if same scale
# plot_beta <- function(models_to_plot, grid_data, meta_info, j,
#                       save_path = NULL, grayscale = FALSE,
#                       show_title = TRUE) {   ### NEW
#   
#   # 加载必要的包
#   library(ggplot2)
#   library(viridis)
#   library(cowplot)
#   library(dplyr)
#   
#   # 提取数据
#   grid_points <- grid_data$grid_points
#   
#   # 准备模型名称映射
#   model_names_map <- list("Lasso" = "BSGL", "GAM" = "GGP-GAM", "MGWR" = "MGWR")
#   
#   # 创建绘图数据列表
#   plot_data_list <- list()
#   
#   # 添加True beta作为第一个
#   true_beta <- calc_true_beta(grid_points, meta_info$p)
#   df_true <- data.frame(
#     x = grid_points$x,
#     y = grid_points$y,
#     beta = true_beta[, j]
#   )
#   plot_data_list[["True"]] <- df_true
#   
#   # 遍历模型并获取 beta 值
#   for (model_name in names(models_to_plot)) {
#     model_data <- models_to_plot[[model_name]]
#     if (model_name == "Lasso") {
#       pred_beta <- get_lasso_betas(model_data, grid_data)
#     } else if (model_name == "GAM") {
#       pred_beta <- get_gam_betas(model_data, grid_points)
#     } else if (model_name == "MGWR") {
#       pred_beta <- get_mgwr_betas(model_data, grid_points)
#     } else {
#       next
#     }
#     
#     df <- data.frame(
#       x = grid_points$x,
#       y = grid_points$y,
#       beta = pred_beta[, j]
#     )
#     plot_data_list[[model_names_map[[model_name]]]] <- df
#   }
#   
#   # 统一范围
#   all_beta_values <- unlist(lapply(plot_data_list, function(df) df$beta))
#   global_limits <- range(all_beta_values, na.rm = TRUE)
#   
#   # 文件名和路径
#   dir.create("plots", showWarnings = FALSE)
#   num_methods <- length(models_to_plot) + 1
#   filename_base <- paste0(
#     "beta", j, "_", num_methods,
#     "_comparison_n", meta_info$n,
#     "_p", meta_info$p,
#     if (grayscale) "_gray" else "",
#     ".png"
#   )
#   if (is.null(save_path)) {
#     save_path <- file.path("plots", filename_base)
#   }
#   
#   # 创建图
#   plot_list <- list()
#   method_order <- c("True", "BSGL", "GGP-GAM", "MGWR")
#   
#   for (i in seq_along(method_order)) {
#     method <- method_order[i]
#     if (!method %in% names(plot_data_list)) next
#     
#     df <- plot_data_list[[method]]
#     
#     p <- ggplot(df, aes(x = x, y = y, fill = beta)) +
#       geom_raster() +
#       coord_fixed() +
#       xlab("") + ylab("") +
#       theme_bw() +
#       theme(
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(), 
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
#         plot.margin = margin(t = 1, r = 1, b = 1, l = 1)
#       )
#     
#     ### NEW：是否显示子图标题
#     if (show_title) {
#       p <- p + ggtitle(method)
#     }
#     
#     # 颜色方案 + 统一范围
#     if (grayscale) {
#       p <- p + scale_fill_gradient(
#         low = "white",
#         high = "black",
#         limits = global_limits,
#         name = NULL,
#         labels = function(x) sprintf("%.1f", x)
#       )
#     } else {
#       p <- p + scale_fill_viridis_c(
#         limits = global_limits,
#         name = NULL,
#         labels = function(x) sprintf("%.1f", x)
#       )
#     }
#     
#     # 只有最后一个图显示图例
#     if (i < length(method_order)) {
#       p <- p + theme(legend.position = "none")
#     } else {
#       p <- p + theme(
#         legend.position = "right",
#         legend.key.width = unit(0.3, "cm"),
#         legend.key.height = unit(1.2, "cm"),
#         legend.text = element_text(size = 12, face = "bold"),
#         legend.margin = margin(l = 1, r = 1),
#         legend.box.margin = margin(l = 1, r = 1)
#       )
#     }
#     
#     plot_list[[method]] <- p
#   }
#   
#   # 组合图形
#   if (length(plot_list) > 0) {
#     combined_plot <- plot_grid(
#       plotlist = plot_list,
#       nrow = 1,
#       align = "hv"
#     )
#     
#     ggsave(
#       save_path, combined_plot,
#       width = 16, height = 4,
#       dpi = 300, bg = "white"
#     )
#     
#     cat("Saved:", save_path, "\n")
#   }
# }
