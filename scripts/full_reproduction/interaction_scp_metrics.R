rm(list = ls())

source(file.path("R", "helpers_gplasso.R"))
source(file.path("R", "helpers_gpprior.R"))

library(ggplot2)

fit_path <- file.path(
  "results/interaction/main_effects_rerun",
  "main_interaction_pilot_fit.rds"
)
out_dir <- "results/interaction/main_effects_rerun"

obj <- readRDS(fit_path)
grid_data <- obj$data$grid
true_beta <- grid_data$true_betas
var_names <- colnames(obj$data$train$X)

grid_coords <- as.matrix(grid_data$grid_points)
grid_df <- as.data.frame(grid_data$grid_points)
colnames(grid_df) <- c("x", "y")

ci_level <- 0.95
alpha <- 1 - ci_level
lower_prob <- alpha / 2
upper_prob <- 1 - alpha / 2

n_grid <- nrow(grid_coords)
calc_metrics_from_ci <- function(method, ci, true_beta, var_names, grid_df) {
  m <- ncol(ci$excludes0)
  metrics <- vector("list", m)
  map_list <- vector("list", m)

  for (j in seq_len(m)) {
    pred_signal <- ci$excludes0[, j]
    true_signal <- abs(true_beta[, j]) > 1e-8

    TP <- sum(pred_signal & true_signal)
    FP <- sum(pred_signal & !true_signal)
    FN <- sum(!pred_signal & true_signal)
    TN <- sum(!pred_signal & !true_signal)

    precision <- ifelse(TP + FP > 0, TP / (TP + FP), NA_real_)
    recall <- ifelse(TP + FN > 0, TP / (TP + FN), NA_real_)
    specificity <- ifelse(TN + FP > 0, TN / (TN + FP), NA_real_)
    fpr <- ifelse(FP + TN > 0, FP / (FP + TN), NA_real_)
    f1 <- ifelse(
      is.na(precision) || is.na(recall) || precision + recall == 0,
      NA_real_,
      2 * precision * recall / (precision + recall)
    )

    metrics[[j]] <- data.frame(
      method = method,
      variable = var_names[j],
      true_status = ifelse(any(true_signal), "Active", "Null"),
      SCP = mean(pred_signal),
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      precision = precision,
      recall = recall,
      specificity = specificity,
      FPR = fpr,
      F1 = f1,
      n_true_signal = sum(true_signal),
      n_pred_signal = sum(pred_signal),
      n_grid = n_grid
    )

    map_list[[j]] <- data.frame(
      grid_df,
      method = method,
      variable = var_names[j],
      true_signal = true_signal,
      pred_signal = pred_signal,
      mean_beta = ci$mean[, j],
      lower = ci$lower[, j],
      upper = ci$upper[, j]
    )
  }

  list(
    metrics = do.call(rbind, metrics),
    maps = do.call(rbind, map_list)
  )
}

bsgl_ci <- {
  fit <- obj$bsgl_fit
  Psi_grid <- make_basis(
    grid_coords,
    fit$L,
    saved_knots = fit$saved_knots
  )
  n_samples <- dim(fit$beta)[3]
  m <- fit$m

  mean_beta <- lower <- upper <- matrix(0, n_grid, m)
  excludes0 <- matrix(FALSE, n_grid, m)

  for (j in seq_len(m)) {
    beta_grid_samples <- matrix(0, nrow = n_grid, ncol = n_samples)
    for (s in seq_len(n_samples)) {
      beta_grid_samples[, s] <- Psi_grid %*% fit$beta[, j, s]
    }

    lower[, j] <- apply(beta_grid_samples, 1, quantile, lower_prob)
    upper[, j] <- apply(beta_grid_samples, 1, quantile, upper_prob)
    mean_beta[, j] <- rowMeans(beta_grid_samples)
    excludes0[, j] <- (lower[, j] > 0) | (upper[, j] < 0)
  }

  list(mean = mean_beta, lower = lower, upper = upper, excludes0 = excludes0)
}

gs_ci <- gs_surface_ci(obj$gs_fit, grid_data, ci_level = ci_level)

bsgl_res <- calc_metrics_from_ci("BSGL", bsgl_ci, true_beta, var_names, grid_df)
gs_res <- calc_metrics_from_ci("Gaussian SVC", gs_ci, true_beta, var_names, grid_df)

metrics_df <- rbind(bsgl_res$metrics, gs_res$metrics)
maps_df <- rbind(bsgl_res$maps, gs_res$maps)
metrics_df$variable <- factor(metrics_df$variable, levels = var_names)
maps_df$variable <- factor(maps_df$variable, levels = var_names)
maps_df$pred_label <- ifelse(maps_df$pred_signal, "CI excludes 0", "CI includes 0")

write.csv(
  subset(metrics_df, method == "BSGL", select = -method),
  file.path(out_dir, "interaction_bsgl_scp_f1_fpr.csv"),
  row.names = FALSE
)

write.csv(
  subset(metrics_df, method == "Gaussian SVC", select = -method),
  file.path(out_dir, "interaction_gs_scp_f1_fpr.csv"),
  row.names = FALSE
)

write.csv(
  metrics_df,
  file.path(out_dir, "interaction_bsgl_gs_scp_f1_fpr.csv"),
  row.names = FALSE
)

write.csv(
  subset(maps_df, method == "BSGL", select = -method),
  file.path(out_dir, "interaction_bsgl_scp_maps_long.csv"),
  row.names = FALSE
)

write.csv(
  subset(maps_df, method == "Gaussian SVC", select = -method),
  file.path(out_dir, "interaction_gs_scp_maps_long.csv"),
  row.names = FALSE
)

scp_map <- ggplot(maps_df, aes(x, y, fill = pred_signal)) +
  geom_raster() +
  coord_fixed(expand = FALSE) +
  facet_wrap(~ variable, ncol = 4) +
  scale_fill_manual(
    values = c(`FALSE` = "grey90", `TRUE` = "#2C7FB8"),
    labels = c(`FALSE` = "CI includes 0", `TRUE` = "CI excludes 0"),
    name = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 11),
    legend.position = "bottom"
  ) +
  labs(
    title = "BSGL SCP Maps: 95% CI Excludes Zero",
    x = NULL,
    y = NULL
  )

ggsave(
  file.path(out_dir, "interaction_bsgl_scp_maps.png"),
  scp_map,
  width = 9.5,
  height = 5.8,
  dpi = 240,
  bg = "white"
)

scp_bar <- ggplot(metrics_df, aes(variable, SCP, fill = true_status)) +
  geom_col(width = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c(Active = "#2C7FB8", Null = "grey70")) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.margin = margin(10, 20, 40, 20)
  ) +
  labs(
    title = "BSGL Spatial Coverage Probability",
    x = NULL,
    y = "SCP",
    fill = "True status"
  )

ggsave(
  file.path(out_dir, "interaction_bsgl_scp_bar.png"),
  scp_bar,
  width = 9,
  height = 4.8,
  dpi = 240,
  bg = "white"
)

print(metrics_df)
cat("Saved metrics and SCP figures in:", normalizePath(out_dir), "\n")
