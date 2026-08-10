source(file.path("R", "helpers_gplasso.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(viridis)
})

fit_path <- file.path(
  "results/interaction/main_effects_rerun",
  "main_interaction_pilot_fit.rds"
)
out_dir <- file.path("results/interaction/main_effects_rerun", "x1_contribution_maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(fit_path)
dat <- obj$data
grid_data <- dat$grid
grid <- as.data.frame(grid_data$grid_points)
colnames(grid) <- c("x", "y")

var_names <- colnames(dat$train$X)
stopifnot(
  var_names[1] == "X1",
  var_names[2] == "X2",
  var_names[6] == "centered_scaled_X1_X2"
)

ci_level <- 0.95
alpha <- 1 - ci_level
lower_prob <- alpha / 2
upper_prob <- 1 - alpha / 2

fit <- obj$bsgl_fit
grid_coords <- as.matrix(grid_data$grid_points)
Psi_grid <- make_basis(grid_coords, fit$L, saved_knots = fit$saved_knots)
n_grid <- nrow(Psi_grid)
n_samples <- dim(fit$beta)[3]

beta_samples_for <- function(j) {
  samples <- matrix(0, nrow = n_grid, ncol = n_samples)
  for (s in seq_len(n_samples)) {
    samples[, s] <- Psi_grid %*% fit$beta[, j, s]
  }
  samples
}

beta1_samples <- beta_samples_for(1)
beta12_samples <- beta_samples_for(6)

summarize_samples <- function(samples) {
  data.frame(
    mean = rowMeans(samples),
    lower = apply(samples, 1, quantile, lower_prob),
    upper = apply(samples, 1, quantile, upper_prob)
  )
}

beta12_df <- cbind(grid, summarize_samples(beta12_samples))
beta12_df$excludes0 <- (beta12_df$lower > 0) | (beta12_df$upper < 0)
beta12_scp <- mean(beta12_df$excludes0)

p_beta12_scp <- ggplot(beta12_df, aes(x, y)) +
  geom_point(
    aes(color = factor(as.numeric(excludes0))),
    alpha = 0.95,
    size = 2.0,
    shape = 15
  ) +
  coord_fixed(expand = FALSE) +
  scale_color_manual(
    values = c("0" = "#DEEBF7", "1" = "#2171B5"),
    labels = c("CI includes 0", "CI excludes 0"),
    name = "",
    guide = guide_legend(
      keywidth = unit(0.8, "cm"),
      keyheight = unit(0.8, "cm"),
      override.aes = list(size = 5, shape = 15)
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "#E8F4F8"),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  ) +
  ggtitle(sprintf("\u03b212 SCP map (SCP = %.3f)", beta12_scp))

ggsave(
  file.path(out_dir, "beta12_bsgl_scp_map.png"),
  p_beta12_scp,
  width = 5.2,
  height = 5.0,
  dpi = 300,
  bg = "white"
)

p_beta12_mean <- ggplot(beta12_df, aes(x, y, fill = mean)) +
  geom_raster() +
  coord_fixed(expand = FALSE) +
  scale_fill_viridis_c(option = "viridis", name = NULL) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  ) +
  ggtitle("Estimated \u03b212")

ggsave(
  file.path(out_dir, "beta12_bsgl_mean_map.png"),
  p_beta12_mean,
  width = 5.2,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

# The fitted design uses centered/scaled X1:X2:
# X12_scaled = scale01((X1 - c1) * (X2 - c2)).
# To show how the X1-related contribution changes with X2, fix X1 and compare
# low vs high X2 values using the exact same coding as the training design.
train_X <- dat$train$X
center_x1 <- mean(train_X[, "X1"])
center_x2 <- mean(train_X[, "X2"])
prod12_train <- (train_X[, "X1"] - center_x1) * (train_X[, "X2"] - center_x2)
prod_min <- min(prod12_train)
prod_max <- max(prod12_train)

scale_prod12 <- function(x1, x2) {
  ((x1 - center_x1) * (x2 - center_x2) - prod_min) / (prod_max - prod_min)
}

x1_fixed <- as.numeric(quantile(train_X[, "X1"], 0.99))
x2_low <- as.numeric(quantile(train_X[, "X2"], 0.01))
x2_high <- as.numeric(quantile(train_X[, "X2"], 0.99))

x12_low <- scale_prod12(x1_fixed, x2_low)
x12_high <- scale_prod12(x1_fixed, x2_high)

contrib_low_samples <- x1_fixed * beta1_samples + x12_low * beta12_samples
contrib_high_samples <- x1_fixed * beta1_samples + x12_high * beta12_samples

contrib_df <- rbind(
  data.frame(
    grid,
    condition = sprintf("low X2 (%.2f)", x2_low),
    summarize_samples(contrib_low_samples)
  ),
  data.frame(
    grid,
    condition = sprintf("high X2 (%.2f)", x2_high),
    summarize_samples(contrib_high_samples)
  )
)
contrib_df$condition <- factor(
  contrib_df$condition,
  levels = unique(contrib_df$condition)
)

contrib_lims <- range(contrib_df$mean, na.rm = TRUE)

p_contrib <- ggplot(contrib_df, aes(x, y, fill = mean)) +
  geom_raster() +
  coord_fixed(expand = FALSE) +
  facet_wrap(~ condition, nrow = 1) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = contrib_lims,
    name = NULL
  ) +
  theme_void(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  ) +
  ggtitle(sprintf("X1 contribution, fixed X1 = %.2f", x1_fixed))

ggsave(
  file.path(out_dir, "x1_contribution_low_high_x2.png"),
  p_contrib,
  width = 8.2,
  height = 4.2,
  dpi = 300,
  bg = "white"
)

contrib_delta_df <- data.frame(
  grid,
  delta = rowMeans(contrib_high_samples) - rowMeans(contrib_low_samples)
)

p_delta <- ggplot(contrib_delta_df, aes(x, y, fill = delta)) +
  geom_raster() +
  coord_fixed(expand = FALSE) +
  scale_fill_viridis_c(option = "viridis", name = NULL) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  ) +
  ggtitle("X1 contribution: high X2 - low X2")

ggsave(
  file.path(out_dir, "x1_contribution_high_minus_low_x2.png"),
  p_delta,
  width = 6.2,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

write.csv(beta12_df, file.path(out_dir, "beta12_bsgl_scp_map.csv"), row.names = FALSE)
write.csv(contrib_df, file.path(out_dir, "x1_contribution_low_high_x2.csv"), row.names = FALSE)
write.csv(contrib_delta_df, file.path(out_dir, "x1_contribution_high_minus_low_x2.csv"), row.names = FALSE)

cat("Saved interaction contribution figures in:", normalizePath(out_dir), "\n")
cat("beta12 SCP:", beta12_scp, "\n")
cat("fixed X1:", x1_fixed, "low X2:", x2_low, "high X2:", x2_high, "\n")
cat("coded X12 low/high:", x12_low, x12_high, "\n")
