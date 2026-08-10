# Regenerate Main Figure 2: simulation spatial SCP maps.
#
# The figure in the paper corresponds to n = 5000, m = 10, rep_id = 2.
# This script refits BSGL for that replicate, computes pointwise 95% credible
# interval exclusion maps, and writes the figure to figures/main/.

source(file.path("R", "beta_functions.R"))
source(file.path("R", "helpers_gplasso.R"))
source(file.path("R", "help_inclu.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

get_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "") return(default)
  as.integer(value)
}

get_env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "") return(default)
  as.numeric(value)
}

n_val <- 5000
m_val <- 10
rep_id <- 2

fit_L <- get_env_int("FIG2_L", 25)
a_lam <- get_env_num("FIG2_A_LAMBDA", 30)
b_lam <- get_env_num("FIG2_B_LAMBDA", 1)
niter <- get_env_int("FIG2_NITER", 5000)
nburn <- get_env_int("FIG2_NBURN", 500)
seed <- get_env_int("FIG2_SEED", 123)

out_dir <- file.path("figures", "main")
res_dir <- file.path("results", "simulation", "main_fig2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

fit_path <- file.path(res_dir, "main_fig2_n5000_m10_rep02_bsgl_fit.rds")

if (file.exists(fit_path)) {
  fit <- readRDS(fit_path)
  grid_data <- fit$grid_data
} else {
  load(file.path("data", "data_rep10", "n5000_p10.RData"))
  dat <- data_list[[rep_id]]

  set.seed(seed)
  fit0 <- fit_bsgl(
    y = dat$train$y,
    X = dat$train$X,
    coords = dat$train$coords,
    L = fit_L,
    niter = niter,
    nburn = nburn,
    a_lam = a_lam,
    b_lam = b_lam,
    verbose = TRUE
  )

  fit <- list(
    bsgl_fit = fit0,
    grid_data = dat$grid,
    settings = list(
      n = n_val, m = m_val, rep_id = rep_id, seed = seed,
      L = fit_L, a_lam = a_lam, b_lam = b_lam,
      niter = niter, nburn = nburn
    )
  )
  saveRDS(fit, fit_path)
  grid_data <- dat$grid
}

surface_ci <- function(fit_obj, grid_data, ci_level = 0.95) {
  coords <- as.matrix(grid_data$grid_points)
  Psi_new <- make_basis(coords, fit_obj$L, saved_knots = fit_obj$saved_knots)
  n_grid <- nrow(Psi_new)
  n_samples <- dim(fit_obj$beta)[3]
  alpha <- 1 - ci_level

  excludes0 <- matrix(FALSE, n_grid, fit_obj$m)
  for (j in seq_len(fit_obj$m)) {
    beta_samples <- matrix(0, n_grid, n_samples)
    for (s in seq_len(n_samples)) {
      beta_samples[, s] <- Psi_new %*% fit_obj$beta[, j, s]
    }
    lo <- apply(beta_samples, 1, quantile, alpha / 2)
    hi <- apply(beta_samples, 1, quantile, 1 - alpha / 2)
    excludes0[, j] <- (lo > 0) | (hi < 0)
  }
  excludes0
}

excludes <- surface_ci(fit$bsgl_fit, grid_data, ci_level = 0.95)
scp <- colMeans(excludes)

plot_vars <- c(1, 2, 3, 5, 6, 7)
grid_df <- as.data.frame(grid_data$grid_points)
colnames(grid_df) <- c("x", "y")

plot_df <- do.call(rbind, lapply(plot_vars, function(j) {
  data.frame(
    x = grid_df$x,
    y = grid_df$y,
    beta = paste0("\u03b2", j, " (SCP = ", sprintf("%.1f", 100 * scp[j]), "%)"),
    excludes0 = factor(as.numeric(excludes[, j]), levels = c(0, 1))
  )
}))

write.csv(
  data.frame(beta = seq_len(ncol(excludes)), SCP = scp),
  file.path(res_dir, "main_fig2_n5000_m10_rep02_scp.csv"),
  row.names = FALSE
)

p <- ggplot(plot_df, aes(x, y)) +
  geom_point(
    aes(color = excludes0, shape = excludes0),
    size = 0.7,
    stroke = 0.15,
    alpha = 0.95
  ) +
  facet_wrap(~ beta, nrow = 2) +
  coord_fixed() +
  scale_color_manual(values = c("0" = "grey85", "1" = "black"), guide = "none") +
  scale_shape_manual(
    values = c("0" = 1, "1" = 16),
    labels = c("No", "Yes"),
    name = "95% CI excludes 0"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  file.path(out_dir, "fig2_simulation_scp_map.png"),
  p,
  width = 7.2,
  height = 5.0,
  dpi = 300,
  bg = "white"
)

print(scp[plot_vars])
cat("Saved:", normalizePath(file.path(out_dir, "fig2_simulation_scp_map.png")), "\n")
