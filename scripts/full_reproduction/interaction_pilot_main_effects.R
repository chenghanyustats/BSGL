rm(list = ls())

source(file.path("R", "beta_functions.R"))
source(file.path("R", "helpers_gplasso.R"))
source(file.path("R", "helpers_gpprior.R"))
source(file.path("R", "helpers_gam.R"))
source(file.path("R", "helpers_mgwr.R"))
source(file.path("R", "help_inclu.R"))

library(ggplot2)

out_dir <- "results/interaction/main_effects_rerun"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

beta_int <- function(x, y) {
  12 * exp(-((x - 10)^2 + (y - 10)^2) / 45)
}

scale01 <- function(x) {
  rng <- range(x)
  (x - rng[1]) / (rng[2] - rng[1])
}

make_interaction_data_main <- function(n = 500, p_main = 5, sigma2 = 0.1, seed = 20260531) {
  set.seed(seed)
  n_test <- floor(0.2 * n)
  grid_size <- 45

  make_x_main <- function(n_now) {
    X <- replicate(p_main, scale01(rnorm(n_now)))
    colnames(X) <- paste0("X", seq_len(p_main))
    X
  }

  train_x_main <- make_x_main(n)
  center <- colMeans(train_x_main[, c(1, 2, 5), drop = FALSE])

  make_block <- function(n_now, X_main) {
    coords_20 <- cbind(runif(n_now, 0, 20), runif(n_now, 0, 20))
    X12 <- scale01((X_main[, 1] - center["X1"]) * (X_main[, 2] - center["X2"]))
    X15 <- scale01((X_main[, 1] - center["X1"]) * (X_main[, 5] - center["X5"]))
    X <- cbind(X_main, X12, X15)
    colnames(X) <- c(
      paste0("X", seq_len(p_main)),
      "centered_scaled_X1_X2",
      "centered_scaled_X1_X5"
    )

    beta <- matrix(0, n_now, p_main + 2)
    beta[, 1] <- beta1(coords_20[, 1], coords_20[, 2])
    beta[, 2] <- beta2(coords_20[, 1], coords_20[, 2])
    beta[, 3] <- beta3(coords_20[, 1], coords_20[, 2])
    beta[, p_main + 1] <- beta_int(coords_20[, 1], coords_20[, 2])
    beta[, p_main + 2] <- 0

    y <- rowSums(X * beta) + rnorm(n_now, 0, sqrt(sigma2))
    list(y = y, X = X, coords = coords_20 / 20, true_beta = beta)
  }

  train <- make_block(n, train_x_main)
  test <- make_block(n_test, make_x_main(n_test))

  grid <- expand.grid(
    x = seq(0, 1, length.out = grid_size),
    y = seq(0, 1, length.out = grid_size)
  )

  grid_20 <- cbind(grid$x * 20, grid$y * 20)
  true_betas <- matrix(0, nrow(grid), p_main + 2)
  true_betas[, 1] <- beta1(grid_20[, 1], grid_20[, 2])
  true_betas[, 2] <- beta2(grid_20[, 1], grid_20[, 2])
  true_betas[, 3] <- beta3(grid_20[, 1], grid_20[, 2])
  true_betas[, p_main + 1] <- beta_int(grid_20[, 1], grid_20[, 2])
  true_betas[, p_main + 2] <- 0

  list(
    train = train[c("y", "X", "coords")],
    test = test[c("y", "X", "coords")],
    grid = list(grid_points = grid, true_betas = true_betas),
    meta = list(n = n, p_main = p_main, p = p_main + 2, sigma2 = sigma2, seed = seed)
  )
}

surface_mse <- function(true_beta, est_beta) {
  data.frame(
    variable = colnames(est_beta),
    mse = colMeans((true_beta - est_beta)^2),
    active = colMeans(abs(true_beta) > 1e-10) > 0
  )
}

method_surface_mse <- function(method, true_beta, est_beta) {
  cbind(method = method, surface_mse(true_beta, est_beta))
}

dat <- make_interaction_data_main(n = 1000, p_main = 5, sigma2 = 0.1, seed = 20260531)
var_names <- colnames(dat$train$X)

cat("Main-effects plus interaction pilot data\n")
cat("n:", dat$meta$n, "p:", dat$meta$p, "\n")
cat("Active variables: X1, X2, X3,", var_names[dat$meta$p - 1], "\n")
cat("Null interaction:", var_names[dat$meta$p], "\n")
cat("Train cor(X1, X2):", round(cor(dat$train$X[, 1], dat$train$X[, 2]), 3), "\n")
cat("Train cor(X1, X1:X2):", round(cor(dat$train$X[, 1], dat$train$X[, dat$meta$p - 1]), 3), "\n")
cat("Train cor(X2, X1:X2):", round(cor(dat$train$X[, 2], dat$train$X[, dat$meta$p - 1]), 3), "\n")
cat("Train cor(X1, X1:X5):", round(cor(dat$train$X[, 1], dat$train$X[, dat$meta$p]), 3), "\n")
cat("Train cor(X5, X1:X5):", round(cor(dat$train$X[, 5], dat$train$X[, dat$meta$p]), 3), "\n\n")

fit_L <- 25
niter <- 1200
nburn <- 300

cat("Fitting BSGL...\n")
bsgl_fit <- fit_bsgl(
  y = dat$train$y,
  X = dat$train$X,
  coords = dat$train$coords,
  L = fit_L,
  niter = niter,
  nburn = nburn,
  a_lam = 20,
  b_lam = 1,
  verbose = TRUE
)

bsgl_pred <- pred_bsgl(bsgl_fit, dat$test$X, dat$test$coords)
bsgl_grid_beta <- betas_bsgl(bsgl_fit, dat$grid)
colnames(bsgl_grid_beta) <- var_names

scp <- calc_scp_bsgl(bsgl_fit, dat$grid, ci_level = 0.95)
names(scp) <- var_names

cat("Fitting Gaussian SVC...\n")
gs_fit <- fit_gs(
  y = dat$train$y,
  X = dat$train$X,
  coords = dat$train$coords,
  L = fit_L,
  niter = niter,
  nburn = nburn,
  a_k = 2,
  b_k = 1,
  a_s = 2,
  b_s = 1,
  verbose = TRUE
)

gs_pred <- pred_gs(gs_fit, dat$test$X, dat$test$coords)
gs_grid_beta <- betas_gs(gs_fit, dat$grid)
colnames(gs_grid_beta) <- var_names
gs_scp <- calc_scp_gs(gs_fit, dat$grid, ci_level = 0.95)
names(gs_scp) <- var_names

cat("Fitting GGP-GAM...\n")
gam_fit <- fit_gam(dat$train$X, dat$train$y, dat$train$coords)
gam_pred <- pred_gam(gam_fit, dat$test$X, dat$test$coords)
gam_grid_beta <- get_gam_betas(gam_fit, as.matrix(dat$grid$grid_points))
colnames(gam_grid_beta) <- var_names

cat("Fitting MGWR...\n")
mgwr_fit <- fit_mgwr(dat$train$X, dat$train$y, dat$train$coords)
mgwr_pred <- pred_mgwr(mgwr_fit, dat$test$X, dat$test$coords)
mgwr_grid_beta <- get_mgwr_betas(mgwr_fit, as.matrix(dat$grid$grid_points))
colnames(mgwr_grid_beta) <- var_names

pred_summary <- data.frame(
  method = c("BSGL", "Gaussian SVC", "GGP-GAM", "MGWR"),
  mspe = c(
    mean((dat$test$y - bsgl_pred$mean)^2),
    mean((dat$test$y - gs_pred$mean)^2),
    mean((dat$test$y - gam_pred)^2),
    mean((dat$test$y - mgwr_pred)^2)
  ),
  coverage = c(
    mean(dat$test$y >= bsgl_pred$lower & dat$test$y <= bsgl_pred$upper),
    mean(dat$test$y >= gs_pred$lower & dat$test$y <= gs_pred$upper),
    NA_real_,
    NA_real_
  )
)

selection_summary <- data.frame(
  variable = var_names,
  active = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
  bsgl_mean_norm = sum_bsgl(bsgl_fit)$mean_norm,
  bsgl_scp_95 = as.numeric(scp),
  gs_mean_norm = sum_gs(gs_fit)$mean_norm,
  gs_scp_95 = as.numeric(gs_scp)
)

mse_summary <- rbind(
  method_surface_mse("BSGL", dat$grid$true_betas, bsgl_grid_beta),
  method_surface_mse("Gaussian SVC", dat$grid$true_betas, gs_grid_beta),
  method_surface_mse("GGP-GAM", dat$grid$true_betas, gam_grid_beta),
  method_surface_mse("MGWR", dat$grid$true_betas, mgwr_grid_beta)
)

write.csv(pred_summary, file.path(out_dir, "main_interaction_pred_summary.csv"), row.names = FALSE)
write.csv(selection_summary, file.path(out_dir, "main_interaction_selection_summary.csv"), row.names = FALSE)
write.csv(mse_summary, file.path(out_dir, "main_interaction_surface_mse.csv"), row.names = FALSE)
saveRDS(
  list(
    data = dat,
    bsgl_fit = bsgl_fit,
    gs_fit = gs_fit,
    gam_fit = gam_fit,
    mgwr_fit = mgwr_fit,
    pred_summary = pred_summary,
    selection_summary = selection_summary,
    bsgl_scp = scp,
    gs_scp = gs_scp,
    mse_summary = mse_summary
  ),
  file.path(out_dir, "main_interaction_pilot_fit.rds")
)

cat("\nPrediction summary\n")
print(pred_summary)
cat("\nSelection / SCP summary\n")
print(selection_summary)
cat("\nSurface MSE summary\n")
print(mse_summary)
cat("\nSaved outputs in:", normalizePath(out_dir), "\n")
