source(file.path("R", "beta_functions.R"))
source(file.path("R", "plot_beta.R"))
source(file.path("R", "load_data.R"))
source(file.path("R", "helpers_ggpgam.R"))
source(file.path("R", "helpers_gplasso.R"))
source(file.path("R", "fit_chain.R"))
source(file.path("R", "helpers_gpprior.R"))
source(file.path("R", "cv_ggp.R"))
source(file.path("R", "calc_beta.R"))
source(file.path("R", "beta_mse.R"))
source(file.path("R", "help_inclu.R"))

library(dplyr)
library(ggplot2)
library(fields)
library(mvtnorm)
library(MCMCpack)
library(splines)
library(GIGrvg)

make_abrupt_main_data <- function(n = 1000, p = 5, sigma2 = 0.1, seed = 123) {
  if (!is.null(seed)) set.seed(seed)

  make_x <- function(n_now) {
    x <- rnorm(n_now)
    (x - min(x)) / (max(x) - min(x))
  }

  coords_raw <- cbind(runif(n, 0, 20), runif(n, 0, 20))
  coords <- coords_raw / 20

  X <- matrix(0, nrow = n, ncol = p)
  for (j in seq_len(p)) X[, j] <- make_x(n)

  beta <- calc_cluster_beta(coords_raw, p)
  y <- rowSums(X * beta) + rnorm(n, 0, sqrt(sigma2))

  n_test <- floor(n * 0.2)
  test_coords_raw <- cbind(runif(n_test, 0, 20), runif(n_test, 0, 20))
  test_coords <- test_coords_raw / 20

  test_X <- matrix(0, nrow = n_test, ncol = p)
  for (j in seq_len(p)) test_X[, j] <- make_x(n_test)

  test_beta <- calc_cluster_beta(test_coords_raw, p)
  test_y <- rowSums(test_X * test_beta) + rnorm(n_test, 0, sqrt(sigma2))

  grid_x <- seq(0, 1, length.out = 50)
  grid_y <- seq(0, 1, length.out = 50)
  grid_points <- as.matrix(expand.grid(grid_x, grid_y))
  grid_raw <- grid_points * 20
  true_betas <- calc_cluster_beta(grid_raw, p)

  list(
    train = list(X = X, y = y, coords = coords, true_beta = beta),
    test = list(X = test_X, y = test_y, coords = test_coords, true_beta = test_beta),
    grid = list(grid_points = grid_points, true_betas = true_betas),
    meta = list(n = n, p = p, sigma2 = sigma2, seed = seed)
  )
}

out_dir <- "results/abrupt_boundary/main_simulation_settings_rerun"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data <- make_abrupt_main_data(n = 1000, p = 5, sigma2 = 0.1, seed = 123)
train_data <- data$train
test_data <- data$test
grid_data <- data$grid

n <- nrow(train_data$X)
p <- ncol(train_data$X)
true_beta <- train_data$true_beta
train_grid_data <- list(grid_points = train_data$coords)

cat("Fitting GGP-GAM...\n")
gam_model <- fit_gam(train_data$X, train_data$y, train_data$coords)
pred_gam <- pred_gam(gam_model, test_data$X, test_data$coords)
mspe_gam <- mean((test_data$y - pred_gam)^2)
gam_betas <- get_gam_betas(gam_model, train_data$coords)
mse_1_gam <- mean((true_beta[, 1:3] - gam_betas[, 1:3])^2)
mse_0_gam <- mean((true_beta[, 4:p] - gam_betas[, 4:p])^2)

cat("Cross-validating BSGL...\n")
cv_bsgl <- cv_bsgl(
  y = train_data$y,
  X = train_data$X,
  coords = train_data$coords,
  L_range = c(25, 36),
  a_lambda_range = c(15, 30, 35),
  b_lambda_range = c(0.1, 1, 2),
  niter = 800,
  nburn = 100,
  k_folds = 5
)

best_bsgl <- cv_bsgl$best_params
print(best_bsgl)

cat("Fitting BSGL...\n")
bsgl_fit <- fit_bsgl(
  y = train_data$y,
  X = train_data$X,
  coords = train_data$coords,
  L = best_bsgl$L,
  a_lam = best_bsgl$a_lambda,
  b_lam = best_bsgl$b_lambda,
  niter = 5000,
  nburn = 500,
  verbose = TRUE
)

bsgl_res <- eval_bsgl(bsgl_fit, test_data = test_data)
bsgl_beta <- betas_bsgl(bsgl_fit, train_grid_data)
mse_1_bsgl <- mean((true_beta[, 1:3] - bsgl_beta[, 1:3])^2)
mse_0_bsgl <- mean((true_beta[, 4:p] - bsgl_beta[, 4:p])^2)
mse_avg_bsgl <- mean((true_beta - bsgl_beta)^2)
bsgl_scp <- calc_scp_bsgl(bsgl_fit, grid_data, ci_level = 0.95)

cat("Fitting GSVC...\n")
gs_fit <- fit_gs(
  y = train_data$y,
  X = train_data$X,
  coords = train_data$coords,
  L = best_bsgl$L,
  niter = 5000,
  nburn = 500,
  a_k = 2,
  b_k = 1,
  a_s = 2,
  b_s = 1,
  verbose = TRUE
)

gs_res <- eval_gs(gs_fit, test_data = test_data)
gs_beta <- betas_gs(gs_fit, train_grid_data)
mse_1_gs <- mean((true_beta[, 1:3] - gs_beta[, 1:3])^2)
mse_0_gs <- mean((true_beta[, 4:p] - gs_beta[, 4:p])^2)
mse_avg_gs <- mean((true_beta - gs_beta)^2)
gs_scp <- calc_scp_gs(gs_fit, grid_data, ci_level = 0.95)

summary_df <- data.frame(
  n = n,
  p = p,
  Method = c("GGP-GAM", "GSVC", "BSGL"),
  MSPE = c(mspe_gam, gs_res$mspe, bsgl_res$mspe),
  Coverage = c(NA, gs_res$cover, bsgl_res$cover),
  MSE_1 = c(mse_1_gam, mse_1_gs, mse_1_bsgl),
  MSE_0 = c(mse_0_gam, mse_0_gs, mse_0_bsgl),
  MSE_avg = c(
    mean((true_beta - gam_betas)^2),
    mse_avg_gs,
    mse_avg_bsgl
  )
)

scp_df <- rbind(
  data.frame(Method = "GSVC", Variable = names(gs_scp), SCP = as.numeric(gs_scp)),
  data.frame(Method = "BSGL", Variable = names(bsgl_scp), SCP = as.numeric(bsgl_scp))
)

write.csv(summary_df, file.path(out_dir, "abrupt_boundary_mainsettings_summary.csv"), row.names = FALSE)
write.csv(scp_df, file.path(out_dir, "abrupt_boundary_mainsettings_scp.csv"), row.names = FALSE)
write.csv(cv_bsgl$cv_results, file.path(out_dir, "abrupt_boundary_mainsettings_bsgl_cv.csv"), row.names = FALSE)
saveRDS(
  list(
    bsgl_fit = bsgl_fit,
    gs_fit = gs_fit,
    gam_model = gam_model,
    grid_data = grid_data,
    meta = data$meta,
    summary = summary_df,
    scp = scp_df,
    cv_bsgl = cv_bsgl
  ),
  file.path(out_dir, "abrupt_boundary_mainsettings_fit.rds")
)

print(summary_df)
print(scp_df)
