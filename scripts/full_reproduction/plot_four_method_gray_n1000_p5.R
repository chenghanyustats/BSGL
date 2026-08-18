source(file.path("R", "beta_functions.R"))
source(file.path("R", "plot_beta.R"))
source(file.path("R", "load_data.R"))
source(file.path("R", "helpers_ggpgam.R"))
source(file.path("R", "helpers_gplasso.R"))
source(file.path("R", "fit_chain.R"))
source(file.path("R", "helpers_gpprior.R"))
source(file.path("R", "helpers_mgwr.R"))
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

data <- load_data(n = 1000, p = 5, save_dir = "data/data_indep")

train_data <- data$train
test_data <- data$test
grid_data <- data$grid
meta_info <- data$meta

n <- nrow(train_data$X)
p <- ncol(train_data$X)
true_beta <- calc_true_beta(train_data$coords, p)
train_grid_data <- list(grid_points = train_data$coords)

gam_model <- fit_gam(train_data$X, train_data$y, train_data$coords)
pred_gam <- pred_gam(gam_model, test_data$X, test_data$coords)
mspe_gam <- mean((test_data$y - pred_gam)^2)
gam_betas <- get_gam_betas(gam_model, train_data$coords)
mse_1_gam <- mean((true_beta[, 1:3] - gam_betas[, 1:3])^2)
mse_0_gam <- mean((true_beta[, 4:p] - gam_betas[, 4:p])^2)

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
scp_95 <- calc_scp_bsgl(bsgl_fit, grid_data, ci_level = 0.95)

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

mgwr_model <- fit_mgwr(train_data$X, train_data$y, train_data$coords)
pred_m <- pred_mgwr(mgwr_model, test_data$X, test_data$coords)
mspe_mgwr <- mean((test_data$y - pred_m)^2)
mgwr_beta <- get_mgwr_betas(mgwr_model, train_data$coords)
mse_1_mgwr <- mean((true_beta[, 1:3] - mgwr_beta[, 1:3])^2)
mse_0_mgwr <- mean((true_beta[, 4:p] - mgwr_beta[, 4:p])^2)

summary_df <- data.frame(
  n = n,
  p = p,
  Method = c("GGP-GAM", "MGWR", "BSGL", "Gaussian SVC"),
  MSPE = c(mspe_gam, mspe_mgwr, bsgl_res$mspe, gs_res$mspe),
  Coverage = c(NA, NA, bsgl_res$cover, gs_res$cover),
  MSE_1 = c(mse_1_gam, mse_1_mgwr, mse_1_bsgl, mse_1_gs),
  MSE_0 = c(mse_0_gam, mse_0_mgwr, mse_0_bsgl, mse_0_gs)
)

out_dir <- "figures/reproduction_outputs/simulation/gray_four_method_n1000_p5_rerun"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(summary_df, file.path(out_dir, "gray_run_summary.csv"), row.names = FALSE)
saveRDS(
  list(
    bsgl_fit = bsgl_fit,
    gs_fit = gs_fit,
    gam_model = gam_model,
    mgwr_model = mgwr_model,
    grid_data = grid_data,
    meta_info = meta_info,
    summary = summary_df,
    scp_95 = scp_95
  ),
  file.path(out_dir, "gray_run_fit.rds")
)

for (j in c(3, 5)) {
  plot_comp(
    bsgl_fit = bsgl_fit,
    gs_fit = gs_fit,
    gam_model = gam_model,
    mgwr_model = mgwr_model,
    grid_data = grid_data,
    meta_info = meta_info,
    j = j,
    n_val = n,
    p_val = p,
    grayscale = TRUE,
    show_title = TRUE,
    save_plot = TRUE,
    save_dir = out_dir
  )
}

print(summary_df)
print(scp_95)
