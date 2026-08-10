source(file.path("R", "beta_functions.R"))
source(file.path("R", "load_data.R"))
source(file.path("R", "helpers_gpprior.R"))
source(file.path("R", "fit_chain.R"))
source(file.path("R", "calc_beta.R"))
source(file.path("R", "beta_mse.R"))

library(dplyr)
library(fields)
library(mvtnorm)
library(MCMCpack)
library(splines)
library(GIGrvg)

# =========================================================
# settings
# =========================================================
n_list <- c(1000, 2000, 5000, 10000)
p_list <- c(5, 7, 10)

L_fixed <- 25   # 如果你之前 Gaussian SVC 用 36，就改成 36

results_gs_all <- data.frame()

# =========================================================
# run Gaussian SVC only
# =========================================================
for (n_now in n_list) {
  for (p_now in p_list) {
    
    cat("\n====================================\n")
    cat("Running Gaussian SVC: n =", n_now, ", p =", p_now, "\n")
    cat("====================================\n")
    
    data <- load_data(n = n_now, p = p_now, save_dir = "data/data_indep")
    
    train_data <- data$train
    test_data  <- data$test
    
    n <- nrow(train_data$X)
    p <- ncol(train_data$X)
    
    true_beta <- calc_true_beta(train_data$coords, p)
    
    train_grid_data <- list(
      grid_points = train_data$coords
    )
    
    gs_fit <- fit_gs(
      y = train_data$y,
      X = train_data$X,
      coords = train_data$coords,
      L = L_fixed,
      niter = 5000,
      nburn = 500,
      a_k = 2,
      b_k = 1,
      a_s = 2,
      b_s = 1,
      verbose = TRUE
    )
    
    gs_res <- eval_gs(
      gs_fit,
      test_data = test_data
    )
    
    gs_beta <- betas_gs(
      gs_fit,
      train_grid_data
    )
    
    mse_1_gs   <- mean((true_beta[, 1:3] - gs_beta[, 1:3])^2)
    mse_0_gs   <- mean((true_beta[, 4:p] - gs_beta[, 4:p])^2)
    mse_avg_gs <- mean((true_beta - gs_beta)^2)
    
    one_res <- data.frame(
      n = n,
      p = p,
      Method = "Gaussian SVC",
      L = L_fixed,
      MSPE = gs_res$mspe,
      Coverage = gs_res$cover,
      MSE_1 = mse_1_gs,
      MSE_0 = mse_0_gs,
      MSE_avg = mse_avg_gs
    )
    
    print(one_res)
    
    results_gs_all <- rbind(results_gs_all, one_res)
    
    dir.create("results_gs", showWarnings = FALSE)
    write.csv(
      results_gs_all,
      file = "results_gs/gaussian_svc_all_results.csv",
      row.names = FALSE
    )
  }
}

cat("\n========== All Gaussian SVC Results ==========\n")
print(results_gs_all)

results_gs_round <- results_gs_all
num_cols <- sapply(results_gs_round, is.numeric)
results_gs_round[num_cols] <- lapply(results_gs_round[num_cols], round, 4)

cat("\n========== Rounded Results ==========\n")
print(results_gs_round)
