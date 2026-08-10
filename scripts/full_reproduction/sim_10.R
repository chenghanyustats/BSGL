rm(list = ls())

# =========================================================
# source files
# =========================================================
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

# =========================================================
# settings
# =========================================================
data_dir <- "data/data_rep10"

out_dir <- "results/simulation/full_rep10_run"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

single_dir <- file.path(out_dir, "single_runs")
if (!dir.exists(single_dir)) dir.create(single_dir, recursive = TRUE)

summary_dir <- file.path(out_dir, "summary")
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)

n_vec <- c(1000, 5000, 10000)
p_vec <- c(5, 10, 20)
rep_vec <- 1:10

# final-quality settings
cv_L_range <- c(25, 36)
cv_a_lambda_range <- c(15, 30, 35)
cv_b_lambda_range <- c(0.1, 1, 2)
cv_niter <- 800
cv_nburn <- 100
cv_k_folds <- 5

mcmc_niter <- 5000
mcmc_nburn <- 500

# set TRUE if you want to skip MGWR first
skip_mgwr <- FALSE

# =========================================================
# helper: safe wrapper
# =========================================================
safe_run <- function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      message("ERROR: ", conditionMessage(e))
      return(NULL)
    }
  )
}

# =========================================================
# helper: run one replicate
# =========================================================
run_one_rep <- function(n_val, p_val, rep_id) {
  
  result_file <- file.path(
    single_dir,
    sprintf("result_n%d_p%d_rep%02d.rds", n_val, p_val, rep_id)
  )
  
  if (file.exists(result_file)) {
    cat(sprintf("Skipping existing result: n=%d, p=%d, rep=%d\n",
                n_val, p_val, rep_id))
    return(readRDS(result_file))
  }
  
  cat("\n=========================================================\n")
  cat(sprintf("Running n=%d, p=%d, rep=%d\n", n_val, p_val, rep_id))
  cat("=========================================================\n")
  
  rep_file <- file.path(data_dir, sprintf("n%d_p%d.RData", n_val, p_val))
  
  if (!file.exists(rep_file)) {
    stop(sprintf("Cannot find data file: %s", rep_file))
  }
  
  load(rep_file)  # loads data_list
  
  if (!exists("data_list")) {
    stop(sprintf("data_list not found in %s", rep_file))
  }
  
  if (length(data_list) < rep_id) {
    stop(sprintf("Requested rep_id=%d but only %d reps found in %s",
                 rep_id, length(data_list), rep_file))
  }
  
  data <- data_list[[rep_id]]
  
  train_data <- data$train
  test_data  <- data$test
  grid_data  <- data$grid
  meta_info  <- data$meta
  
  n <- nrow(train_data$X)
  p <- ncol(train_data$X)
  
  true_beta <- calc_true_beta(train_data$coords, p)
  
  train_grid_data <- list(
    grid_points = train_data$coords
  )
  
  seed_now <- if (!is.null(meta_info$seed)) meta_info$seed else NA
  
  # =========================================================
  # GGP-GAM
  # =========================================================
  cat("\n===== GGP-GAM =====\n")
  
  gam_out <- safe_run({
    gam_model <- fit_gam(
      train_data$X,
      train_data$y,
      train_data$coords
    )
    
    pred_gam <- pred_gam(
      gam_model,
      test_data$X,
      test_data$coords
    )
    
    mspe_gam <- mean((test_data$y - pred_gam)^2)
    
    gam_betas <- get_gam_betas(
      gam_model,
      train_data$coords
    )
    
    mse_1_gam <- mean((true_beta[, 1:3] - gam_betas[, 1:3])^2)
    mse_0_gam <- mean((true_beta[, 4:p] - gam_betas[, 4:p])^2)
    
    data.frame(
      n = n,
      p = p,
      rep_id = rep_id,
      seed = seed_now,
      Method = "GGP-GAM",
      MSPE = mspe_gam,
      Coverage = NA_real_,
      MSE_1 = mse_1_gam,
      MSE_0 = mse_0_gam
    )
  })
  
  # =========================================================
  # BSGL CV and fit
  # =========================================================
  cat("\n===== BSGL CV =====\n")
  
  best_bsgl <- safe_run({
    cv_bsgl <- cv_bsgl(
      y = train_data$y,
      X = train_data$X,
      coords = train_data$coords,
      L_range = cv_L_range,
      a_lambda_range = cv_a_lambda_range,
      b_lambda_range = cv_b_lambda_range,
      niter = cv_niter,
      nburn = cv_nburn,
      k_folds = cv_k_folds
    )
    cv_bsgl$best_params
  })
  
  if (is.null(best_bsgl)) {
    stop(sprintf("BSGL CV failed for n=%d, p=%d, rep=%d",
                 n_val, p_val, rep_id))
  }
  
  print(best_bsgl)
  
  cat("\n===== BSGL fit =====\n")
  
  bsgl_out <- safe_run({
    bsgl_fit <- fit_bsgl(
      y = train_data$y,
      X = train_data$X,
      coords = train_data$coords,
      L = best_bsgl$L,
      a_lam = best_bsgl$a_lambda,
      b_lam = best_bsgl$b_lambda,
      niter = mcmc_niter,
      nburn = mcmc_nburn,
      verbose = TRUE
    )
    
    bsgl_res <- eval_bsgl(
      bsgl_fit,
      test_data = test_data
    )
    
    bsgl_beta <- betas_bsgl(
      bsgl_fit,
      train_grid_data
    )
    
    mse_1_bsgl <- mean((true_beta[, 1:3] - bsgl_beta[, 1:3])^2)
    mse_0_bsgl <- mean((true_beta[, 4:p] - bsgl_beta[, 4:p])^2)
    
    scp_95 <- calc_scp_bsgl(
      bsgl_fit,
      grid_data,
      ci_level = 0.95
    )
    
    scp_file <- file.path(
      single_dir,
      sprintf("scp_n%d_p%d_rep%02d.rds", n_val, p_val, rep_id)
    )
    saveRDS(scp_95, scp_file)
    
    data.frame(
      n = n,
      p = p,
      rep_id = rep_id,
      seed = seed_now,
      Method = "BSGL",
      MSPE = bsgl_res$mspe,
      Coverage = bsgl_res$cover,
      MSE_1 = mse_1_bsgl,
      MSE_0 = mse_0_bsgl
    )
  })
  
  # save intermediate after BSGL, in case later methods fail
  partial_results <- bind_rows(gam_out, bsgl_out)
  partial_file <- file.path(
    single_dir,
    sprintf("partial_n%d_p%d_rep%02d_after_bsgl.rds", n_val, p_val, rep_id)
  )
  saveRDS(partial_results, partial_file)
  
  # =========================================================
  # Gaussian SVC
  # =========================================================
  cat("\n===== Gaussian SVC =====\n")
  
  gs_out <- safe_run({
    gs_fit <- fit_gs(
      y = train_data$y,
      X = train_data$X,
      coords = train_data$coords,
      L = best_bsgl$L,
      niter = mcmc_niter,
      nburn = mcmc_nburn,
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
    
    mse_1_gs <- mean((true_beta[, 1:3] - gs_beta[, 1:3])^2)
    mse_0_gs <- mean((true_beta[, 4:p] - gs_beta[, 4:p])^2)
    
    data.frame(
      n = n,
      p = p,
      rep_id = rep_id,
      seed = seed_now,
      Method = "Gaussian SVC",
      MSPE = gs_res$mspe,
      Coverage = gs_res$cover,
      MSE_1 = mse_1_gs,
      MSE_0 = mse_0_gs
    )
  })
  
  partial_results <- bind_rows(gam_out, bsgl_out, gs_out)
  partial_file <- file.path(
    single_dir,
    sprintf("partial_n%d_p%d_rep%02d_after_gs.rds", n_val, p_val, rep_id)
  )
  saveRDS(partial_results, partial_file)
  
  # =========================================================
  # MGWR
  # =========================================================
  mgwr_out <- NULL
  
  if (!skip_mgwr) {
    cat("\n===== MGWR =====\n")
    
    mgwr_out <- safe_run({
      mgwr_model <- fit_mgwr(
        train_data$X,
        train_data$y,
        train_data$coords
      )
      
      pred_m <- pred_mgwr(
        mgwr_model,
        test_data$X,
        test_data$coords
      )
      
      mspe_mgwr <- mean((test_data$y - pred_m)^2)
      
      mgwr_beta <- get_mgwr_betas(
        mgwr_model,
        train_data$coords
      )
      
      mse_1_mgwr <- mean((true_beta[, 1:3] - mgwr_beta[, 1:3])^2)
      mse_0_mgwr <- mean((true_beta[, 4:p] - mgwr_beta[, 4:p])^2)
      
      data.frame(
        n = n,
        p = p,
        rep_id = rep_id,
        seed = seed_now,
        Method = "MGWR",
        MSPE = mspe_mgwr,
        Coverage = NA_real_,
        MSE_1 = mse_1_mgwr,
        MSE_0 = mse_0_mgwr
      )
    })
  }
  
  results_compare <- bind_rows(
    gam_out,
    mgwr_out,
    bsgl_out,
    gs_out
  )
  
  saveRDS(results_compare, result_file)
  
  csv_file <- sub("\\.rds$", ".csv", result_file)
  write.csv(results_compare, csv_file, row.names = FALSE)
  
  cat("\n===== Finished one replicate =====\n")
  print(results_compare)
  cat("Saved:", result_file, "\n")
  
  return(results_compare)
}

# =========================================================
# main loop
# =========================================================
all_results <- list()
counter <- 1

for (n_val in n_vec) {
  for (p_val in p_vec) {
    for (rep_id in rep_vec) {
      res <- safe_run({
        run_one_rep(n_val, p_val, rep_id)
      })
      
      if (!is.null(res)) {
        all_results[[counter]] <- res
        counter <- counter + 1
      }
      
      # update combined file after every replicate
      if (length(all_results) > 0) {
        combined_now <- bind_rows(all_results)
        saveRDS(
          combined_now,
          file.path(summary_dir, "rep10_results_combined_running.rds")
        )
        write.csv(
          combined_now,
          file.path(summary_dir, "rep10_results_combined_running.csv"),
          row.names = FALSE
        )
      }
    }
  }
}

# =========================================================
# collect all finished results from disk
# =========================================================
result_files <- list.files(
  single_dir,
  pattern = "^result_n.*_p.*_rep.*\\.rds$",
  full.names = TRUE
)

all_done <- bind_rows(lapply(result_files, readRDS))

saveRDS(
  all_done,
  file.path(summary_dir, "rep10_results_all_done.rds")
)

write.csv(
  all_done,
  file.path(summary_dir, "rep10_results_all_done.csv"),
  row.names = FALSE
)

# =========================================================
# summary mean and sd
# =========================================================
summary_mean_sd <- all_done %>%
  group_by(n, p, Method) %>%
  summarise(
    n_rep_done = n(),
    MSPE_mean = mean(MSPE, na.rm = TRUE),
    MSPE_sd = sd(MSPE, na.rm = TRUE),
    Coverage_mean = mean(Coverage, na.rm = TRUE),
    Coverage_sd = sd(Coverage, na.rm = TRUE),
    MSE_1_mean = mean(MSE_1, na.rm = TRUE),
    MSE_1_sd = sd(MSE_1, na.rm = TRUE),
    MSE_0_mean = mean(MSE_0, na.rm = TRUE),
    MSE_0_sd = sd(MSE_0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n, p, Method)

saveRDS(
  summary_mean_sd,
  file.path(summary_dir, "rep10_summary_mean_sd.rds")
)

write.csv(
  summary_mean_sd,
  file.path(summary_dir, "rep10_summary_mean_sd.csv"),
  row.names = FALSE
)

cat("\n=========================================================\n")
cat("All available results summarized.\n")
cat("Summary saved to:\n")
cat(file.path(summary_dir, "rep10_summary_mean_sd.csv"), "\n")
cat("=========================================================\n")

print(summary_mean_sd)
