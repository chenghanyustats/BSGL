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
single_dir <- file.path(out_dir, "single_runs")
summary_dir <- file.path(out_dir, "summary")

dir.create(single_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

n_vec <- c(1000, 5000, 10000)
p_vec <- c(5, 10, 20)
rep_vec <- 1:10

# =========================================================
# helper
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
# run MGWR for one replicate
# =========================================================
run_mgwr_one <- function(n_val, p_val, rep_id) {
  
  full_file <- file.path(
    single_dir,
    sprintf("result_full_n%d_p%d_rep%02d.rds", n_val, p_val, rep_id)
  )
  
  old_full_file <- file.path(
    single_dir,
    sprintf("result_n%d_p%d_rep%02d.rds", n_val, p_val, rep_id)
  )
  
  nomgwr_file <- file.path(
    single_dir,
    sprintf("result_nomgwr_n%d_p%d_rep%02d.rds", n_val, p_val, rep_id)
  )
  
  partial_after_gs_file <- file.path(
    single_dir,
    sprintf("partial_n%d_p%d_rep%02d_after_gs.rds", n_val, p_val, rep_id)
  )
  
  if (file.exists(full_file)) {
    cat(sprintf(
      "Skipping existing full result: n=%d, p=%d, rep=%d\n",
      n_val, p_val, rep_id
    ))
    return(readRDS(full_file))
  }
  
  # If old full result already exists and contains MGWR, reuse it
  if (file.exists(old_full_file)) {
    old_res <- readRDS(old_full_file)
    
    if ("MGWR" %in% old_res$Method) {
      cat(sprintf(
        "Reusing old full result with MGWR: n=%d, p=%d, rep=%d\n",
        n_val, p_val, rep_id
      ))
      
      saveRDS(old_res, full_file)
      write.csv(old_res, sub("\\.rds$", ".csv", full_file), row.names = FALSE)
      return(old_res)
    }
  }
  
  # Find base result without MGWR
  if (file.exists(nomgwr_file)) {
    base_res <- readRDS(nomgwr_file)
    cat(sprintf(
      "Using no-MGWR result: n=%d, p=%d, rep=%d\n",
      n_val, p_val, rep_id
    ))
  } else if (file.exists(partial_after_gs_file)) {
    base_res <- readRDS(partial_after_gs_file)
    base_res <- base_res[base_res$Method != "MGWR", ]
    cat(sprintf(
      "Using partial after GS: n=%d, p=%d, rep=%d\n",
      n_val, p_val, rep_id
    ))
  } else if (file.exists(old_full_file)) {
    base_res <- readRDS(old_full_file)
    base_res <- base_res[base_res$Method != "MGWR", ]
    cat(sprintf(
      "Using old result after dropping MGWR: n=%d, p=%d, rep=%d\n",
      n_val, p_val, rep_id
    ))
  } else {
    message(sprintf(
      "No base result found for n=%d, p=%d, rep=%d. Skipping.",
      n_val, p_val, rep_id
    ))
    return(NULL)
  }
  
  if ("MGWR" %in% base_res$Method) {
    results_full <- base_res
    saveRDS(results_full, full_file)
    write.csv(results_full, sub("\\.rds$", ".csv", full_file), row.names = FALSE)
    return(results_full)
  }
  
  # Load data
  rep_file <- file.path(data_dir, sprintf("n%d_p%d.RData", n_val, p_val))
  
  if (!file.exists(rep_file)) {
    message("Cannot find data file: ", rep_file)
    return(NULL)
  }
  
  load(rep_file)  # loads data_list
  
  if (!exists("data_list")) {
    message("data_list not found in ", rep_file)
    return(NULL)
  }
  
  data <- data_list[[rep_id]]
  
  train_data <- data$train
  test_data  <- data$test
  meta_info  <- data$meta
  
  n <- nrow(train_data$X)
  p <- ncol(train_data$X)
  
  true_beta <- calc_true_beta(train_data$coords, p)
  seed_now <- if (!is.null(meta_info$seed)) meta_info$seed else NA
  
  cat("\n=========================================================\n")
  cat(sprintf("Running MGWR only: n=%d, p=%d, rep=%d\n", n_val, p_val, rep_id))
  cat("=========================================================\n")
  
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
    
    mgwr_beta <- get_mgwr_betas(
      mgwr_model,
      train_data$coords
    )
    
    data.frame(
      n = n,
      p = p,
      rep_id = rep_id,
      seed = seed_now,
      Method = "MGWR",
      MSPE = mean((test_data$y - pred_m)^2),
      Coverage = NA_real_,
      MSE_1 = mean((true_beta[, 1:3] - mgwr_beta[, 1:3])^2),
      MSE_0 = mean((true_beta[, 4:p] - mgwr_beta[, 4:p])^2)
    )
  })
  
  if (is.null(mgwr_out)) {
    message(sprintf(
      "MGWR failed for n=%d, p=%d, rep=%d",
      n_val, p_val, rep_id
    ))
    return(base_res)
  }
  
  results_full <- bind_rows(
    base_res[base_res$Method != "MGWR", ],
    mgwr_out
  )
  
  saveRDS(results_full, full_file)
  write.csv(results_full, sub("\\.rds$", ".csv", full_file), row.names = FALSE)
  
  cat("\n===== Finished MGWR replicate =====\n")
  print(results_full)
  cat("Saved:", full_file, "\n")
  
  return(results_full)
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
        run_mgwr_one(n_val, p_val, rep_id)
      })
      
      if (!is.null(res)) {
        all_results[[counter]] <- res
        counter <- counter + 1
      }
      
      if (length(all_results) > 0) {
        combined_now <- bind_rows(all_results)
        
        saveRDS(
          combined_now,
          file.path(summary_dir, "rep10_results_combined_running_full.rds")
        )
        
        write.csv(
          combined_now,
          file.path(summary_dir, "rep10_results_combined_running_full.csv"),
          row.names = FALSE
        )
      }
    }
  }
}

# =========================================================
# collect all full results
# =========================================================
result_files <- list.files(
  single_dir,
  pattern = "^result_full_n.*_p.*_rep.*\\.rds$",
  full.names = TRUE
)

all_done <- bind_rows(lapply(result_files, readRDS))

saveRDS(
  all_done,
  file.path(summary_dir, "rep10_results_all_done_full.rds")
)

write.csv(
  all_done,
  file.path(summary_dir, "rep10_results_all_done_full.csv"),
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
  file.path(summary_dir, "rep10_summary_mean_sd_full.rds")
)

write.csv(
  summary_mean_sd,
  file.path(summary_dir, "rep10_summary_mean_sd_full.csv"),
  row.names = FALSE
)

cat("\n=========================================================\n")
cat("All available full results summarized.\n")
cat("Summary saved to:\n")
cat(file.path(summary_dir, "rep10_summary_mean_sd_full.csv"), "\n")
cat("=========================================================\n")

print(summary_mean_sd)
