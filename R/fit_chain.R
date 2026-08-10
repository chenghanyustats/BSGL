# Multi-chain fitting function with R-hat diagnostics
fit_lggp_multichain <- function(y, X, coords, L = 50, 
                                niter = 3000, nburn = 500, 
                                nchains = 4,
                                a_lambda = 0.1, b_lambda = 5,
                                parallel = TRUE) {
  
  cat("=== Multi-chain LGGP Fitting ===\n")
  cat("Chains:", nchains, "| Iterations per chain:", niter, "| Burn-in:", nburn, "\n")
  
  # Run multiple chains
  if (parallel && require(parallel)) {
    cat("Running chains in parallel...\n")
    cl <- makeCluster(min(nchains, detectCores() - 1))
    
    # Export ALL necessary objects and functions
    clusterExport(cl, c("fit_lggp", "make_basis"), 
                  envir = environment())
    
    # Load required libraries on each worker
    clusterEvalQ(cl, {
      library(mvtnorm)
      library(MCMCpack)
      library(splines)
      library(GIGrvg)
    })
    
    # Run chains with explicit parameter passing
    chain_results <- parLapply(cl, 1:nchains, function(chain_id, 
                                                       y_data, X_data, coords_data,
                                                       L_param, niter_param, nburn_param,
                                                       a_lambda_param, b_lambda_param) {
      set.seed(12345 + chain_id * 100)
      fit_lggp(y_data, X_data, coords_data, L_param, niter_param, nburn_param, 
               a_lambda_param, b_lambda_param)
    }, y_data = y, X_data = X, coords_data = coords,
    L_param = L, niter_param = niter, nburn_param = nburn,
    a_lambda_param = a_lambda, b_lambda_param = b_lambda)
    
    stopCluster(cl)
  } else {
    cat("Running chains sequentially...\n")
    chain_results <- vector("list", nchains)
    for (chain_id in 1:nchains) {
      set.seed(12345 + chain_id * 100)
      cat("\n--- Chain", chain_id, "---\n")
      chain_results[[chain_id]] <- fit_lggp(y, X, coords, L, niter, nburn, 
                                            a_lambda, b_lambda)
    }
  }
  
  # Combine results
  combined <- combine_chains(chain_results)
  
  # Compute R-hat diagnostics (only key parameters by default)
  cat("\n=== Computing R-hat diagnostics ===\n")
  diagnostics <- compute_rhat(chain_results, check_all_tau2 = FALSE)
  
  # Print diagnostics
  print_diagnostics(diagnostics)
  
  # Add diagnostics to results
  combined$diagnostics <- diagnostics
  combined$chain_results <- chain_results
  
  return(combined)
}


# Combine chains into single result
combine_chains <- function(chain_results) {
  nchains <- length(chain_results)
  
  # Dimensions
  L <- dim(chain_results[[1]]$beta_samples)[1]
  m <- dim(chain_results[[1]]$beta_samples)[2]
  nsamples_per_chain <- dim(chain_results[[1]]$beta_samples)[3]
  
  # Combine all chains
  beta_combined <- array(0, dim = c(L, m, nsamples_per_chain * nchains))
  tau2_combined <- matrix(0, nrow = nsamples_per_chain * nchains, ncol = m)
  sigma2_combined <- numeric(nsamples_per_chain * nchains)
  lambda2_combined <- numeric(nsamples_per_chain * nchains)
  
  for (chain in 1:nchains) {
    idx_start <- (chain - 1) * nsamples_per_chain + 1
    idx_end <- chain * nsamples_per_chain
    
    beta_combined[, , idx_start:idx_end] <- chain_results[[chain]]$beta_samples
    tau2_combined[idx_start:idx_end, ] <- chain_results[[chain]]$tau2_samples
    sigma2_combined[idx_start:idx_end] <- chain_results[[chain]]$sigma2_samples
    lambda2_combined[idx_start:idx_end] <- chain_results[[chain]]$lambda2_samples
  }
  
  return(list(
    beta_samples = beta_combined,
    tau2_samples = tau2_combined,
    sigma2_samples = sigma2_combined,
    lambda2_samples = lambda2_combined,
    Psi = chain_results[[1]]$Psi,
    saved_knots = chain_results[[1]]$saved_knots,
    n = chain_results[[1]]$n,
    m = chain_results[[1]]$m,
    L = chain_results[[1]]$L,
    nchains = nchains
  ))
}


# Compute R-hat for KEY parameters only
compute_rhat <- function(chain_results, check_all_tau2 = FALSE) {
  nchains <- length(chain_results)
  L <- dim(chain_results[[1]]$beta_samples)[1]
  m <- dim(chain_results[[1]]$beta_samples)[2]
  nsamples <- dim(chain_results[[1]]$beta_samples)[3]
  
  # 1. R-hat for beta (L2 norms) - ESSENTIAL
  beta_norms_chains <- array(0, dim = c(nsamples, m, nchains))
  for (chain in 1:nchains) {
    for (j in 1:m) {
      beta_norms_chains[, j, chain] <- apply(
        chain_results[[chain]]$beta_samples[, j, ], 
        2, 
        function(x) sqrt(sum(x^2))
      )
    }
  }
  
  rhat_beta <- numeric(m)
  for (j in 1:m) {
    rhat_beta[j] <- compute_rhat_single(beta_norms_chains[, j, ])
  }
  
  # 2. R-hat for sigma2 - ESSENTIAL
  sigma2_chains <- matrix(0, nrow = nsamples, ncol = nchains)
  for (chain in 1:nchains) {
    sigma2_chains[, chain] <- chain_results[[chain]]$sigma2_samples
  }
  rhat_sigma2 <- compute_rhat_single(sigma2_chains)
  
  # 3. R-hat for lambda2 - IMPORTANT (global shrinkage)
  lambda2_chains <- matrix(0, nrow = nsamples, ncol = nchains)
  for (chain in 1:nchains) {
    lambda2_chains[, chain] <- chain_results[[chain]]$lambda2_samples
  }
  rhat_lambda2 <- compute_rhat_single(lambda2_chains)
  
  # 4. R-hat for tau2 - OPTIONAL (only if user requests)
  rhat_tau2 <- NULL
  tau2_chains <- NULL
  
  if (check_all_tau2) {
    tau2_chains <- array(0, dim = c(nsamples, m, nchains))
    for (chain in 1:nchains) {
      tau2_chains[, , chain] <- chain_results[[chain]]$tau2_samples
    }
    
    rhat_tau2 <- numeric(m)
    for (j in 1:m) {
      rhat_tau2[j] <- compute_rhat_single(tau2_chains[, j, ])
    }
  }
  
  return(list(
    rhat_beta = rhat_beta,
    rhat_sigma2 = rhat_sigma2,
    rhat_lambda2 = rhat_lambda2,
    rhat_tau2 = rhat_tau2,  # NULL if not checked
    beta_norms_chains = beta_norms_chains,
    sigma2_chains = sigma2_chains,
    lambda2_chains = lambda2_chains,
    tau2_chains = tau2_chains
  ))
}


# Compute R-hat for a single parameter (Gelman-Rubin statistic)
compute_rhat_single <- function(chains_matrix) {
  # chains_matrix: nsamples x nchains
  nsamples <- nrow(chains_matrix)
  nchains <- ncol(chains_matrix)
  
  # Chain means
  chain_means <- colMeans(chains_matrix)
  
  # Overall mean
  overall_mean <- mean(chains_matrix)
  
  # Between-chain variance
  B <- nsamples * var(chain_means)
  
  # Within-chain variance
  W <- mean(apply(chains_matrix, 2, var))
  
  # Variance estimate
  var_plus <- ((nsamples - 1) / nsamples) * W + (1 / nsamples) * B
  
  # R-hat
  rhat <- sqrt(var_plus / W)
  
  return(rhat)
}


# Print diagnostics summary (simplified)
print_diagnostics <- function(diagnostics) {
  cat("\n=== Convergence Diagnostics (R-hat) ===\n")
  cat("Threshold: R-hat < 1.1 indicates good convergence\n\n")
  
  # Essential: Beta coefficients
  cat("** Beta coefficients (L2 norms) **\n")
  for (j in 1:length(diagnostics$rhat_beta)) {
    status <- ifelse(diagnostics$rhat_beta[j] < 1.1, "✓", "✗")
    cat(sprintf("  Variable %d: %.4f %s\n", j, diagnostics$rhat_beta[j], status))
  }
  
  # Essential: Noise variance
  cat("\n** Sigma2 (noise variance) **\n")
  status <- ifelse(diagnostics$rhat_sigma2 < 1.1, "✓", "✗")
  cat(sprintf("  %.4f %s\n", diagnostics$rhat_sigma2, status))
  
  # Important: Global shrinkage
  cat("\n** Lambda2 (global shrinkage) **\n")
  status <- ifelse(diagnostics$rhat_lambda2 < 1.1, "✓", "✗")
  cat(sprintf("  %.4f %s\n", diagnostics$rhat_lambda2, status))
  
  # Optional: Variable-specific shrinkage
  if (!is.null(diagnostics$rhat_tau2)) {
    cat("\n** Tau2 (variable-specific shrinkage) **\n")
    for (j in 1:length(diagnostics$rhat_tau2)) {
      status <- ifelse(diagnostics$rhat_tau2[j] < 1.1, "✓", "✗")
      cat(sprintf("  Variable %d: %.4f %s\n", j, diagnostics$rhat_tau2[j], status))
    }
  }
  
  # Summary
  all_rhats <- c(diagnostics$rhat_beta, diagnostics$rhat_sigma2, 
                 diagnostics$rhat_lambda2)
  if (!is.null(diagnostics$rhat_tau2)) {
    all_rhats <- c(all_rhats, diagnostics$rhat_tau2)
  }
  
  max_rhat <- max(all_rhats)
  n_converged <- sum(all_rhats < 1.1)
  n_total <- length(all_rhats)
  
  cat("\n=== Summary ===\n")
  cat(sprintf("Max R-hat: %.4f\n", max_rhat))
  cat(sprintf("Converged: %d/%d (%.1f%%)\n", 
              n_converged, n_total, 100 * n_converged / n_total))
  
  if (max_rhat < 1.1) {
    cat("✓ All key parameters converged!\n")
  } else {
    cat("✗ Some parameters have not converged.\n")
    cat("  Recommendations:\n")
    cat("  - Run longer chains (increase niter)\n")
    cat("  - Check trace plots for poor mixing\n")
    cat("  - Consider thinning if autocorrelation is high\n")
  }
  
  # Practical guidance
  cat("\n=== Practical Notes ===\n")
  cat("Priority for convergence checking:\n")
  cat("  1. Beta (your main effects) - MUST converge\n")
  cat("  2. Sigma2 (uncertainty) - MUST converge\n")
  cat("  3. Lambda2 (shrinkage) - Should converge\n")
  cat("  4. Tau2 (optional) - Less critical\n")
}


# Plot trace plots for diagnostics
plot_trace <- function(diagnostics, param = "sigma2") {
  if (param == "sigma2") {
    chains <- diagnostics$sigma2_chains
    main_title <- "Trace Plot: sigma2"
    ylab <- "sigma2"
  } else if (param == "lambda2") {
    chains <- diagnostics$lambda2_chains
    main_title <- "Trace Plot: lambda2"
    ylab <- "lambda2"
  } else if (grepl("^beta", param)) {
    j <- as.numeric(sub("beta", "", param))
    chains <- diagnostics$beta_norms_chains[, j, ]
    main_title <- paste0("Trace Plot: beta", j, " (L2 norm)")
    ylab <- "L2 norm"
  }
  
  nchains <- ncol(chains)
  colors <- rainbow(nchains)
  
  plot(chains[, 1], type = "l", col = colors[1], 
       main = main_title, xlab = "Iteration", ylab = ylab,
       ylim = range(chains))
  
  for (chain in 2:nchains) {
    lines(chains[, chain], col = colors[chain])
  }
  
  legend("topright", legend = paste("Chain", 1:nchains), 
         col = colors, lty = 1, cex = 0.8)
}


# Plot all trace plots
plot_all_traces <- function(diagnostics, m) {
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  
  plot_trace(diagnostics, "sigma2")
  plot_trace(diagnostics, "lambda2")
  
  for (j in 1:min(4, m)) {
    plot_trace(diagnostics, paste0("beta", j))
  }
}


# Modified CV function with multi-chain
cv_lggp_multichain <- function(y, X, coords, 
                               L_range, a_lambda_range, b_lambda_range,
                               niter = 800, nburn = 100, 
                               nchains = 3,
                               k_folds = 5) {
  
  cat("=== Cross-validation with Multi-chain MCMC ===\n")
  
  # Create folds
  n <- length(y)
  fold_ids <- sample(rep(1:k_folds, length.out = n))
  
  # Grid search
  param_grid <- expand.grid(
    L = L_range,
    a_lambda = a_lambda_range,
    b_lambda = b_lambda_range
  )
  
  cv_scores <- numeric(nrow(param_grid))
  max_rhats <- numeric(nrow(param_grid))
  
  for (i in 1:nrow(param_grid)) {
    params <- param_grid[i, ]
    cat("\n--- Testing: L =", params$L, ", a_lambda =", params$a_lambda, 
        ", b_lambda =", params$b_lambda, "---\n")
    
    fold_mspes <- numeric(k_folds)
    fold_max_rhats <- numeric(k_folds)
    
    for (k in 1:k_folds) {
      cat("Fold", k, "... ")
      
      # Split data
      test_idx <- which(fold_ids == k)
      train_idx <- setdiff(1:n, test_idx)
      
      # Fit model with multiple chains
      fit <- fit_lggp_multichain(
        y = y[train_idx],
        X = X[train_idx, ],
        coords = coords[train_idx, ],
        L = params$L,
        niter = niter,
        nburn = nburn,
        nchains = nchains,
        a_lambda = params$a_lambda,
        b_lambda = params$b_lambda,
        parallel = FALSE  # Sequential for CV
      )
      
      # Check convergence (only key parameters)
      all_rhats <- c(fit$diagnostics$rhat_beta, 
                     fit$diagnostics$rhat_sigma2, 
                     fit$diagnostics$rhat_lambda2)
      fold_max_rhats[k] <- max(all_rhats)
      
      # Predict
      pred <- predict_lggp(fit, X[test_idx, ], coords[test_idx, ])
      fold_mspes[k] <- mean((y[test_idx] - pred$mean)^2)
      
      cat("MSPE =", round(fold_mspes[k], 3), 
          ", Max R-hat =", round(fold_max_rhats[k], 4), "\n")
    }
    
    cv_scores[i] <- mean(fold_mspes)
    max_rhats[i] <- mean(fold_max_rhats)
    
    cat("Average MSPE:", round(cv_scores[i], 3), 
        ", Average Max R-hat:", round(max_rhats[i], 4), "\n")
  }
  
  # Find best parameters (prioritize convergence)
  converged <- max_rhats < 1.1
  if (any(converged)) {
    best_idx <- which.min(cv_scores[converged])
    best_idx <- which(converged)[best_idx]
  } else {
    cat("Warning: No parameters achieved full convergence!\n")
    best_idx <- which.min(cv_scores)
  }
  
  return(list(
    best_params = param_grid[best_idx, ],
    cv_scores = cv_scores,
    max_rhats = max_rhats,
    param_grid = param_grid,
    best_idx = best_idx
  ))
}