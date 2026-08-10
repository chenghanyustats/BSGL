# Calculate beta estimation accuracy
beta_mse <- function(results, grid_data, coords_train) {
  
  # Get posterior mean of beta coefficients
  beta_mean <- apply(results$beta_samples, c(1, 2), mean)
  
  # Generate basis functions for grid points
  grid_coords <- as.matrix(grid_data$grid_points)
  Psi_grid <- make_basis(grid_coords, results$L)
  
  # Predict beta coefficients on grid
  pred_betas <- Psi_grid %*% beta_mean
  true_betas <- grid_data$true_betas
  
  # Calculate MSE for non-zero variables (variables 1-3)
  mse_1_vars <- numeric(3)
  for (j in 1:3) {
    mse_1_vars[j] <- mean((pred_betas[, j] - true_betas[, j])^2)
  }
  mse_beta_1 <- mean(mse_1_vars)
  
  # Calculate MSE for zero variables (variables 4-5)
  mse_0_vars <- numeric(2)
  for (j in 4:5) {
    mse_0_vars[j-3] <- mean((pred_betas[, j] - true_betas[, j])^2)
  }
  mse_beta_0 <- mean(mse_0_vars)
  
  # Calculate average MSE
  mse_beta_avg <- mean(c(mse_1_vars, mse_0_vars))
  
  # Display results
  cat("Beta Estimation Accuracy:\n")
  cat("MSE for non-zero variables (1-3):", round(mse_beta_1, 6), "\n")
  cat("MSE for zero variables (4-5):", round(mse_beta_0, 6), "\n")
  cat("Average MSE:", round(mse_beta_avg, 6), "\n")
  
  cat("\nIndividual MSEs:\n")
  for (j in 1:5) {
    if (j <= 3) {
      cat(sprintf("Variable %d (non-zero): %.6f\n", j, mse_1_vars[j]))
    } else {
      cat(sprintf("Variable %d (zero): %.6f\n", j, mse_0_vars[j-3]))
    }
  }
  
  return(list(
    mse_beta_1 = mse_beta_1,
    mse_beta_0 = mse_beta_0, 
    mse_beta_avg = mse_beta_avg,
    individual_mse = c(mse_1_vars, mse_0_vars),
    pred_betas = pred_betas,
    true_betas = true_betas
  ))
}