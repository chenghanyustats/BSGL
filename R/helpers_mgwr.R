# Main MGWR function definitions
library(GWmodel)
library(sp)

#' Fit the MGWR comparison model
#'
#' Fits a multiscale geographically weighted regression model using `GWmodel`.
#'
#' @param X Numeric covariate matrix.
#' @param y Numeric response vector.
#' @param coords Numeric coordinate matrix with two columns.
#'
#' @param ... Additional arguments passed to `GWmodel::gwr.multiscale()`.
#'
#' @return A fitted MGWR model object from `GWmodel`.
fit_mgwr <- function(X, y, coords, ...) {
  # Create SpatialPointsDataFrame
  data_df <- data.frame(y = y)
  for(i in 1:ncol(X)) {
    data_df[[paste0("X", i)]] <- X[, i]
  }
  coordinates(data_df) <- coords

  # Construct model formula
  formula_str <- paste("y ~", paste(paste0("X", 1:ncol(X)), collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # Fit MGWR model
  mgwr_model <- gwr.multiscale(formula_obj, data = data_df, adaptive = TRUE,
                               kernel = "bisquare", max.iterations = 1000,
                               criterion = "CVR", verbose = FALSE, ...)

  return(mgwr_model)
}

#' Predict from a fitted MGWR model
#'
#' @param model A fitted model returned by [fit_mgwr()].
#' @param X_new Numeric matrix of new covariates.
#' @param coords_new Numeric matrix of new coordinates.
#'
#' @return Numeric vector of predicted responses.
pred_mgwr <- function(model, X_new, coords_new) {
  # Extract training-data information
  train_coords <- coordinates(model$SDF)
  coeffs_data <- model$SDF@data

  # Extract coefficients
  intercept <- coeffs_data$Intercept
  x_cols <- grep("^X[0-9]+$", names(coeffs_data))
  train_betas <- as.matrix(coeffs_data[, x_cols, drop = FALSE])

  # Extract bandwidths
  bws <- model$GW.arguments$bws
  adaptive <- model$GW.arguments$adaptive

  X_new <- as.matrix(X_new)
  coords_new <- as.matrix(coords_new)
  n_new <- nrow(X_new)
  pred <- numeric(n_new)

  # Prediction
  for(i in 1:n_new) {
    distances <- sqrt(rowSums((train_coords - matrix(coords_new[i,],
                                                     nrow = nrow(train_coords),
                                                     ncol = ncol(coords_new),
                                                     byrow = TRUE))^2))

    if(adaptive) {
      sorted_dist <- sort(distances)
      bandwidth <- sorted_dist[bws[1]]
    } else {
      bandwidth <- bws[1]
    }

    weights <- ifelse(distances <= bandwidth,
                      (1 - (distances/bandwidth)^2)^2,
                      0)
    weights <- weights / sum(weights)

    weighted_intercept <- sum(weights * intercept)
    weighted_betas <- colSums(weights * train_betas)

    pred[i] <- weighted_intercept + sum(X_new[i, ] * weighted_betas)
  }

  return(pred)
}

get_mgwr_betas <- function(model, grid_points) {
  # Extract coefficients at training locations
  coeffs_data <- model$SDF@data
  train_coords <- coordinates(model$SDF)

  x_cols <- grep("^X[0-9]+$", names(coeffs_data))
  train_betas <- as.matrix(coeffs_data[, x_cols, drop = FALSE])

  # Interpolate at grid locations
  grid_points <- as.matrix(grid_points)
  n_grid <- nrow(grid_points)
  p <- ncol(train_betas)
  grid_betas <- matrix(0, nrow = n_grid, ncol = p)

  for(j in 1:p) {
    for(i in 1:n_grid) {
      distances <- sqrt(rowSums((train_coords - matrix(grid_points[i,],
                                                       nrow = nrow(train_coords),
                                                       ncol = ncol(grid_points),
                                                       byrow = TRUE))^2))

      distances[distances < 1e-10] <- 1e-10
      weights <- 1 / distances^2
      weights <- weights / sum(weights)

      grid_betas[i, j] <- sum(weights * train_betas[, j])
    }
  }

  return(grid_betas)
}

# # Example call
# data <- load_data(n = 5000, p = 5, save_dir = "new_data")
#
# train_data <- data$train
# test_data <- data$test
# grid_data <- data$grid
# meta_info <- data$meta
# grid_points <- grid_data$grid_points
# n <- nrow(train_data$X)
# p <- ncol(train_data$X)
# true_beta <- calc_true_beta(train_data$coords, p)

# mgwr_model <- fit_mgwr(train_data$X, train_data$y, train_data$coords)
# pred_mgwr_test <- pred_mgwr(mgwr_model, test_data$X, test_data$coords)
# mgwr_beta <- get_mgwr_betas(mgwr_model, train_data$coords)
#
# # Calculate evaluation metrics
# mspe_mgwr <- mean((test_data$y - pred_mgwr_test)^2)
# mgwr_mse <- calc_mse(true_beta, mgwr_beta, p)
# mse_1_mgwr = mgwr_mse[1]
# mse_0_mgwr = mgwr_mse[2]
# mse_avg_mgwr = mgwr_mse[3]
#
# cat("\nResults with MGWR:\n")
# cat("MSPE:", mspe_mgwr, "\n")
# cat("MSE_beta_1:", mse_1_mgwr, "\n")
# cat("MSE_beta_0:", mse_0_mgwr, "\n")
# cat("MSE_avg:", mse_avg_mgwr, "\n")
