# 2. GAM (Generalized Additive Model)
library(mgcv)
# GAM fitting
#' Fit the GGP-GAM comparison model
#'
#' Fits a generalized additive model with Gaussian-process smooth terms by
#' covariate over spatial coordinates.
#'
#' @param X Numeric covariate matrix.
#' @param y Numeric response vector.
#' @param coords Numeric coordinate matrix with two columns.
#'
#' @return A fitted `mgcv::gam` object.
fit_gam <- function(X, y, coords) {

  data <- as.data.frame(X)
  colnames(data) <- paste0("X", 1:ncol(X))
  data$y <- y
  data$x <- coords[,1]
  data$y_coord <- coords[,2]

  linear_terms <- paste0("X", 1:ncol(X), collapse = " + ")
  smooth_terms <- paste0("s(x, y_coord, bs='gp', by=X", 1:ncol(X), ")", collapse = " + ")

  formula_str <- paste0("y ~ 0 + ", linear_terms, " + ", smooth_terms)
  formula_obj <- as.formula(formula_str)
  #
  # cat("GAM formula:", deparse(formula_obj), "\n")

  gam.fit <- gam(formula_obj, data = data)
  return(gam.fit)
}

# GAM prediction
#' Predict from a fitted GGP-GAM model
#'
#' @param fit A fitted model returned by [fit_gam()].
#' @param X_new Numeric matrix of new covariates.
#' @param coords_new Numeric matrix of new coordinates.
#'
#' @return Numeric vector of predicted responses.
pred_gam <- function(fit, X_new, coords_new) {
  data_new <- as.data.frame(X_new)
  colnames(data_new) <- paste0("X", 1:ncol(X_new))
  data_new$x <- coords_new[,1]
  data_new$y_coord <- coords_new[,2]

  pred <- predict(fit, newdata = data_new)
  return(pred)
}

# Extract varying coefficients
get_gam_betas <- function(model, grid_points) {

  var_names <- names(model$model)
  x_vars <- var_names[grepl("^X[0-9]+$", var_names)]
  p <- length(x_vars)
  betas <- matrix(0, nrow(grid_points), p)

  for(i in 1:p) {

    newdata <- data.frame(
      x = grid_points[,1],
      y_coord = grid_points[,2]
    )

    for(j in 1:p) {
      newdata[,paste0("X",j)] <- 0
    }

    newdata[,paste0("X",i)] <- 1

    betas[,i] <- predict(model, newdata = newdata)
  }

  return(betas)
}


metrics <- function(true_betas, pred_betas) {
  rsq = function(x, y) cor(x, y)^2
  rmse = function(x, y) sqrt(mean((x - y)^2))

  metrics <- data.frame(
    R_squared = sapply(1:ncol(true_betas), function(i) rsq(true_betas[,i], pred_betas[,i])),
    RMSE = sapply(1:ncol(true_betas), function(i) rmse(true_betas[,i], pred_betas[,i]))
  )

  return(metrics)
}


# Visualization comparison
library(cols4all)

# Example call
beta_comp <- function(true_beta, recon_beta, grid_size = 50, var_index) {
  # Create sequence for x and y
  x_seq <- seq(0, 20, length.out = grid_size)
  y_seq <- seq(0, 20, length.out = grid_size)

  # Convert vector to matrix if needed
  if(is.vector(true_beta)) {
    true_beta <- matrix(true_beta, grid_size, grid_size)
  }
  if(is.vector(recon_beta)) {
    recon_beta <- matrix(recon_beta, grid_size, grid_size)
  }

  # Check data validity
  if (all(is.na(true_beta)) || all(is.na(recon_beta))) {
    stop("Invalid beta values: all values are NA.")
  }

  # Set up plotting area - 1 row, 2 columns
  par(mfrow = c(1, 2),
      mar = c(4, 4, 3, 6),  # Increase right margin
      oma = c(0, 0, 2, 0))  # Outer margin

  # Plot true beta
  fields::image.plot(x_seq, y_seq, true_beta,
             main = paste("True beta", var_index),
             xlab = "x", ylab = "y",
             col = hcl.colors(50))

  # Plot reconstructed beta
  fields::image.plot(x_seq, y_seq, recon_beta,
             main = paste("Ggp_Gam beta", var_index),
             xlab = "x", ylab = "y",
             col = hcl.colors(50))

  # Reset plot parameters
  par(mfrow = c(1, 1))
}
