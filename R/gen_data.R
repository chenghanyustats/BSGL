
library(splines)

gen_data <- function(n, p, sigma2 = 1, seed = NULL, 
                     beta_functions = list(beta1, beta2, beta3),
                     save_dir = "original_data",
                     cor = FALSE) {  # Added parameter

  if (!dir.exists(save_dir)) dir.create(save_dir)

  exp_dir <- file.path(save_dir, sprintf("n%d_p%d", n, p))
  if (!dir.exists(exp_dir)) dir.create(exp_dir)

  if (!is.null(seed)) set.seed(seed)

  coords_x <- runif(n, 0, 20)
  coords_y <- runif(n, 0, 20)

  make_x <- function(n) {
    x <- rnorm(n)
    x <- (x - min(x))/(max(x) - min(x))
    return(x)
  }

  # Generate predictor matrix
  X <- matrix(0, n, p)

  if(cor) {
    S <- matrix(0.6, 3, 3)
    diag(S) <- 1
    X3 <- rmvnorm(n, sigma = S)
    for(j in 1:3) {
      X[,j] <- (X3[,j] - min(X3[,j])) / (max(X3[,j]) - min(X3[,j]))
    }
    if(p > 3) {
      for(j in 4:p) X[,j] <- make_x(n)
    }
  } else {
    for(j in 1:p) {
      X[,j] <- make_x(n)
    }
  }

  y <- numeric(n)
  for(j in 1:min(p, length(beta_functions))) {
    beta_values <- beta_functions[[j]](coords_x, coords_y)
    y <- y + beta_values * X[,j]
  }

  e <- rnorm(n, 0, sqrt(sigma2))
  y <- y + e

  train_coords <- cbind(coords_x / 20, coords_y / 20)

  # Test data
  n_test <- floor(n * 0.2)
  coords_x_test <- runif(n_test, 0, 20)
  coords_y_test <- runif(n_test, 0, 20)

  X_test <- matrix(0, n_test, p)

  if(cor) {
    X3_test <- rmvnorm(n_test, sigma = S)
    for(j in 1:3) {
      X_test[,j] <- (X3_test[,j] - min(X3_test[,j])) / (max(X3_test[,j]) - min(X3_test[,j]))
    }
    if(p > 3) {
      for(j in 4:p) X_test[,j] <- make_x(n_test)
    }
  } else {
    for(j in 1:p) {
      X_test[,j] <- make_x(n_test)
    }
  }

  y_test <- numeric(n_test)
  for(j in 1:min(p, length(beta_functions))) {
    beta_values_test <- beta_functions[[j]](coords_x_test, coords_y_test)
    y_test <- y_test + beta_values_test * X_test[,j]
  }

  e_test <- rnorm(n_test, 0, sqrt(sigma2))
  y_test <- y_test + e_test

  test_coords <- cbind(coords_x_test / 20, coords_y_test / 20)

  # Grid data
  grid_size <- 50
  x_grid <- seq(0, 1, length.out = grid_size)
  y_grid <- seq(0, 1, length.out = grid_size)

  grid_points <- expand.grid(x = x_grid, y = y_grid)

  true_betas <- matrix(0, nrow = grid_size^2, ncol = p)
  for(j in 1:min(p, length(beta_functions))) {
    grid_coords_20 <- expand.grid(x = x_grid * 20, y = y_grid * 20)
    true_betas[,j] <- beta_functions[[j]](grid_coords_20$x, grid_coords_20$y)
  }

  train_data <- list(y = y, X = X, coords = train_coords)
  test_data <- list(y = y_test, X = X_test, coords = test_coords)
  grid_data <- list(grid_points = grid_points, true_betas = true_betas)
  meta_info <- list(n = n, p = p, sigma2 = sigma2, seed = seed)

  save(train_data, file = file.path(exp_dir, "train_data.RData"))
  save(test_data, file = file.path(exp_dir, "test_data.RData"))
  save(grid_data, file = file.path(exp_dir, "grid_data.RData"))
  save(meta_info, file = file.path(exp_dir, "meta_info.RData"))

  return(exp_dir)
}


# Other correlation structures
# # Option 1: stronger correlation (0.8)
# S <- matrix(0.8, 3, 3)
# diag(S) <- 1
#
# # Option 2: AR(1) structure - stronger correlation between neighboring variables
# rho <- 0.6
# S <- matrix(0, 3, 3)
# for(i in 1:3) {
#   for(j in 1:3) {
#     S[i,j] <- rho^abs(i-j)
#   }
# }
# # Result:
# #      [,1] [,2] [,3]
# # [1,]  1.0  0.6  0.36
# # [2,]  0.6  1.0  0.6
# # [3,] 0.36  0.6  1.0
#
# # Option 3: only X1 and X2 are correlated; X3 is independent
# S <- diag(3)
# S[1,2] <- S[2,1] <- 0.7
# # Result:
# #      [,1] [,2] [,3]
# # [1,]  1.0  0.7  0.0
# # [2,]  0.7  1.0  0.0
# # [3,]  0.0  0.0  1.0
