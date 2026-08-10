
library(splines)

gen_rep <- function(n, p, sigma2 = 0.1, seed = NULL,
                    beta_functions = list(beta1, beta2, beta3),
                    cor = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  make_x <- function(n) {
    x <- rnorm(n)
    x <- (x - min(x)) / (max(x) - min(x))
    return(x)
  }
  
  coords_x <- runif(n, 0, 20)
  coords_y <- runif(n, 0, 20)
  
  X <- matrix(0, n, p)
  
  if (cor) {
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      stop("Package 'mvtnorm' is required when cor = TRUE.")
    }
    
    S <- matrix(0.6, 3, 3)
    diag(S) <- 1
    
    X3 <- mvtnorm::rmvnorm(n, sigma = S)
    for (j in 1:min(3, p)) {
      X[, j] <- (X3[, j] - min(X3[, j])) / (max(X3[, j]) - min(X3[, j]))
    }
    
    if (p > 3) {
      for (j in 4:p) X[, j] <- make_x(n)
    }
  } else {
    for (j in 1:p) X[, j] <- make_x(n)
  }
  
  y <- numeric(n)
  for (j in 1:min(p, length(beta_functions))) {
    beta_values <- beta_functions[[j]](coords_x, coords_y)
    y <- y + beta_values * X[, j]
  }
  
  y <- y + rnorm(n, 0, sqrt(sigma2))
  
  train_coords <- cbind(coords_x / 20, coords_y / 20)
  
  n_test <- floor(n * 0.2)
  coords_x_test <- runif(n_test, 0, 20)
  coords_y_test <- runif(n_test, 0, 20)
  
  X_test <- matrix(0, n_test, p)
  
  if (cor) {
    S <- matrix(0.6, 3, 3)
    diag(S) <- 1
    
    X3_test <- mvtnorm::rmvnorm(n_test, sigma = S)
    for (j in 1:min(3, p)) {
      X_test[, j] <- (X3_test[, j] - min(X3_test[, j])) /
        (max(X3_test[, j]) - min(X3_test[, j]))
    }
    
    if (p > 3) {
      for (j in 4:p) X_test[, j] <- make_x(n_test)
    }
  } else {
    for (j in 1:p) X_test[, j] <- make_x(n_test)
  }
  
  y_test <- numeric(n_test)
  for (j in 1:min(p, length(beta_functions))) {
    beta_values_test <- beta_functions[[j]](coords_x_test, coords_y_test)
    y_test <- y_test + beta_values_test * X_test[, j]
  }
  
  y_test <- y_test + rnorm(n_test, 0, sqrt(sigma2))
  
  test_coords <- cbind(coords_x_test / 20, coords_y_test / 20)
  
  grid_size <- 50
  x_grid <- seq(0, 1, length.out = grid_size)
  y_grid <- seq(0, 1, length.out = grid_size)
  grid_points <- expand.grid(x = x_grid, y = y_grid)
  
  true_betas <- matrix(0, nrow = grid_size^2, ncol = p)
  grid_coords_20 <- expand.grid(x = x_grid * 20, y = y_grid * 20)
  
  for (j in 1:min(p, length(beta_functions))) {
    true_betas[, j] <- beta_functions[[j]](
      grid_coords_20$x,
      grid_coords_20$y
    )
  }
  
  return(list(
    train = list(
      y = y,
      X = X,
      coords = train_coords
    ),
    test = list(
      y = y_test,
      X = X_test,
      coords = test_coords
    ),
    grid = list(
      grid_points = grid_points,
      true_betas = true_betas
    ),
    meta = list(
      n = n,
      p = p,
      sigma2 = sigma2,
      seed = seed,
      cor = cor
    )
  ))
}
