load_data_cluster <- function(n = 1000, p = 15, sigma = 1, seed = 123) {
  
  set.seed(seed)
  
  coords <- cbind(
    runif(n, 0, 20),
    runif(n, 0, 20)
  )
  
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta <- calc_cluster_beta(coords, p)
  
  y <- rowSums(X * beta) + rnorm(n, sd = sigma)
  
  n_test <- n
  test_coords <- cbind(
    runif(n_test, 0, 20),
    runif(n_test, 0, 20)
  )
  test_X <- matrix(rnorm(n_test * p), nrow = n_test, ncol = p)
  test_beta <- calc_cluster_beta(test_coords, p)
  test_y <- rowSums(test_X * test_beta) + rnorm(n_test, sd = sigma)
  
  grid_x <- seq(0, 20, length.out = 50)
  grid_y <- seq(0, 20, length.out = 50)
  grid_points <- as.matrix(expand.grid(grid_x, grid_y))
  true_betas <- calc_cluster_beta(grid_points, p)
  
  train_data <- list(
    X = X,
    y = y,
    coords = coords,
    true_beta = beta
  )
  
  test_data <- list(
    X = test_X,
    y = test_y,
    coords = test_coords,
    true_beta = test_beta
  )
  
  grid_data <- list(
    grid_points = grid_points,
    true_betas = true_betas
  )
  
  return(list(
    train = train_data,
    test = test_data,
    grid = grid_data
  ))
}