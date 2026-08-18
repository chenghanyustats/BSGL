calc_true_beta <- function(coords, p) {
  n <- nrow(coords)
  true_beta <- matrix(0, n, p)

  # coords are in [0,1]; rescale to [0,20] to calculate beta
  coords_20 <- coords * 20

  for (i in 1:n) {
    x <- coords_20[i, 1]
    y <- coords_20[i, 2]
    true_beta[i, 1] <- beta1(x, y)
    true_beta[i, 2] <- beta2(x, y)
    true_beta[i, 3] <- beta3(x, y)
  }
  return(true_beta)
}

# calc_true_beta <- function(coords, p) {
#   n <- nrow(coords)
#   true_beta <- matrix(0, n, p)
#   for (i in 1:n) {
#     x <- coords[i, 1] * 20
#     y <- coords[i, 2] * 20
#     true_beta[i, 1] <- beta1(x, y)
#     true_beta[i, 2] <- beta2(x, y)
#     true_beta[i, 3] <- beta3(x, y)
#   }
#   return(true_beta)
# }

calc_mse <- function(true_beta, pred_beta, p) {
  mse_1 <- mean((true_beta[,1:3] - pred_beta[,1:3])^2)
  mse_0 <- mean((true_beta[,4:p] - pred_beta[,4:p])^2)
  mse_avg <- mean((true_beta - pred_beta)^2)
  return(c(mse_1, mse_0, mse_avg))
}
