# Define spatially varying coefficient functions
beta1 <- function(x, y) {
  return(20 * cos(pi*x/20) * cos(pi*y/20))
}

beta2 <- function(x, y) {
  return(18 * cos(pi*x/18) * sin(pi*y/18))
}

beta3 <- function(x, y) {
  return(20 * exp(-((x-10)^2 + (y-10)^2)/50))  # Gaussian surface with a high central value
}


# calc_cluster_beta <- function(coords, p) {
#
#   x <- coords[, 1]
#   y <- coords[, 2]
#   n <- nrow(coords)
#
#   beta <- matrix(0, nrow = n, ncol = p)
#
#   region <- rep(NA, n)
#   region[x <= 10 & y <= 10] <- 1
#   region[x > 10  & y <= 10] <- 2
#   region[x <= 10 & y > 10]  <- 3
#   region[x > 10  & y > 10]  <- 4
#
#   beta[, 1] <- c(1.5, -1.0, 0.8, -0.6)[region]
#   beta[, 2] <- c(-1.2, 0.7, -0.8, 1.3)[region]
#   beta[, 3] <- c(1.0, 1.0, -1.0, -1.0)[region]
#
#   return(beta)
# }


calc_cluster_beta <- function(coords, p) {

  x <- coords[, 1]
  y <- coords[, 2]
  n <- nrow(coords)

  beta <- matrix(0, nrow = n, ncol = p)

  region <- rep(NA, n)
  region[x <= 10 & y <= 10] <- 1
  region[x > 10  & y <= 10] <- 2
  region[x <= 10 & y > 10]  <- 3
  region[x > 10  & y > 10]  <- 4

  # beta1: left half constant, right half smooth bump
  beta[, 1] <- NA

  beta[x <= 10 & y <= 10, 1] <- 1.5
  beta[x <= 10 & y > 10, 1]  <- 0.8

  beta[x > 10, 1] <- 1.5 * exp(-((x[x > 10] - 15)^2 + (y[x > 10] - 10)^2) / 30)

  # beta2 and beta3 remain piecewise constant
  beta[, 2] <- c(-1.2, 0.7, -0.8, 1.3)[region]
  beta[, 3] <- c(1.0, 1.0, -1.0, -1.0)[region]

  return(beta)
}
