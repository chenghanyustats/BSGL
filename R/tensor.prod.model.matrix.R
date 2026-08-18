
tensor.prod.model.matrix <- function(Bx, By) {
  nx <- nrow(Bx)
  ny <- nrow(By)
  px <- ncol(Bx)
  py <- ncol(By)

  # Create tensor product
  B <- matrix(0, nx, px * py)
  for (i in 1:px) {
    for (j in 1:py) {
      B[, (i-1)*py + j] <- Bx[,i] * By[,j]
    }
  }
  return(B)
}
