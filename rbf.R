# RBF kernel

ker <- function(x, l, sigf) {
  rbf <- rbfdot(sigma = 1/l)
  return(sigf * kernelMatrix(rbf, x = x))
}

ker2 <- function(x, y, l, sigf) {
  rbf <- rbfdot(sigma = 1/l)
  return(sigf * kernelMatrix(rbf, x = x, y = y))
}