# RBF kernel

ker <- function(x, l, sigf) {
  rbf <- laplacedot(sigma = 1/l)
  return(sigf * kernelMatrix(rbf, x = x))
}

ker2 <- function(x, y, l, sigf) {
  rbf <- laplacedot(sigma = 1/l)
  return(sigf * kernelMatrix(rbf, x = x, y = y))
}