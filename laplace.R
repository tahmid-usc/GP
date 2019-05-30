# RBF kernel

ker <- function(x, l, sigf) {
  rbf <- laplacedot(sigma = 1/l)
  return(sigf * kernelMatrix(rbf, x = x))
}