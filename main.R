library(mvtnorm)

Hyper <- function(x, y) {
  
  marlik <- function(theta) {
    n <- length(x)
    theta <- theta^2
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, sigma = k + theta[3] * diag(n), log = T)
  }
  
  hyp <- optim(par=rep(1, 3), fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 1000))
  print(hyp)
  return(hyp$par^2)
  
}

gpfit <- function(x, y, ker = 'rbf') {
  
  if(ker == 'rbf') source("rbf.R")
  if(ker == 'laplace') source("laplace.R")
  
  n <- length(x)
  theta <- Hyper(x, y)
  newX <- seq(min(x), max(x), length.out = 10)
  kx <- ker2(x = newX, y = x, l = theta[1], sigf = theta[2])
  k <- ker(x = x, l = theta[1], sigf = theta[2]) + theta[3] * diag(n)
  kinv <- chol2inv(chol(k))
  
  mu <- kx %*% (kinv %*% y)
  
  return(list(x = x, y=y, theta = theta, fitted = mu, ker = ker))
}

predict.gp <- function(gp, x) {
  
  if(gp$ker == 'rbf') source("rbf.R")
  if(gp$ker == 'laplace') source("laplace.R")
  
  theta <- gp$theta
  n <- length(gp$x)
  kx <- ker2(x = x, y = gp$x, l = theta[1], sigf = theta[2])
  k <- ker(x = gp$x, l = theta[1], sigf = theta[2]) + theta[3] * diag(n)
  kinv <- chol2inv(chol(k))
  mu <- kx %*% (kinv %*% as.matrix(gp$y))
  return(mu)
  
}

fit <- gpfit(x, y, ker = 'rbf')

predict.gp(gp = fit, x = 2:10)



plot(fit$x, fit$y, type = 'l')
