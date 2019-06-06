library(mvtnorm)
library(kernlab)

Hyper <- function(x, y) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    theta <- theta^2
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, sigma = k + theta[3] * diag(n), log = T)
  }
  
  hyp <- optim(par=rep(1, 3), fn = marlik, method = 'BFGS',
               control=list(maxit = 10000))
  print(hyp)
  return(hyp$par^2)
  
}

gpfit <- function(x, y, ker = 'rbf') {
  
  if(ker == 'rbf') source("rbf.R")
  if(ker == 'laplace') source("laplace.R")
  
  x <- as.matrix(x)
  n <- dim(x)[1]
  
  #mx <- apply(x, 2, mean)
  #sdx <- apply(x, 2, sd)
  my <- mean(y)
  sdy <- sd(y)
  
  #x <- t((t(as.matrix(x)) - mx) / sdx)
  y <- (y- my)/sdy
  
  theta <- Hyper(x, y)
  
  return(list(x = x, y=y, theta = theta, ker = ker, my = my, sdy = sdy))
}

predict.gp <- function(gp, x) {
  
  if(gp$ker == 'rbf') source("rbf.R")
  if(gp$ker == 'laplace') source("laplace.R")
  
  x <- as.matrix(x)
  theta <- gp$theta
  n <- dim(gp$x)[1]
  nx <- dim(x)[1]
  #x <- t((t(as.matrix(x)) - gp$mx) / gp$sdx)

  kx <- ker2(x = x, y = gp$x, l = theta[1], sigf = theta[2])
  kxx <- ker(x = x, l = theta[1], sigf = theta[2]) + theta[3] * diag(nx)
  k <- ker(x = gp$x, l = theta[1], sigf = theta[2]) + theta[3] * diag(n)
  kinv <- chol2inv(chol(k))
  mu <- kx %*% (kinv %*% as.matrix(gp$y))
  mu <- mu * gp$sdy + gp$my
  sigma <- kxx - kx %*% (kinv %*% t(kx))
  #diag(sigma)[diag(sigma)<0] <- 0
  ll <- mu - 1.96 * (sqrt(diag(sigma)) * gp$sdy )
  ul <- mu + 1.96 * (sqrt(diag(sigma)) * gp$sdy)
  return(list(mu = mu, sigma = sigma, ll = ll, ul = ul))
  
}
