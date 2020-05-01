library(mvtnorm)
library(kernlab)
library(optimx)





# RBF kernel

ker <- function(x, l, sigf) {
  rbf <- rbfdot(sigma = 1/l^2)
  return(sigf^2 * kernelMatrix(rbf, x = x))
}

ker2 <- function(x, y, l, sigf) {
  rbf <- rbfdot(sigma = 1/l^2)
  return(sigf^2 * kernelMatrix(rbf, x = x, y = y))
}




mu <- function(t, a = 1, b0 = 1, b1 = 1) {
  a / (1 + b0 * exp(- b1 * t))^2
}

t <- seq(-10,10, length.out = 100)
mut <- mu(t,a = 1, b0 = 2, b1 = .7)
plot(t, mut, type = 'l')

covmat <- ker(t, l = 2, sigf = 1)
ft <- rmvnorm(1, mut, covmat)
plot(t, ft, type = 'l')

covmat <- ker(t, l = 10, sigf = .01)
gt <- rmvnorm(1, mean = rep(0, 100), sigma =  covmat)
plot(t, gt, type = 'l')

genY <- mu(t,a = 1, b0 = 2, b1 = .7) + gt + rnorm(100, mean  = 0 , sd = .01)
plot(t, genY)





Hyper <- function(x, y) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    #theta <- theta^2
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, mean = mu(x, a = theta[4], b0 = theta[5], b1 = theta[6]),  sigma = k + theta[3]^2 * diag(n), log = T)
  }
  
  hyp <- optim(par=c(1,1,1,1,1,1), fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 10000))
  print(hyp)
  return(hyp$par)
  
}

#--- Multistart
parseq <- c(.01,.1,1,2,3,5,10)
parmat <- expand.grid(parseq, parseq, parseq, parseq, parseq, parseq)

Hyper.ms <- function(x, y) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    #theta <- theta^2
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, mean = mu(x, a = theta[4], b0 = theta[5], b1 = theta[6]),  sigma = k + theta[3]^2 * diag(n), log = T)
  }
  
  hyp <- multistart(parmat=parmat, fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 10000))
  #print(hyp)
  return(hyp)
  
}

#---


theta <- Hyper(x = t, y = genY)
theta <- Hyper.ms(x = t, y = genY)




tstar <- seq(-15,15, length.out = 100)
mutstar <- mu(tstar,a = 1, b0 = 2, b1 = .7)
pred <- mu(t,a = theta[4], b0 = theta[5], b1 = theta[6])
plot(t, pred, type = 'l', lwd = 2, col = 2, ylim = c(0,1), main = 'Mean function estimation')
lines(t, mutstar, lwd = 2)

n <- length(t)
nx <- length(tstar)
kx <- ker2(x = tstar, y = t, l = theta[1], sigf = theta[2])
kxx <- ker(x = tstar, l = theta[1], sigf = theta[2]) + theta[3]^2 * diag(nx)
k <- ker(x = t, l = theta[1], sigf = theta[2]) + theta[3]^2 * diag(n)
kinv <- chol2inv(chol(k))
posmu <- kx %*% (kinv %*% matrix(genY - pred, ncol = 1))
posmu <- pred + posmu

possigma <- kxx - kx %*% (kinv %*% t(kx))
ll <- posmu + 1.96 * sqrt(diag(possigma))
ul <- posmu - 1.96 * sqrt(diag(possigma))


plot(tstar, posmu, type = 'l', lwd = 2, col = 2, main = 'Posterior mean')
lines(tstar, ll, lty = 2)
lines(tstar, ul, lty = 2)
points(t, genY)
lines(t, mutstar, lwd = 2)
