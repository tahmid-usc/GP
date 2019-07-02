
fit.gGP <- function(x, tlim = c(0, 1), maxit = 1000) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    m <- dim(x)[2]
    t <- seq(tlim[1], tlim[2], length.out = m)
    theta <- theta^2
    ms <- n/theta[3]
    B <- as.matrix(apply(x, 2, sum))
    k <- ker(t, l = theta[1], sigf = theta[2]) 
    kinv <- chol2inv(chol(k))
    C <- chol2inv(chol(kinv + ms * diag(m)))
    
    logl <- (1/(2*theta[3])) 
    logl <- logl *  (sum(x^2) - (1/theta[3]) * (t(B) %*% (C %*% B)))
    dM <- (ms * k + diag(m))
    logl <- logl + .5 * log((1/theta[3]^n)) - .5 * determinant(dM, logarithm =  T)$modulus
    logl <- as.numeric(logl)
    return(logl)
  }
  
  hyp <- optim(par=rep(1, 3), fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = maxit))
  hyp <- hyp$par^2
  print(hyp)
  
  x <- as.matrix(x)
  n <- dim(x)[1]
  m <- dim(x)[2]
  t <- seq(tlim[1], tlim[2], length.out = m)
  theta <- hyp
  ms <- n/theta[3]
  B <- as.matrix(apply(x, 2, sum))
  k <- ker(t, l = theta[1], sigf = theta[2]) 
  kinv <- chol2inv(chol(k))
  C <- chol2inv(chol(kinv + ms * diag(m)))
  
    return(list(x = x, t = t, theta = theta, n = n, m = m, ms = ms, B = B, C = C))
  
}


predict.gp <- function(gp, xs) {
  
  x <- gp$x
  theta <- gp$theta
  n <- gp$n
  m <- gp$m
  ns <- length(xs)
  t <- gp$t
  ks <- ker2(xs, t, l = theta[1], sigf = theta[2])
  kss <- ker(xs, l = theta[1], sigf = theta[2])
  B <- gp$B
  ms <- gp$ms
  C <- gp$C
  mu <- (1/ theta[3]) * ks %*% (((diag(m) - ms * C) %*% B))
  sigma <- kss + theta[3] * diag(ns) - ms * ks %*% ((diag(m) - ms * C) %*% t(ks))
  
  ll <- mu - 1.96 * sqrt(diag(sigma))
  ul <- mu + 1.96 * sqrt(diag(sigma))
  
  return(list(mu = mu, sigma = sigma, ul = ul, ll = ll))
  
}
