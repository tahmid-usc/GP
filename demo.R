# Demonestration

source('main.R')

x <- seq(0,1, length.out = 100)
y <- 10 * exp(- x * 4) * sin(x*8) + rnorm(100)
plot(x, y, pch = 16)

fit <- gpfit(x, y, ker = 'rbf')

plot(x, y, pch = 16)
lines(fit$newX, fit$fitted)

test <- seq(0,1, length.out = 100)
pred <- predict.gp(gp = fit, x = test)
plot(x, y, pch = 16)
lines(test, pred$mu)
lines(test,pred$ll, lty = 2)
lines(test, pred$ul, lty = 2)


# Multiple regression

x1 <- runif(100, 0, 1)
x2 <- runif(100, 0, 1)

x <- cbind(x1,x2) 
y <- 10 * exp(- x1 * 4) * sin(x2*8) + rnorm(100) 

fit <- gpfit(x, y, ker = 'rbf')

tx1 <- seq(0,1, length.out = 100)
tx2 <- seq(0,1, length.out = 100)
test <- cbind(tx1, tx2)

pred <- predict.gp(gp = fit, x = test)

ty <- 10 * exp(- tx1 * 4) * sin(tx2*8)
mean((ty - pred$mu)^2)
plot((ty - pred$mu))


library(mgcv)

x <- data.frame(x)
test <- data.frame(test)
names(test) <- c('x1', 'x2')
gamfit <- gam(y~ s(x1,x2))
pred.gam <- predict(gamfit, newdata = test)
mean((ty - pred.gam)^2)
plot((ty - pred.gam))
