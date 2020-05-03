rm(list=ls())
library(readr)
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)
library(readxl)

source('repgrid.R')
source('rbf.R')

#yeast data example
yeast <- read.delim("data/yeast.txt")
gene <- read_excel("data/gene.xlsx")
names(gene)[1] <- "X"
mdt <- merge(yeast, gene, by = "X")
phase <- ifelse(mdt$Peak == 'G1', 'G1', 'NonG1')
mdt <- cbind(mdt, phase)
mdt <- mdt[,-(2:7)]
mdt <- mdt[, -(20:78)]
mdt <- na.omit(mdt)

alpha <- mdt[,2:19]
alpha <- matrix(t(alpha), ncol = 1)
mdt$X <- as.character(mdt$X)
gene <- rep(mdt$X, each = 18)
phase <- rep(mdt$phase, each = 18)
time <- seq(0,100, length.out = 18)
time <- rep(time, times = dim(mdt)[1])

ydt <- data.frame(gene, time, phase, alpha)
names(ydt) <- c('gene', 'time', 'phase', 'alpha')

plot(ydt[ydt$phase == 'G1',]$time, ydt[ydt$phase == 'G1',]$alpha)
plot(ydt[ydt$phase == 'NonG1',]$time, ydt[ydt$phase == 'NonG1',]$alpha)

x <- mdt[mdt$phase == 'G1',]
x <- x[,2:19]

gGP <- fit.gGP(x, tlim = c(0, 100))

xs <- seq(0, 100, length.out = 100)
gppred <- predict.gp(gGP, xs)


plot(ydt[ydt$phase == 'G1',]$time, ydt[ydt$phase == 'G1',]$alpha, pch = 16, cex = .5)
lines(xs, gppred$mu, type = 'l', lwd = 3, col = 2)
ll <- gppred$ll
ul <- gppred$ul
lines(xs, ll, type = 'l', lwd = 2, lty = 2)
lines(xs, ul, type = 'l', lwd = 2, lty = 2)


#------------------------------------------------------

# Phoneme data

phoneme <- read.csv("data/phoneme.data")

x <- phoneme[phoneme$g == 'aa',2:257]
#x <- (x - mean(x)) / sd(x)

gGP <- fit.gGP(x, tlim = c(1, 256), maxit = 1000)

xs <- seq(1, 256, length.out = 100)
gppred <- predict.gp(gGP, xs)

t <- seq(1, 256, length.out = dim(x)[2])

plot(t, x[1,], pch = 16, cex = .3, ylim = c(min(x), max(x)))
for (i in 1:dim(x)[1]) { points(t, x[i,], pch = 16, cex = .3)  }
lines(xs, gppred$mu, type = 'l', lwd = 3, col = 2)
lines(xs, gppred$ll, type = 'l', lwd = 2, col = 2, lty = 2)
lines(xs, gppred$ul, type = 'l', lwd = 2, col = 2, lty = 2)
