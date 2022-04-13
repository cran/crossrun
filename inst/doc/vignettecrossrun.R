## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE, 
  comment   = "#>",
  fig.width = 7.15,
  fig.height = 3.5,
  echo = FALSE,
  message = FALSE)

library(crossrun)

## ----fig1, fig.cap="Figure 1"-------------------------------------------------
set.seed(32)
n <- 24
y <- rnorm(n)
x <- seq(n)

op <- par(mar = c(bottom = 0, left   = 0, top    = 0, right  = 0))

plot(x, y,
     axes = FALSE,
     type = "b",
     pch  = '',
     lwd  = 1.5,
     ylab = '',
     xlab = '')

lines(x, rep(median(y), n),
      col = 'grey40')

text(x, y)

par(op)

## ----symm14, echo=FALSE, message=FALSE----------------------------------------
library(crossrun)
j14s <- joint14symm
rownames(j14s)[1] <- "C = 0"
colnames(j14s)[1] <- "L = 1"
knitr::kable(j14s, caption = 'Table 1')

## ----p14.6, echo=FALSE, message=FALSE, fig.width=15---------------------------
p14.6<- joint14.6
rownames(p14.6)[1] <- "C = 0"
colnames(p14.6)[1] <- "L = 1"
knitr::kable(round(p14.6,1), caption = 'Table 2')

## ----em14tab, echo=FALSE, message=FALSE---------------------------------------
library(crossrun)
em14 <- joint14em
knitr::kable(em14, caption = 'Table 3')

## ----emcomp14, message=FALSE--------------------------------------------------
library(crossrun)
j14s <- joint14symm
j14em <- joint14em
plot(x      = 0:13,
     y      = (cumsum(cumsumm(j14s)[,14])/(2^13)),
     type   = "l",
     xlab   = "Number of crossings, n=14",
     ylab   = "CDF",
     las    = 1)
points(x    = 1:13,
       y    = cumsum(cumsumm(j14em)[,7])/(choose(14,7)/2),
       type = "l",
       col  = "red",
       lty  = "dotted")
lines(x     = c(8, 9),
      y     = c(0.1, 0.1),
      col   = "red")
text(x      = 9,
     y      = 0.1,
     pos    = 4,
     labels = "red: empirical median",
     col    = "red")

plot(x      = 1:14,
     y      = cumsum(cumsummcol(j14s)[14,])/(2^13),
     type   = "l",
     xlab   = "Longest run, n=14",
     ylab   = "CDF",
     las    = 1)
points(x    = 1:7,
       y    = cumsum(cumsummcol(j14em)[13,])/(choose(14,7)/2),
       type = "l",
       col  = "red",
       lty  = "dotted")
lines(x   = c(8, 9),
      y   = c(0.1, 0.1),
      col = "red")
text(x      = 9,
     y      = 0.1,
     pos    = 4,
     labels = "red: empirical median",
     col    = "red")

## ----emcomp60, message=FALSE--------------------------------------------------
library(crossrun)
j60s <- joint60symm
j60em <- joint60em
plot(x      = 0:59,
     y      = (cumsum(cumsumm(j60s)[,60])/(2^59)),
     type   = "l",
     xlab   = "Number of crossings, n=60",
     ylab   = "CDF",
     las    = 1)
points(x    = 1:59,
       y    = cumsum(cumsumm(j60em)[,30])/(choose(60,30)/2),
       type = "l",
       col  = "red",
       lty  = "dotted")
lines(x     = c(30, 35),
      y     = c(0.1, 0.1),
      col   = "red")
text(x      = 35,
     y      = 0.1,
     pos    = 4,
     labels = "red: empirical median",
     col    = "red")

plot(x      = 1:60,
     y      = cumsum(cumsummcol(j60s)[60,])/(2^59),
     type   = "l",
     xlab   = "Longest run, n=60",
     ylab   = "CDF",
     las    = 1)
points(x    = 1:30,
       y    = cumsum(cumsummcol(j60em)[59,])/(choose(60,30)/2),
       type = "l",
       col  = "red",
       lty  = "dotted")
lines(x   = c(30, 35),
      y   = c(0.1, 0.1),
      col = "red")
text(x      = 35,
     y      = 0.1,
     pos    = 4,
     labels = "red: empirical median",
     col    = "red")

## ----exact1, message=FALSE, echo=T--------------------------------------------
library(crossrun)
library(Rmpfr)

exact1   <- asNumeric((2^4) * exactbin(n = 5, p = 0.6))
iter1    <- asNumeric(crossrunbin(nmax = 5, prob = 0.6)$pt[[5]])
compare1 <- cbind(exact1,iter1)

compare1

## ----bincoeff14, message=FALSE, echo=T----------------------------------------
library(crossrun)
library(Rmpfr)

bincoeff13           <- Rmpfr::chooseMpfr.all(13) # binomial coefficients, n - 1 = 13
bincoeff13iter       <- cumsumm(j14s)[-1, 14]     # row sums, n - 1 = 13
compare13            <- rbind(asNumeric(bincoeff13), bincoeff13iter)
row.names(compare13) <- c("exact","iter")
compare13
max(abs(bincoeff13 - bincoeff13iter))

## ----sim14, message=FALSE, echo=T---------------------------------------------
library(crossrun)
set.seed(83938487)
sim14 <- simclbin(nser = 14, nsim = 10000)
(matrix(0:13, nrow = 1) %*% p14.6%*% matrix(1:14, ncol = 1)) / 2^13
mean(sim14$nc0.6 * sim14$lr0.6)

## ----sim14plot, message=FALSE-------------------------------------------------
library(crossrun)
plot(x      = as.numeric(names(table(sim14$nc0.6))),
     y      = (cumsum(cumsumm(p14.6)[,14]) /
                 (2^13))[as.numeric(names(table(sim14$nc0.6))) + 1],
     type   = "l",
     xlab   = "Number of crossings",
     ylab   = "CDF",
     las    = 1)
points(x    = as.numeric(names(table(sim14$nc0.6))),
       y    = cumsum(table(sim14$nc0.6))/sum(table(sim14$nc0.6)),
       type = "l",
       col  = "red",
       lty  = "dotted")
lines(x     = c(0, 1),
      y     = c(0.9, 0.9),
      col   = "red")
text(x      = 1,
     y      = 0.9,
     pos    = 4,
     labels = "red: simulations",
     col    = "red")

plot(x      = as.numeric(names(table(sim14$lr0.6))),
     y      = as.numeric(cumsum(cumsummcol(p14.6)[14,]) /
                           sum(cumsummcol(p14.6)[14,]))[
                             as.numeric(names(table(sim14$lr0.6)))],
     type   = "l",
     xlab   = "Longest run",
     ylab   = "CDF",
     las    = 1)
points(x    = as.numeric(names(table(sim14$lr0.6))),
       y    = cumsum(table(sim14$lr0.6))/sum(table(sim14$lr0.6)),
       type = "l",
       col  = "red",
       lty  = "dotted")
lines(x   = c(2, 3),
      y   = c(0.9, 0.9),
      col = "red")
text(x      = 3,
     y      = 0.9,
     pos    = 4,
     labels = "red: simulations",
     col    = "red")

## ----em14sim, echo=FALSE, message=FALSE---------------------------------------
library(crossrun)
set.seed(32)
em14s <- simclem(m1=7, nsim = 171600)
em14stab <- round(table(em14s)/100,1)
rownames(em14stab) <- paste0("C=",1:13)
colnames(em14stab) <- paste0("L=",1:7)
knitr::kable(em14stab, caption = 'Table 4')

## ----em60sim, echo=FALSE, message=FALSE---------------------------------------
library(crossrun)
set.seed(32)
em60s <- simclem(m1=30)
em60c <- as.numeric(names(table(em60s$cs)))
em60l <- as.numeric(names(table(em60s$ls)))
minc <- min(em60c)
minl <- min(em60l)
plot(x      = em60c, 
     y      = (cumsum(cumsumm(joint60em)[,30])/(choose(60,30)/2))[em60c],
     type   = "l",
     xlab   = "Number of crossings",
     ylab   = "CDF",
     las    = 1)
points(x    = em60c,
       y    = cumsum(table(em60s$cs))/sum(table(em60s$cs)),
       type = "l",
       col  = "red",
       lty  = "dotted")
lines(x     = c(minc, minc+1),
      y     = c(0.9, 0.9),
      col   = "red")
text(x      = minc+1,
     y      = 0.9,
     pos    = 4,
     labels = "red: simulations",
     col    = "red")
plot(x      = em60l, 
     y      = (cumsum(cumsummcol(joint60em)[59,])/(choose(60,30)/2))[em60l],
     type   = "l",
     xlab   = "Number of crossings",
     ylab   = "CDF",
     las    = 1)
points(x    = em60l,
       y    = cumsum(table(em60s$ls))/sum(table(em60s$ls)),
       type = "l",
       col  = "red",
       lty  = "dotted")
lines(x     = c(minl, minl+1),
      y     = c(0.9, 0.9),
      col   = "red")
text(x      = minl+1,
     y      = 0.9,
     pos    = 4,
     labels = "red: simulations",
     col    = "red")

