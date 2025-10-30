rm(list=ls(all=TRUE))

library(MASS)
library(tidyr)
library("lmtest")
library("sandwich")
library(ggplot2)
library(rdrobust)
library(latex2exp)
library(sf)
library(dplyr)  # For sampling
library(progressr)
library(doParallel)
library(foreach)
library(expm)
library(rdrobust)
library(rd2d)

set.seed(3)

################################## Load Data ###################################

x <- read.csv("spp.csv")
x$X <- NULL
colnames(x) <- c("x.1", "x.2","y","d")
na.ok <- complete.cases(x$x.1) & complete.cases(x$x.2)
x <- x[na.ok,]

neval <- 40
eval <- matrix(nrow = neval, ncol = 2)
for (i in 1: ceiling(neval * 0.5)){
  eval[i,] <- c(0, 40 - (i-1) * 40 / ceiling(neval * 0.5))
}
for (i in (ceiling(neval * 0.5)+1): neval){
  eval[i,] <- c((i - ceiling(neval * 0.5) - 1) *56 / (ceiling(neval * 0.5)),0)
}
eval <- data.frame(eval)
colnames(eval) <- c("x.1", "x.2")

# subsetting eval
eval <- eval[c(10:30),]
neval <- 21

dist_to_kink <- rep(0, neval)
for (i in 1:neval){
  dist_to_kink[i] <- min(sqrt(eval[i,1]^2 + eval[i,2]^2), 75 - eval[i,1], 75 - eval[i,2])
}

########################### Kink Adaptive Bandwidth ############################

adaptive_bandwidth <- function(h_smooth, h_kink, dist_to_kink){
  out <- pmax(pmin(h_smooth, dist_to_kink), h_kink)
  return(out)
}

######################## Polynomial Fits to the Data ###########################

DGP <- 1
# DGP <- 2

formula <- y ~ x.1 + x.2

# Fit models separately for d = 0 and d = 1
model_d0 <- lm(formula, data = x[x$d == 0, ])
model_d1 <- lm(formula, data = x[x$d == 1, ])

# Summary of the models
summary(model_d0)
summary(model_d1)

# Simulate data with independent Beta-distributed coordinates
m <- 2000

alpha_1 <- 3
beta_1 <- 4
alpha_2 <- 3
beta_2 <- 4

n <- 20 * 20 * 25 * 2
# n <- 40 * 40 * 25
# n <- 10 * 10 * 25

num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores - 3)
registerDoParallel(cl)  # Register parallel backend

system.time({
  foreach(j = 1:m) %dopar% {
    library(MASS)
    library(tidyr)
    library("lmtest")
    library("sandwich")
    library(ggplot2)
    library(rdrobust)
    library(latex2exp)
    library(sf)
    library(dplyr)  # For sampling
    library(progressr)
    library(doParallel)
    library(foreach)
    library(expm)
    library(rdrobust)
    library(rd2d)
    
    for (DGP in c(1)){
      eval.subset <- eval
      neval <- nrow(eval)
      
      # generate data set for the ith monte carlo
      
      x.1 <- 100 * rbeta(n, alpha_1, beta_1) - 25
      x.2 <- 100 * rbeta(n, alpha_2, beta_2) - 25
      
      X <- cbind(x.1, x.2)
      X <- as.data.frame(X)
      d <- x.1 >= 0 & x.2 >= 0
      y <- ifelse(d == 0, predict(model_d0, newdata = X), predict(model_d1, newdata = X)) 
      # y <- ifelse(d == 0, predict(model_d0, newdata = X), predict(model_d1, newdata = X)) * 3 # DGP == 2
      if (DGP == 1){
        sd0 <- 0.3321
        sd1 <- 0.4352
      } else {
        sd0 <- 0.1
        sd1 <- 0.1
      }
      y <- y + rnorm(n,mean = 0, sd = sd0) * (1 - d) + rnorm(n,mean = 0, sd = sd1) * d
      
      D <- proxy::dist(X, eval.subset, method = "euclidean")  # Use "euclidean" for Euclidean distances
      d_expanded <- matrix(rep(2 * d - 1, times = ncol(D)), nrow = nrow(D), ncol = ncol(D))
      D <- D * d_expanded
      
      # smooth boundary
      result.kinkoff <- rd2d.dist(y,D,kink = "off", repp = 5000, vce = "hc1")
      out.kinkoff <- cbind(result.kinkoff$results$h0, result.kinkoff$results$Est.p, result.kinkoff$results$CI.lower,
                           result.kinkoff$results$CI.upper, result.kinkoff$results$CB.lower,
                           result.kinkoff$results$CB.upper, result.kinkoff$results$Est.q)
      
      # unknown kink location
      result.kinkon <- rd2d.dist(y,D,kink = "on", repp = 5000, vce = "hc1", rbc = "off")
      out.kinkon <- cbind(result.kinkon$results$h0, result.kinkon$results$Est.p, result.kinkon$results$CI.lower,
                          result.kinkon$results$CI.upper, result.kinkon$results$CB.lower,
                          result.kinkon$results$CB.upper, result.kinkon$results$Est.q)
      
      # kink adaptive
      h.adaptive <- adaptive_bandwidth(result.kinkoff$results$h0, result.kinkon$results$h0, dist_to_kink)
      h.adaptive <- cbind(h.adaptive, h.adaptive)
      result.adaptive <- rd2d.dist(y,D,kink = "off", repp = 5000, vce = "hc1", h = h.adaptive)
      out.adaptive <- cbind(result.adaptive$results$h0, result.adaptive$results$Est.p, result.adaptive$results$CI.lower,
                            result.adaptive$results$CI.upper, result.adaptive$results$CB.lower,
                            result.adaptive$results$CB.upper, result.adaptive$results$Est.q)
      
      # rdrobust
      bws <- matrix(0, nrow = neval, ncol = 2)
      for (i in 1:neval){
        out <- rdbwselect(y, D[,i], vce = "hc1")
        bws[i,1] <- out$bws[1]
        bws[i,2] <- out$bws[2]
      }
      result.rdrobust <-  rd2d.dist(y,D, h = bws, kink = "off", repp = 5000, vce = "hc1")
      out.rdrobust <- cbind(result.rdrobust$results$h0, result.rdrobust$results$Est.p, result.rdrobust$results$CI.lower,
                            result.rdrobust$results$CI.upper, result.rdrobust$results$CB.lower,
                            result.rdrobust$results$CB.upper, result.rdrobust$results$Est.q)
      
      write.csv(out.kinkoff, file = sprintf("monte_outputs/kinkoff_sim%d_DGP%d_new.csv", j, DGP), row.names = FALSE)
      write.csv(out.kinkon, file = sprintf("monte_outputs/kinkon_sim%d_DGP%d_new.csv", j, DGP), row.names = FALSE)
      write.csv(out.adaptive, file = sprintf("monte_outputs/adaptive_sim%d_DGP%d_new.csv", j, DGP), row.names = FALSE)
      write.csv(out.rdrobust, file = sprintf("monte_outputs/rdrobustadj_sim%d_DGP%d_new.csv", j, DGP), row.names = FALSE)
    }
  }})

stopCluster(cl)
