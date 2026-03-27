################################################################################################
# Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods
# Simulations
# Authors: M. D. Cattaneo, R. Titiunik, R. R. Yu
################################################################################################

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

#################### Polynomial / Variance Fits to the Data ####################

# Mean models
formula_linear <- y ~ x.1 + x.2
formula_quad   <- y ~ x.1 + x.2 + I(x.1^2) + I(x.2^2) + I(x.1 * x.2)

mean_models <- list(
  linear = list(
    d0 = lm(formula_linear, data = x[x$d == 0, ]),
    d1 = lm(formula_linear, data = x[x$d == 1, ])
  ),
  quadratic = list(
    d0 = lm(formula_quad, data = x[x$d == 0, ]),
    d1 = lm(formula_quad, data = x[x$d == 1, ])
  )
)

# Homoskedastic variance: constant SD on each side
sd_homo_d0 <- sd(resid(mean_models$linear$d0))
sd_homo_d1 <- sd(resid(mean_models$linear$d1))

# Heteroskedastic variance model: log Var(Y|X) = a + b1 x1 + b2 x2 + c1 x1^2 + c2 x2^2 + c12 x1 x2
fit_var_model <- function(dat, mean_fit) {
  rr <- resid(mean_fit)
  dat2 <- dat
  dat2$log_r2 <- log(pmax(rr^2, 1e-8))
  lm(log_r2 ~ x.1 + x.2 + I(x.1^2) + I(x.2^2) + I(x.1 * x.2), data = dat2)
}

# Better to pair variance fit with the richer quadratic mean fit
var_models <- list(
  d0 = fit_var_model(x[x$d == 0, ], mean_models$quadratic$d0),
  d1 = fit_var_model(x[x$d == 1, ], mean_models$quadratic$d1)
)

predict_mu <- function(newdata, d_vec, mean_type) {
  mu0 <- predict(mean_models[[mean_type]]$d0, newdata = newdata)
  mu1 <- predict(mean_models[[mean_type]]$d1, newdata = newdata)
  ifelse(d_vec == 0, mu0, mu1)
}

predict_sigma <- function(newdata, d_vec, var_type) {
  if (var_type == "homoskedastic") {
    ifelse(d_vec == 0, sd_homo_d0, sd_homo_d1)
  } else {
    logv0 <- predict(var_models$d0, newdata = newdata)
    logv1 <- predict(var_models$d1, newdata = newdata)
    sd0 <- sqrt(pmax(exp(logv0), 1e-12))
    sd1 <- sqrt(pmax(exp(logv1), 1e-12))
    ifelse(d_vec == 0, sd0, sd1)
  }
}

dgp_grid <- expand.grid(
  mean_type = c("linear", "quadratic"),
  var_type  = c("homoskedastic", "heteroskedastic"),
  stringsAsFactors = FALSE
)


# write into tables

dir.create("Results", showWarnings = FALSE)

fmt <- function(x) {
  if (is.na(x) || abs(x) < 1e-12) return("$0$")
  s <- formatC(x, format = "e", digits = 2)
  p <- strsplit(s, "e")[[1]]
  a <- p[1]
  b <- as.integer(p[2])
  if (b == 0) {
    paste0("$", a, "$")
  } else {
    paste0("$", a, " \\times 10^{", b, "}$")
  }
}

get <- function(m, name) {
  if (is.null(m)) return("$0$")
  if (name %in% names(coef(m))) fmt(coef(m)[name]) else "$0$"
}

rows <- c("(Intercept)", "x.1", "x.2", "I(x.1^2)", "I(x.2^2)", "I(x.1 * x.2)")

labels_mean <- c("$\\beta_{t,0}$","$\\beta_{t,11}$","$\\beta_{t,12}$",
                 "$\\beta_{t,21}$","$\\beta_{t,22}$","$\\beta_{t,23}$")

labels_var <- c("$\\gamma_{t,0}$","$\\gamma_{t,11}$","$\\gamma_{t,12}$",
                "$\\gamma_{t,21}$","$\\gamma_{t,22}$","$\\gamma_{t,23}$")

# ---------- mean table ----------
lines <- c(
  "\\begin{tabular}{@{}rrrrr@{}}",
  "\\hline\\hline",
  "    & \\multicolumn{2}{c}{Linear} & \\multicolumn{2}{c}{Quadratic} \\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "    & \\multicolumn{1}{c}{$t=0$} & \\multicolumn{1}{c}{$t=1$} & \\multicolumn{1}{c}{$t=0$} & \\multicolumn{1}{c}{$t=1$} \\\\ \\midrule"
)

for (i in seq_along(rows)) {
  r <- rows[i]
  vals <- c(
    get(mean_models$linear$d0, r),
    get(mean_models$linear$d1, r),
    get(mean_models$quadratic$d0, r),
    get(mean_models$quadratic$d1, r)
  )
  lines <- c(lines, paste0(labels_mean[i], " & ", paste(vals, collapse = " & "), " \\\\"))
}

lines <- c(lines, "\\hline\\hline", "\\end{tabular}")
writeLines(lines, "Results/dgp-mean.tex")


# ---------- variance table ----------
lines <- c(
  "\\begin{tabular}{@{}rrrrr@{}}",
  "\\hline\\hline",
  "    & \\multicolumn{2}{c}{Homoskedastic} & \\multicolumn{2}{c}{Heteroskedastic} \\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "    & \\multicolumn{1}{c}{$t=0$} & \\multicolumn{1}{c}{$t=1$} & \\multicolumn{1}{c}{$t=0$} & \\multicolumn{1}{c}{$t=1$} \\\\ \\midrule"
)

for (i in seq_along(rows)) {
  r <- rows[i]
  
  if (r == "(Intercept)") {
    homo <- c(fmt(log(sd_homo_d0^2)), fmt(log(sd_homo_d1^2)))
  } else {
    homo <- c("$0$", "$0$")
  }
  
  hetero <- c(
    get(var_models$d0, r),
    get(var_models$d1, r)
  )
  
  vals <- c(homo, hetero)
  
  lines <- c(lines, paste0(labels_var[i], " & ", paste(vals, collapse = " & "), " \\\\"))
}

lines <- c(lines, "\\hline\\hline", "\\end{tabular}")
writeLines(lines, "Results/dgp-var.tex")

m <- 2000

alpha_1 <- 3
beta_1 <- 4
alpha_2 <- 3
beta_2 <- 4

n <- 20 * 20 * 25 * 2

#################### Polynomial / Variance Fits to the Data ####################

num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores - 3)
registerDoParallel(cl)  # Register parallel backend

system.time({
  foreach(
    j = 1:m,
    .packages = c("MASS","tidyr","lmtest","sandwich","ggplot2","rdrobust",
                  "latex2exp","sf","dplyr","progressr","doParallel",
                  "foreach","expm","rd2d"),
    .export = c("eval", "adaptive_bandwidth", "dist_to_kink",
                "alpha_1", "beta_1", "alpha_2", "beta_2", "n",
                "mean_models", "var_models",
                "sd_homo_d0", "sd_homo_d1",
                "predict_mu", "predict_sigma", "dgp_grid")
  ) %dopar% {
    
    eval.subset <- eval
    neval <- nrow(eval)
    
    for (g in 1:nrow(dgp_grid)) {
      
      mean_type_now <- dgp_grid$mean_type[g]
      var_type_now  <- dgp_grid$var_type[g]
      DGP <- g
      
      # generate new covariates only
      x.1 <- 100 * rbeta(n, alpha_1, beta_1) - 25
      x.2 <- 100 * rbeta(n, alpha_2, beta_2) - 25
      X <- data.frame(x.1 = x.1, x.2 = x.2)
      
      d <- as.numeric(x.1 >= 0 & x.2 >= 0)
      
      # generate new outcomes from the already-fitted models
      mu    <- predict_mu(X, d, mean_type_now)
      sigma <- predict_sigma(X, d, var_type_now)
      y     <- mu + rnorm(n, mean = 0, sd = sigma)
      
      D <- proxy::dist(X, eval.subset, method = "euclidean")
      d_expanded <- matrix(rep(2 * d - 1, times = ncol(D)),
                           nrow = nrow(D), ncol = ncol(D))
      D <- D * d_expanded
      
      # smooth boundary
      result.kinkoff <- rd2d.dist(y, D, kink = "off", repp = 5000, vce = "hc1")
      out.kinkoff <- cbind(result.kinkoff$results$h0, result.kinkoff$results$Est.p,
                           result.kinkoff$results$CI.lower, result.kinkoff$results$CI.upper,
                           result.kinkoff$results$CB.lower, result.kinkoff$results$CB.upper,
                           result.kinkoff$results$Est.q)
      
      # unknown kink location
      result.kinkon <- rd2d.dist(y, D, kink = "on", repp = 5000, vce = "hc1", rbc = "off")
      out.kinkon <- cbind(result.kinkon$results$h0, result.kinkon$results$Est.p,
                          result.kinkon$results$CI.lower, result.kinkon$results$CI.upper,
                          result.kinkon$results$CB.lower, result.kinkon$results$CB.upper,
                          result.kinkon$results$Est.q)
      
      # adaptive
      h.adaptive <- adaptive_bandwidth(result.kinkoff$results$h0,
                                       result.kinkon$results$h0,
                                       dist_to_kink)
      h.adaptive <- cbind(h.adaptive, h.adaptive)
      result.adaptive <- rd2d.dist(y, D, kink = "off", repp = 5000, vce = "hc1", h = h.adaptive)
      out.adaptive <- cbind(result.adaptive$results$h0, result.adaptive$results$Est.p,
                            result.adaptive$results$CI.lower, result.adaptive$results$CI.upper,
                            result.adaptive$results$CB.lower, result.adaptive$results$CB.upper,
                            result.adaptive$results$Est.q)
      
      # rdrobust
      bws <- matrix(0, nrow = neval, ncol = 2)
      for (i in 1:neval) {
        out <- rdbwselect(y, D[, i], vce = "hc1")
        bws[i, 1] <- out$bws[1]
        bws[i, 2] <- out$bws[2]
      }
      result.rdrobust <- rd2d.dist(y, D, h = bws, kink = "off", repp = 5000, vce = "hc1")
      out.rdrobust <- cbind(result.rdrobust$results$h0, result.rdrobust$results$Est.p,
                            result.rdrobust$results$CI.lower, result.rdrobust$results$CI.upper,
                            result.rdrobust$results$CB.lower, result.rdrobust$results$CB.upper,
                            result.rdrobust$results$Est.q)
      
      tag <- paste0(mean_type_now, "_", var_type_now)
      
      write.csv(out.kinkoff,
                file = sprintf("monte_variants_outputs/kinkoff_sim%d_DGP%d_%s.csv", j, DGP, tag),
                row.names = FALSE)
      write.csv(out.kinkon,
                file = sprintf("monte_variants_outputs/kinkon_sim%d_DGP%d_%s.csv", j, DGP, tag),
                row.names = FALSE)
      write.csv(out.adaptive,
                file = sprintf("monte_variants_outputs/adaptive_sim%d_DGP%d_%s.csv", j, DGP, tag),
                row.names = FALSE)
      write.csv(out.rdrobust,
                file = sprintf("monte_variants_outputs/rdrobustadj_sim%d_DGP%d_%s.csv", j, DGP, tag),
                row.names = FALSE)
    }
  }
})

on.exit(stopCluster(cl), add = TRUE)