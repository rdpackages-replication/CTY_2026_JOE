# rd2d: post analysis of Monte Carlo
# Authors: M. D. Cattaneo, R. Titiunik, R. R. Yu
# Last update: 
# rm(list=ls(all=TRUE))
# library(binsreg); library(ggplot2)

rm(list=ls(all=TRUE))

library(MASS)
library(ggplot2)
library(rdrobust)
library(latex2exp)
library(tidyr)
library(dplyr)
library(haven)
library(xtable)
library(expm)
library(rd2d)

# helper: format 1.000 -> "$\\bb_{1}$", 21.000 -> "$\\bb_{21}$", etc.
fmt_bb <- function(x) {
  s <- format(x, trim = TRUE, scientific = FALSE)           # e.g., "5.000"
  s <- sub("(\\.0+)?$", "", s)                               # -> "5"
  sprintf("$\\bb_{%s}$", s)
}

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

eval <- eval[c(10:30),]
neval <- nrow(eval)

# Polynomial Fits to the Data

formula <- y ~ x.1 + x.2 # DGP 1
# formula <- y ~ x.1 + x.2 + I(x.1^2) + I(x.1 * x.2) + I(x.2^2) # DGP 2

# Fit models separately for d = 0 and d = 1
model_d0 <- lm(formula, data = x[x$d == 0, ])
model_d1 <- lm(formula, data = x[x$d == 1, ])

# Calculate true treatment effects
true_effect <- data.frame(eval) # create copy of eval
true_effect$mu.0 <- predict(model_d0, newdata = true_effect)
true_effect$mu.1 <- predict(model_d1, newdata = true_effect)
true_effect$tau <- (true_effect$mu.1 - true_effect$mu.0)

m <- 2000

# DGP <- 1
# DGP <- 2

# subset <- c(1,5,10,15,21,25,30,35,40)
# 
# neval_subset <- length(subset)
# slope_vec <- c(2.4, 2.8, 3.2, 3.6, 4, 4.4, 4.8)
# sd_vec <- c(0.07, 0.08, 0.09, 0.1)

for (DGP in c(1)){
    
    subset <- c(1:neval)
    neval_subset <- neval
    
    result <- list("kinkoff" = matrix(NA, nrow = neval_subset, ncol = 12),
                   "kinkon" = matrix(NA, nrow = neval_subset, ncol = 12),
                   "adaptive" = matrix(NA, nrow = neval_subset, ncol = 12),
                   "rdrobustadj" = matrix(NA, nrow = neval_subset, ncol = 12))
        
    for (method in c("kinkoff", "kinkon", "adaptive", "rdrobustadj")){
      
      subset <- c(1:neval)
      neval_subset <- neval
      
      h <- rep(0,neval_subset)
      tau.hat <- rep(0,neval_subset)
      ci.lower <- rep(0,neval_subset)
      ci.upper <- rep(0,neval_subset)
      cb.lower <- rep(0,neval_subset)
      cb.upper <- rep(0,neval_subset)
      mse <- rep(0,neval_subset)
      ec.ptw <- rep(0,neval_subset)
      il.ptw <- rep(0,neval_subset)
      ec.unif <- 0
      il.unif <- rep(0,neval_subset)
      
      tau.hat.full <- matrix(0, m, neval_subset)
      tau.hat.q.full <- matrix(0, m, neval_subset)
      
      shift <- 1
      
      for (i in c(1: m)){
        # dat.method <- as.matrix(read.csv(sprintf("monte_outputs/%s_sim%d_DGP%d_slope%d_sd%d.csv", method, i, DGP, slopeid, sdid)))
        # dat.method <- as.matrix(read.csv(sprintf("rd2d_monte/%s_sim%d_m1000_sd01_n20000_hc3.csv", method, i)))
        dat.method <- as.matrix(read.csv(sprintf("monte_outputs/%s_sim%d_DGP%d_new.csv", method, i, DGP)))
        h <- h + dat.method[subset,1]
        tau.hat <- tau.hat + dat.method[subset,1 + shift]
        tau.hat.full[i,] <- dat.method[subset,1 + shift]
        ci.lower <- ci.lower + dat.method[subset,2 + shift]
        ci.upper <- ci.upper + dat.method[subset,3 + shift]
        cb.lower <- cb.lower + dat.method[subset,4 + shift]
        cb.upper <- cb.upper + dat.method[subset,5 + shift]
        tau.hat.q.full[i,] <- (dat.method[subset,5 + shift] - dat.method[subset,4 + shift]) / (dat.method[subset,3 + shift] - dat.method[subset,2 + shift])

        mse <- mse + (dat.method[subset,1 + shift] - true_effect$tau[subset])^2
        ec.ptw <- ec.ptw + as.integer((dat.method[subset,2 + shift] <= true_effect$tau[subset]) & (true_effect$tau[subset] <= dat.method[subset,3 + shift]))
        il.ptw <- il.ptw + dat.method[subset,3 + shift] - dat.method[subset,2 + shift]
        ec.unif <- ec.unif + as.integer(sum((dat.method[subset,4 + shift] <= true_effect$tau[subset]) & (true_effect$tau[subset] <= dat.method[subset,5 + shift])) == neval_subset)
        il.unif <- il.unif + dat.method[subset,5 + shift] - dat.method[subset,4 + shift]
      }
      
      h <- h / m
      tau.hat <- tau.hat / m
      ci.lower <- ci.lower / m
      ci.upper <- ci.upper / m
      cb.lower <- cb.lower / m
      cb.upper <- cb.upper / m
      bias <- tau.hat - true_effect$tau[subset]
      mse <- mse / m
      rmse <- sqrt(mse)
      se <- apply(tau.hat.full, 2, sd)
      ec.ptw <- ec.ptw / m
      il.ptw <- il.ptw / m
      ec.unif <- rep(ec.unif/m, neval_subset)
      il.unif <- il.unif / m
      se.q <- apply(tau.hat.q.full, 2, sd)

      # tau.full.df <- as.data.frame(tau.hat.q)
      # se <- apply(tau.full.df, 2, sd)

      result[[method]]  <- cbind(c(1:neval_subset), h, bias, se, rmse, ec.ptw, il.ptw, ec.unif, il.unif,se.q)

      if (method == "kinkoff"){
        bands <- cbind(ci.lower, ci.upper, cb.lower, cb.upper)
      }

      subset <- c(1:neval)
      # subset <- c(1:neval)
      # subset <- c(1:neval)
      out <- result[[method]]

      # Extract the desired subset and columns
      ec.unif <- out[1,8]
      il.unif <- mean(out[,9])
      
      # # === Build reshaped table ===
      # # out columns: 1=index, 2=h, 3=bias, 4=se, 5=rmse, 6=ec.ptw, 7=il.ptw, (8,9 uniform metrics)
      # idx_vals <- out[subset, 1]
      # h_vals   <- out[subset, 2]
      # stats    <- out[subset, 3:7]  # bias, se, rmse, ec.ptw, il.ptw
      # 
      # fmt_bb_idx <- function(i) sprintf("$\\bb_{%s}$", i)
      # 
      # tab <- data.frame(
      #   " "  = rep("", length(idx_vals)),                # leading blank so rows start with '&'
      #   bb   = vapply(idx_vals, fmt_bb_idx, character(1)),
      #   h    = h_vals,                                   # keep h as numeric column
      #   stats,
      #   check.names = FALSE
      # )
      # 
      # # Digits vector: length must be ncol(tab)+1 (first element for rownames)
      # # Two text cols (blank, bb), then 6 numeric cols → set to 3 decimals
      # digits_vector <- c(0, 0, 0, rep(3, 6))
      # 
      # xt <- xtable(tab, digits = digits_vector)
      # 
      # cat(sprintf("Monte Carlo outputs for %s with DGP %d\n", method, DGP))
      # print(xt,
      #       type = "latex",
      #       include.rownames = FALSE,
      #       include.colnames = FALSE,
      #       sanitize.text.function = identity,
      #       hline.after = NULL)
      # 
      # # Optional end-of-block cline (columns 2–8 correspond to bb through last stat)
      # cat("\\cline{2-8}\n")
      
      # === Build reshaped table ===
      idx_vals <- out[subset, 1]
      h_vals   <- out[subset, 2]
      bias     <- out[subset, 3]
      sd_vals  <- out[subset, 4]
      rmse     <- out[subset, 5]
      ec.ptw   <- out[subset, 6]
      il.ptw   <- out[subset, 7]
      
      fmt_bb_idx <- function(i) sprintf("$\\bb_{%s}$", i)
      
      tab <- data.frame(
        bb   = vapply(idx_vals, fmt_bb_idx, character(1)),
        h    = h_vals,
        Bias = bias,
        SD   = sd_vals,
        RMSE = rmse,
        EC   = ec.ptw,
        IL   = il.ptw,
        check.names = FALSE
      )
      
      # Add final uniform row
      uniform_row <- data.frame(
        bb   = "Uniform",
        h    = NA,
        Bias = NA,
        SD   = NA,
        RMSE = NA,
        EC   = round(ec.unif[1], 4),
        IL   = round(il.unif, 3),
        check.names = FALSE
      )
      tab <- rbind(tab, uniform_row)
      
      # Digits: 1 (rownames) + 7 columns = length 8
      digits_vector <- c(0, 0, rep(3, 6))
      
      # Proper alignment: first entry for rownames (usually 'l')
      xt <- xtable(tab, digits = digits_vector, align = c("l", "c", "c", "r", "r", "r", "r", "r"))
      
      # Print LaTeX table
      cat("\\begin{tabular}{ccrrrrr} \\toprule \\toprule\n")
      cat("$\\bb\\in\\B$ & $h$ & Bias & SD & RMSE & EC & IL\\\\ \\midrule\n")
      
      print(xt,
            type = "latex",
            include.rownames = FALSE,
            include.colnames = FALSE,
            sanitize.text.function = identity,
            hline.after = NULL)
      
      cat("\\cline{1-7}\n")
      cat(sprintf("Uniform &  &  &  &  &  %.4f & %.3f \\\\\n", ec.unif[1], il.unif))
      cat("\\bottomrule\\bottomrule\n\\end{tabular}\n")
      

      # Keep uniform summary prints
      cat(sprintf("Uniform coverage for %s is:\n", method)); print(ec.unif)
      cat(sprintf("Uniform confidence interval length for %s is:\n", method)); print(il.unif)
  }
}

if (DGP == 1){
  result_1 <- result
  bands_1 <- bands
} else {
  result_2 <- result
  bands_2 <- bands
}

##################### Plot: Simulation Estimation and Inference ################

DGP <- 1
# DGP <- 2

if (DGP == 1){
  result <- result_1
  bands <- bands_1
} else {
  result <- result_2
  bands <- bands_2
}

# Four point estimates

tau.true <- (true_effect$mu.1 - true_effect$mu.0) * DGP
tau.kinkoff <- tau.true + result$kinkoff[,3]
tau.kinkon <- tau.true + result$kinkon[,3]
tau.adaptive <- tau.true + result$adaptive[,3]
tau.rdrobustadj <- tau.true + result$rdrobustadj[,3]

# Confidence interval and bands from adaptive

CI.lower <- bands[,1]
CI.upper <- bands[,2]
CB.lower <- bands[,3]
CB.upper <- bands[,4]

indx <- c(1:neval)

# Create a data frame for plotting
indx <- c(1:neval)
df <- data.frame(
  indx = rep(indx, 4),
  y = c(tau.kinkoff, tau.adaptive, tau.kinkon, tau.rdrobustadj),
  label = rep(c("Smooth","Adaptive", "Unknown Kink", "Rdrobust"), each = length(indx))
)

df.true <- data.frame(
  indx = indx,
  y = tau.true,
  label = rep(c("Population"), each = length(indx))
)

temp_plot <- ggplot() + theme_bw()

# scatter plot of point estimations
temp_plot <- temp_plot + geom_point(data = df, aes(x = indx, y = y, color = label, shape = label, fill = label, linetype = label))

# line plot of true effects
temp_plot <- temp_plot + geom_line(data = df.true,
                                   aes(x = indx, y = y, color = label, shape = label, fill = label, linetype = label),
                                   size = 0.5, show.legend = TRUE)

temp_plot <- temp_plot + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15,face = "bold"),  # X-axis label size
    axis.title.y = element_text(size = 15,face = "bold"),  # Y-axis label size
    plot.title = element_text(size = 20, hjust = 0.5),  # Title size and centering
    text=element_text(family="Times New Roman", face="bold"),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = c(0.76, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),  # White background, no border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

legend_order <- c(
  "Smooth",
  "Adaptive",
  "Unknown Kink",
  "Rdrobust",
  "Population",
  "Pointwise CI",
  "Uniform CB"
)

temp_plot <- temp_plot + 
  scale_color_manual(
    values = c(
      "Smooth" = "blueviolet",
      "Adaptive" = "grey",
      "Unknown Kink" = "blue",
      "Rdrobust" = "darkgreen",
      "Population" = "black",
      "Pointwise CI" = "burlywood4",
      "Uniform CB" = "dodgerblue4"
    ),
    name = NULL,
    breaks = legend_order
  ) +
  scale_shape_manual(
    values = c(
      "Smooth" = 16,         # filled circle
      "Adaptive" = 15,      # tick/cross
      "Unknown Kink" = 4,   # no shape
      "Rdrobust" = 17,       # filled triangle
      "Population" = 15,     # filled square
      "Pointwise CI" = NA,   # no shape
      "Uniform CB" = 0       # open square
    ),
    name = NULL,
    breaks = legend_order
  ) + 
  scale_fill_manual(
    values = c(
      "Smooth" = NA,
      "Adaptive" = NA,
      "Unknown Kink" = NA,
      "Rdrobust" = NA,
      "Population" = NA,
      "Pointwise CI" = NA,
      "Uniform CB" = "dodgerblue4"
    ),
    name = NULL,
    breaks = legend_order
  ) +
  scale_linetype_manual(
    values = c(
      "Smooth" = 1,         # solid
      "Adaptive" = 5,       # longdash
      "Unknown Kink" = 2,   # dashed
      "Rdrobust" = 10,       # dotted
      "Population" = 1,     # solid
      "Pointwise CI" = 6,   # twodash
      "Uniform CB" = 0      # blank
    ),
    name = NULL,
    breaks = legend_order
  ) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 1),
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 1)
  )

temp_plot <- temp_plot + xlab("Cutoffs on the Boundary") + ylab("Treatment Effect")

print(temp_plot)

ggsave("Results/bias.png", temp_plot, width = 6, height = 5)


