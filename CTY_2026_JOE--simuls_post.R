################################################################################################
# Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods
# Simulations -- Tables
# Authors: M. D. Cattaneo, R. Titiunik, R. R. Yu
################################################################################################

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
subset <- c(1:neval)
neval_subset <- neval

############################
# Simulation settings
############################
m <- 2000

methods <- c("kinkoff", "kinkon", "adaptive", "rdrobustadj")

############################
# Four DGPs
############################
dgp_grid <- data.frame(
  mean_type = c("linear", "quadratic", "linear", "quadratic"),
  var_type  = c("homoskedastic", "homoskedastic", "heteroskedastic", "heteroskedastic"),
  stringsAsFactors = FALSE
)

dgp_grid$tag <- paste0(dgp_grid$mean_type, "_", dgp_grid$var_type)

dgp_grid$label <- c(
  "Linear + Homoskedastic",
  "Quadratic + Homoskedastic",
  "Linear + Heteroskedastic",
  "Quadratic + Heteroskedastic"
)
############################
# Fit mean models once
############################
formula_linear <- y ~ x.1 + x.2
formula_quad   <- y ~ x.1 + x.2 + I(x.1^2) + I(x.1 * x.2) + I(x.2^2)

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

############################
# True treatment effect for each DGP
# Note: true effect depends only on mean model
############################
true_effect_list <- list()

for (g in 1:nrow(dgp_grid)) {
  mean_type_now <- dgp_grid$mean_type[g]
  tag_now <- dgp_grid$tag[g]
  
  true_effect <- data.frame(eval)
  true_effect$mu.0 <- predict(mean_models[[mean_type_now]]$d0, newdata = true_effect)
  true_effect$mu.1 <- predict(mean_models[[mean_type_now]]$d1, newdata = true_effect)
  true_effect$tau  <- true_effect$mu.1 - true_effect$mu.0
  
  true_effect_list[[tag_now]] <- true_effect
}

############################
# Containers
############################
results_all <- list()
bands_all   <- list()

############################
# Monte Carlo summaries
############################
for (g in 1:nrow(dgp_grid)) {
  
  DGP <- g
  tag <- dgp_grid$tag[g]
  label_now <- dgp_grid$label[g]
  true_effect <- true_effect_list[[tag]]
  
  result <- list(
    kinkoff     = matrix(NA, nrow = neval_subset, ncol = 10),
    kinkon      = matrix(NA, nrow = neval_subset, ncol = 10),
    adaptive    = matrix(NA, nrow = neval_subset, ncol = 10),
    rdrobustadj = matrix(NA, nrow = neval_subset, ncol = 10)
  )
  
  bands <- NULL
  
  for (method in methods) {
    
    h        <- rep(0, neval_subset)
    tau.hat  <- rep(0, neval_subset)
    ci.lower <- rep(0, neval_subset)
    ci.upper <- rep(0, neval_subset)
    cb.lower <- rep(0, neval_subset)
    cb.upper <- rep(0, neval_subset)
    mse      <- rep(0, neval_subset)
    ec.ptw   <- rep(0, neval_subset)
    il.ptw   <- rep(0, neval_subset)
    ec.unif  <- 0
    il.unif  <- rep(0, neval_subset)
    
    tau.hat.full   <- matrix(0, nrow = m, ncol = neval_subset)
    tau.hat.q.full <- matrix(0, nrow = m, ncol = neval_subset)
    
    shift <- 1
    
    for (i in 1:m) {
      
      file_now <- sprintf(
        "simuls/%s_sim%d_DGP%d_%s.csv",
        method, i, DGP, tag
      )
      
      if (!file.exists(file_now)) {
        stop(sprintf("Missing file: %s", file_now))
      }
      
      dat.method <- as.matrix(read.csv(file_now))
      
      h <- h + dat.method[subset, 1]
      tau.hat <- tau.hat + dat.method[subset, 1 + shift]
      tau.hat.full[i, ] <- dat.method[subset, 1 + shift]
      ci.lower <- ci.lower + dat.method[subset, 2 + shift]
      ci.upper <- ci.upper + dat.method[subset, 3 + shift]
      cb.lower <- cb.lower + dat.method[subset, 4 + shift]
      cb.upper <- cb.upper + dat.method[subset, 5 + shift]
      
      tau.hat.q.full[i, ] <- (dat.method[subset, 5 + shift] - dat.method[subset, 4 + shift]) /
        (dat.method[subset, 3 + shift] - dat.method[subset, 2 + shift])
      
      mse <- mse + (dat.method[subset, 1 + shift] - true_effect$tau[subset])^2
      
      ec.ptw <- ec.ptw + as.integer(
        (dat.method[subset, 2 + shift] <= true_effect$tau[subset]) &
          (true_effect$tau[subset] <= dat.method[subset, 3 + shift])
      )
      
      il.ptw <- il.ptw + dat.method[subset, 3 + shift] - dat.method[subset, 2 + shift]
      
      ec.unif <- ec.unif + as.integer(
        all(
          (dat.method[subset, 4 + shift] <= true_effect$tau[subset]) &
            (true_effect$tau[subset] <= dat.method[subset, 5 + shift])
        )
      )
      
      il.unif <- il.unif + dat.method[subset, 5 + shift] - dat.method[subset, 4 + shift]
    }
    
    h        <- h / m
    tau.hat  <- tau.hat / m
    ci.lower <- ci.lower / m
    ci.upper <- ci.upper / m
    cb.lower <- cb.lower / m
    cb.upper <- cb.upper / m
    
    bias   <- tau.hat - true_effect$tau[subset]
    mse    <- mse / m
    rmse   <- sqrt(mse)
    se     <- apply(tau.hat.full, 2, sd)
    ec.ptw <- ec.ptw / m
    il.ptw <- il.ptw / m
    ec.unif.scalar <- ec.unif / m
    il.unif <- il.unif / m
    se.q   <- apply(tau.hat.q.full, 2, sd)
    
    result[[method]] <- cbind(
      Index = 1:neval_subset,
      h = h,
      Bias = bias,
      SE = se,
      RMSE = rmse,
      EC.ptw = ec.ptw,
      IL.ptw = il.ptw,
      EC.unif = rep(ec.unif.scalar, neval_subset),
      IL.unif = il.unif,
      SE.q = se.q
    )
    
    if (method == "kinkoff") {
      bands <- cbind(
        CI.lower = ci.lower,
        CI.upper = ci.upper,
        CB.lower = cb.lower,
        CB.upper = cb.upper
      )
    }
    
    cat("\n")
    cat(sprintf("DGP %d: %s\n", DGP, label_now))
    cat(sprintf("Method: %s\n", method))
    cat(sprintf("Uniform coverage: %.4f\n", ec.unif.scalar))
    cat(sprintf("Average uniform band length: %.4f\n", mean(il.unif)))
  }
  
  results_all[[tag]] <- result
  bands_all[[tag]]   <- bands
}

########################################################
# Save 16 LaTeX tabulars: one for each DGP × method
########################################################

dir.create("tables", showWarnings = FALSE)

methods <- c("kinkoff", "kinkon", "adaptive", "rdrobustadj")
method_names <- c(
  kinkoff     = "hsmooth",
  adaptive    = "hkink",
  kinkon      = "unknown kink",
  rdrobustadj = "rdrobust"
)

for (g in 1:nrow(dgp_grid)) {
  
  tag <- dgp_grid$tag[g]
  
  for (method in methods) {
    
    out <- as.data.frame(results_all[[tag]][[method]])
    
    tab <- data.frame(
      bb   = sprintf("$\\bb_{%s}$", out$Index),
      h    = out$h,
      Bias = out$Bias,
      SD   = out$SE,
      RMSE = out$RMSE,
      EC   = out$EC.ptw,
      IL   = out$IL.ptw,
      check.names = FALSE
    )
    
    xt <- xtable(
      tab,
      digits = c(0, 0, rep(3, 6)),
      align = c("l", "c", "c", "r", "r", "r", "r", "r")
    )
    
    lines <- c(
      "\\begin{tabular}{ccrrrrr}",
      "\\toprule \\toprule",
      "\\multicolumn{1}{c}{$\\bb\\in\\B$} & \\multicolumn{1}{c}{$h$} & \\multicolumn{1}{c}{Bias} & \\multicolumn{1}{c}{SD} & \\multicolumn{1}{c}{RMSE} & \\multicolumn{1}{c}{EC} & \\multicolumn{1}{c}{IL} \\\\",
      "\\midrule"
    )
    
    body <- print(
      xt,
      type = "latex",
      include.rownames = FALSE,
      include.colnames = FALSE,
      sanitize.text.function = identity,
      hline.after = NULL,
      comment = FALSE,
      only.contents = TRUE,
      print.results = FALSE
    )
    
    ec_unif <- out$EC.unif[1]
    il_unif <- mean(out$IL.unif, na.rm = TRUE)
    
    lines <- c(
      lines,
      body,
      "\\cline{1-7}",
      sprintf("Uniform &  &  &  &  & %.3f & %.3f \\\\", ec_unif, il_unif),
      "\\bottomrule \\bottomrule",
      "\\end{tabular}"
    )
    
    file_name <- sprintf("tables/tab-simuls_%s_%s.tex", tag, method)
    writeLines(lines, file_name)
  }
}

##################### Plot: Bias Figure for Each DGP #####################

dir.create("figures", showWarnings = FALSE)

methods_plot <- c("kinkoff", "adaptive", "kinkon", "rdrobustadj")

method_labels <- c(
  kinkoff     = "Smooth",
  adaptive    = "Adaptive",
  kinkon      = "Unknown Kink",
  rdrobustadj = "Rdrobust"
)

for (g in 1:nrow(dgp_grid)) {
  
  tag <- dgp_grid$tag[g]
  label_now <- dgp_grid$label[g]
  
  out_kinkoff     <- as.data.frame(results_all[[tag]][["kinkoff"]])
  out_adaptive    <- as.data.frame(results_all[[tag]][["adaptive"]])
  out_kinkon      <- as.data.frame(results_all[[tag]][["kinkon"]])
  out_rdrobustadj <- as.data.frame(results_all[[tag]][["rdrobustadj"]])
  true_effect     <- true_effect_list[[tag]]
  
  indx <- out_kinkoff$Index
  
  df <- rbind(
    data.frame(
      indx  = indx,
      y     = true_effect$tau,
      label = "Population"
    ),
    data.frame(
      indx  = indx,
      y     = out_kinkoff$Bias + true_effect$tau,
      label = method_labels["kinkoff"]
    ),
    data.frame(
      indx  = indx,
      y     = out_adaptive$Bias + true_effect$tau,
      label = method_labels["adaptive"]
    ),
    data.frame(
      indx  = indx,
      y     = out_kinkon$Bias + true_effect$tau,
      label = method_labels["kinkon"]
    ),
    data.frame(
      indx  = indx,
      y     = out_rdrobustadj$Bias + true_effect$tau,
      label = method_labels["rdrobustadj"]
    )
  )
  
  legend_order <- c(
    "Smooth",
    "Adaptive",
    "Unknown Kink",
    "Rdrobust",
    "Population"
  )
  
  temp_plot <- ggplot(
    df,
    aes(x = indx, y = y, color = label, shape = label, linetype = label)
  ) +
    geom_point(size = 1.6) +
    geom_line(linewidth = 0.5) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      plot.title   = element_text(size = 20, hjust = 0.5),
      text         = element_text(family = "Times New Roman", face = "bold"),
      axis.text.x  = element_text(face = "bold", size = 15),
      axis.text.y  = element_text(face = "bold", size = 12),
      legend.position = c(0.76, 1),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_color_manual(
      values = c(
        "Smooth"       = "blueviolet",
        "Adaptive"     = "grey",
        "Unknown Kink" = "blue",
        "Rdrobust"     = "darkgreen",
        "Population"   = "black"
      ),
      name = NULL,
      breaks = legend_order
    ) +
    scale_shape_manual(
      values = c(
        "Smooth"       = 16,
        "Adaptive"     = 15,
        "Unknown Kink" = 4,
        "Rdrobust"     = 17,
        "Population"   = NA
      ),
      name = NULL,
      breaks = legend_order
    ) +
    scale_linetype_manual(
      values = c(
        "Smooth"       = 1,
        "Adaptive"     = 5,
        "Unknown Kink" = 2,
        "Rdrobust"     = 10,
        "Population"   = 1
      ),
      name = NULL,
      breaks = legend_order
    ) +
    guides(
      color    = guide_legend(order = 1),
      shape    = guide_legend(order = 1),
      linetype = guide_legend(order = 1)
    ) +
    coord_cartesian(ylim = c(0.31, 0.41)) +
    xlab("Cutoffs on the Boundary") +
    ylab("Treatment Effect") 
  
  print(temp_plot)
  
  ggsave(
    filename = sprintf("figures/fig3-bias_%s.png", tag),
    plot = temp_plot,
    width = 6,
    height = 5
  )
}
