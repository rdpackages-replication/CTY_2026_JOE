################################################################################################
# Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods
# Empirical Application
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
library(rdrobust)
library(grid)
library(rd2d)

################################## Load Data ###################################

data <- read.csv("spp.csv")
data$X <- NULL
colnames(data) <- c("x.1", "x.2","y","d")
na.ok <- complete.cases(data$x.1) & complete.cases(data$x.2)
data <- data[na.ok,]

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

# subsetting
eval <- eval[c(10:30),]
neval <- 21

dist_to_kink <- rep(0, neval)
for (i in 1:neval){
  dist_to_kink[i] <- sqrt(eval[i,1]^2 + eval[i,2]^2)
}

########################### Kink Adaptive Bandwidth ############################

adaptive_bandwidth <- function(h_smooth, h_kink, dist_to_kink){
  out <- pmax(pmin(h_smooth, dist_to_kink), h_kink)
  return(out)
}

################################################################################
### FIGURE 1: Scatter plot score and RDPLOTs
################################################################################

## Scatter Plot and Boundary

eval_labeled <- eval %>%
  mutate(
    index = row_number(),
    label = as.list(if_else(
      index %in% c(4, 7, 10, 15, 18),
      paste0("$\\textbf{b}_{", index, "}$"),  # e.g. "$x_{10}$"
      NA_character_
    ))
  )


bound <- 21
point_1 <- eval[1,]
x_1 <- point_1$x.1
y_1 <- point_1$x.2

# point_12 <- eval[12,]
# x_12 <- point_12$x.1
# y_12 <- point_12$x.2

point_12 <- eval[12,]
x_12 <- point_12$x.1
y_12 <- point_12$x.2

point_21 <- eval[21,]
x_21 <- point_21$x.1
y_21 <- point_21$x.2

p1.1 <- ggplot() +
  geom_point(
    data = data  %>% sample_frac(0.3),
    aes(x = x.1, y = x.2, color = factor(d)),
    alpha = 0.5, size = 0.5
  ) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 60),
               linetype = "solid", color = "grey", alpha = 1, size = 3) +
  geom_segment(aes(x = 0, xend = 100, y = 0, yend = 0),
               linetype = "solid", color = "grey", alpha = 1, size = 3) +
  
  geom_point(
    data = eval_labeled,
    aes(x = x.1, y = x.2),  # Map color to factor(d)
    alpha = 1,
    size = 0.5
  ) +
  
  # Use scale_color_manual to define your own colors for d=0 and d=1
  scale_color_manual(
    name = NULL,  # Legend title (optional)
    values = c("0" = "#619CFF",  # Color for d=0
               "1" = "#F8766D"), # Color for d=1
    labels = c("0" = "Control", "1" = "Treatment")
  )

annotation_data <- eval_labeled %>%
  filter(!is.na(label))


for(i in seq_len(nrow(annotation_data))) {
  if (i <= 3){
    hjust <- -0.3
    vjust <- 0
  } else {
    hjust <- 0.2
    vjust <- 1.5
  }
  p1.1 <- p1.1 + annotate("text",
                          x = annotation_data$x.1[i],
                          y = annotation_data$x.2[i],
                          label = TeX(as.character(annotation_data$label[i])),
                          hjust = hjust,
                          vjust = vjust,
                          size = 4,
                          color = "black",
                          fontface = "bold")}

# Arrow pointing to point 21
p1.1 <- p1.1 +  geom_segment(
  aes(x = x_12 - 8, y = y_12 - 8, xend = x_12, yend = y_12),
  arrow = arrow(length = unit(0.2, "cm")),
  color = "black",
  size = 1
) +
  
  # Label next to the arrow
  annotate(
    "text",
    x = x_12 - 10,
    y = y_12 - 12,
    label = TeX("$\\textbf{b}_{12}$"),
    color = "black",
    size = 6,
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = eval$x.1[40] + 10,
    y = eval$x.2[40] - 5,
    label = "Boundary",
    color = "black",
    size = 6,
    fontface = "bold"
  ) +
  
  labs(x = "Saber 11 Score", y = "Sisben Score") +
  coord_cartesian(xlim = c(-80, 100), ylim = c(-40, 60)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title   = element_text(size = 20, hjust = 0.5),
    text = element_text(family="Times New Roman", face="bold"),
    axis.text.x  = element_text(face = "bold", size = 12),
    axis.text.y  = element_text(face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "black", linetype = "solid")
  ) +
  xlab("Saber 11") +
  ylab("Sisben")

# --- arrow to b10: start from NW of the point, label at the tail ---
x0_1 <- x_1 - 8
y0_1 <- y_1 + 0
p1.1 <- p1.1 +
  geom_segment(
    aes(x = x0_1, y = y0_1, xend = x_1, yend = y_1),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black", size = 1
  ) +
  annotate(
    "text",
    x = x0_1 - 2, y = y0_1,               # near the arrow tail
    label = latex2exp::TeX("$\\textbf{b}_{1}$"),
    color = "black", size = 6, fontface = "bold",
    hjust = 1, vjust = 0
  )

# --- arrow to b30: start from NE of the point, label at the tail ---
x0_21 <- x_21 + 0
y0_21 <- y_21 - 8
p1.1 <- p1.1 +
  geom_segment(
    aes(x = x0_21, y = y0_21, xend = x_21, yend = y_21),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black", size = 1
  ) +
  annotate(
    "text",
    x = x0_21, y = y_21 - 15,               # near the arrow tail
    label = latex2exp::TeX("$\\textbf{b}_{21}$"),
    color = "black", size = 6, fontface = "bold",
    hjust = 0, vjust = 0
  )

print(p1.1)

ggsave("figures/fig2-a.png", p1.1, width = 6, height = 5)

################################ Scale Data ####################################

scale <- TRUE
if (scale){
  x.1.scale <- 1
  x.2.scale <- sd(data$x.1) / sd(data$x.2)
  data$x.1 <- data$x.1 * x.1.scale
  data$x.2 <- data$x.2 * x.2.scale
  eval$x.1 <- eval$x.1 * x.1.scale
  eval$x.2 <- eval$x.2 * x.2.scale
} else{
  x.1.scale <- 1
  x.2.scale <- 1
}

############################# SPP: Using Distance Method #######################

Y <- data$y
X <- cbind(data$x.1, data$x.2)
t <- data$d
b <- eval

D <- proxy::dist(X, eval, method = "euclidean")  # Use "euclidean" for Euclidean distances
t_expanded <- matrix(rep(2 * t - 1, times = ncol(D)), nrow = nrow(D), ncol = ncol(D))
D <- D * t_expanded

# smooth boundary
result.kinkoff <- rd2d.dist(Y, D, kink = "off", repp = 5000, vce = "hc1")
out.kinkoff <- cbind(result.kinkoff$results$h0, result.kinkoff$results$Est.p,
                     result.kinkoff$results$CI.lower, result.kinkoff$results$CI.upper,
                     result.kinkoff$results$CB.lower, result.kinkoff$results$CB.upper,
                     result.kinkoff$results$Est.q, result.kinkoff$results$z, result.kinkoff$results$`P>|z|`)

# unknown kink location
result.kinkon <- rd2d.dist(Y, D, kink = "on", repp = 5000, vce = "hc1", rbc = "off")
out.kinkon <- cbind(result.kinkon$results$h0, result.kinkon$results$Est.p,
                    result.kinkon$results$CI.lower, result.kinkon$results$CI.upper,
                    result.kinkon$results$CB.lower, result.kinkon$results$CB.upper,
                    result.kinkon$results$Est.q, result.kinkon$results$z, result.kinkon$results$`P>|z|`)

# adaptive
h.adaptive <- adaptive_bandwidth(result.kinkoff$results$h0,
                                 result.kinkon$results$h0,
                                 dist_to_kink)
h.adaptive <- cbind(h.adaptive, h.adaptive)
result.adaptive <- rd2d.dist(Y, D, kink = "off", repp = 5000, vce = "hc1", h = h.adaptive)
out.adaptive <- cbind(result.adaptive$results$h0, result.adaptive$results$Est.p,
                      result.adaptive$results$CI.lower, result.adaptive$results$CI.upper,
                      result.adaptive$results$CB.lower, result.adaptive$results$CB.upper,
                      result.adaptive$results$Est.q, result.adaptive$results$z, result.adaptive$results$`P>|z|`)

# rdrobust
bws <- matrix(0, nrow = neval, ncol = 2)
for (i in 1:neval) {
  out <- rdbwselect(Y, D[, i], vce = "hc1")
  bws[i, 1] <- out$bws[1]
  bws[i, 2] <- out$bws[2]
}
result.rdrobust <- rd2d.dist(Y, D, h = bws, kink = "off", repp = 5000, vce = "hc1")
out.rdrobust <- cbind(result.rdrobust$results$h0, result.rdrobust$results$Est.p,
                      result.rdrobust$results$CI.lower, result.rdrobust$results$CI.upper,
                      result.rdrobust$results$CB.lower, result.rdrobust$results$CB.upper,
                      result.rdrobust$results$Est.q, result.rdrobust$results$z, result.rdrobust$results$`P>|z|`)


################################################################################
### TABLE 6: BATEC Estimation and Inference (SPP Empirical Application)
################################################################################


dir.create("Results", showWarnings = FALSE)

library(xtable)

out_list <- list(
  smooth       = out.kinkoff,
  adaptive     = out.adaptive,
  unknown_kink = out.kinkon,
  rdrobust     = out.rdrobust
)

keep_idx <- 1:neval

fmt_ci  <- function(lo, hi) sprintf("$(%.3f,\\, %.3f)$", lo, hi)
fmt_tau <- function(j) sprintf("$\\mathbf{b}_{%d}$", j)

for (nm in names(out_list)) {
  
  out <- out_list[[nm]]
  
  tab <- data.frame(
    bb    = vapply(keep_idx, fmt_tau, character(1)),
    h     = out[keep_idx, 1],
    tau   = out[keep_idx, 2],
    pval  = out[keep_idx, 9],
    ci    = mapply(fmt_ci, out[keep_idx, 3], out[keep_idx, 4]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  xt <- xtable(
    tab,
    digits = c(0, 0, 3, 3, 3, 0),
    align  = c("l", "c", "r", "r", "r", "c")
  )
  
  lines <- c(
    "\\begin{tabular}{@{}crrrc@{}}",
    "\\toprule\\toprule",
    "\\multicolumn{1}{c}{$\\bb\\in\\B$} & \\multicolumn{1}{c}{$h$} & \\multicolumn{1}{c}{$\\tau(\\bb)$} & \\multicolumn{1}{c}{p-value} & \\multicolumn{1}{c}{95\\% RBC CI} \\\\",
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
  
  lines <- c(
    lines,
    body,
    "\\bottomrule\\bottomrule",
    "\\end{tabular}"
  )
  
  writeLines(lines, sprintf("tables/emp_app_%s.tex", nm))
}

################################################################################
### FIGURE 4: BATEC Estimation and Inference (SPP Empirical Application)
################################################################################

plot_rd2d_method <- function(out,
                             method_name,
                             show_idx = c(1, 5, 9, 11, 13, 17, 21),
                             kink_index = 11,
                             kink_x_text = 10,
                             kink_y = 0.05,
                             ylim = c(0.05, 0.55),
                             save_path = NULL) {
  
  indx <- 1:nrow(out)
  
  df <- data.frame(
    indx  = indx,
    y     = out[, 2],
    label = "BATEC"
  )
  
  df_ribbon <- data.frame(
    indx  = indx,
    ymin  = out[, 5],
    ymax  = out[, 6],
    label = "CB"
  )
  
  df_errorbar <- data.frame(
    indx  = indx,
    ymin  = out[, 3],
    ymax  = out[, 4],
    label = "CI"
  )
  
  temp_plot <- ggplot() + theme_bw()
  
  ## point estimates
  temp_plot <- temp_plot +
    geom_point(
      data = df,
      aes(x = indx, y = y,
          color = label, shape = label, fill = label, linetype = label)
    )
  
  ## confidence band
  temp_plot <- temp_plot +
    geom_ribbon(
      data = df_ribbon,
      aes(x = indx, ymin = ymin, ymax = ymax,
          color = label, fill = label, linetype = label),
      alpha = 0.1
    )
  
  ## confidence interval
  temp_plot <- temp_plot +
    geom_errorbar(
      data = df_errorbar,
      aes(x = indx, ymin = ymin, ymax = ymax,
          color = label, shape = label, fill = label, linetype = label)
    )
  
  temp_plot <- temp_plot +
    xlab("Cutoffs on the Boundary") +
    ylab("Treatment Effect") 
  
  legend_order <- c("BATEC", "CI", "CB")
  
  temp_plot <- temp_plot +
    scale_color_manual(
      values = c("BATEC" = "black", "CI" = "black", "CB" = "dodgerblue4"),
      name = NULL,
      breaks = legend_order
    ) +
    scale_shape_manual(
      values = c("BATEC" = 16, "CI" = 124, "CB" = 0),
      name = NULL,
      breaks = legend_order
    ) +
    scale_fill_manual(
      values = c("BATEC" = NA, "CI" = NA, "CB" = "dodgerblue4"),
      name = NULL,
      breaks = legend_order
    ) +
    scale_linetype_manual(
      values = c("BATEC" = 0, "CI" = 5, "CB" = 0),
      name = NULL,
      breaks = legend_order
    ) +
    guides(
      shape = "none",
      fill = "none",
      linetype = "none",
      color = guide_legend(
        order = 1,
        override.aes = list(
          shape    = c(16, NA, 22),
          linetype = c(0, 5, 0),
          fill     = c(NA, NA, "dodgerblue4"),
          alpha    = c(1, 1, 0.1),
          size     = c(2.5, 1, 5)
        )
      )
    )
  
  temp_plot <- temp_plot +
    geom_vline(
      xintercept = kink_index,
      color = "lightgrey",
      size = 1,
      linetype = "dotted"
    )
  
  temp_plot <- temp_plot +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      plot.title   = element_text(size = 20, hjust = 0.5),
      text         = element_text(family = "Times New Roman", face = "bold"),
      axis.text.x  = element_text(face = "bold", size = 15),
      axis.text.y  = element_text(face = "bold", size = 12),
      legend.position = c(0.8, 1),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  temp_plot <- temp_plot +
    scale_x_continuous(
      breaks = show_idx,
      labels = TeX(paste0("$\\textbf{b}_{", show_idx, "}$"))
    )
  
  temp_plot <- temp_plot +
    coord_cartesian(xlim = c(1, neval), ylim = ylim)
  
  print(temp_plot)
  
  if (!is.null(save_path)) {
    ggsave(save_path, temp_plot, width = 6, height = 5)
  }
  
  invisible(temp_plot)
}

################################## Produce Four Plots ##################################

show_idx <- c(1, 4, 7, 10, 12, 14, 17, 20)

p_smooth <- plot_rd2d_method(
  out = out.kinkoff,
  method_name = "Smooth Boundary",
  show_idx = show_idx,
  kink_index = 12,
  kink_x_text = 11,
  kink_y = -0.3,
  ylim = c(-0.3, 0.8),
  save_path = "figures/fig4-smooth.png"
)

p_adaptive <- plot_rd2d_method(
  out = out.adaptive,
  method_name = "Adaptive",
  show_idx = show_idx,
  kink_index = 12,
  kink_x_text = 11,
  kink_y = -0.3,
  ylim = c(-0.3, 0.8),
  save_path = "figures/fig4-adaptive.png"
)

p_kinkon <- plot_rd2d_method(
  out = out.kinkon,
  method_name = "Unknown Kink Location",
  show_idx = show_idx,
  kink_index = 12,
  kink_x_text = 11,
  kink_y = -0.3,
  ylim = c(-0.3, 0.8),
  save_path = "figures/fig4-unknown_kink.png"
)

p_rdrobust <- plot_rd2d_method(
  out = out.rdrobust,
  method_name = "Rdrobust",
  show_idx = show_idx,
  kink_index = 12,
  kink_x_text = 11,
  kink_y = -0.3,
  ylim = c(-0.3, 0.8),
  save_path = "figures/fig4-rdrobust.png"
)

################################################################################
### FIGURE 2: (b)-(d)
################################################################################

for (element in c(1:3)){
  idx <- c(1,12,21)[element]

  h <- result.kinkoff$results$h0[idx]
  hnumber <- idx
  
  ## 1) Run rdplot and extract bin means + smooth curves
  rdout <- rdplot(Y, D[,idx], nbins = c(500,500), p = 3)
  
  databins <- rdout$vars_bins %>%
    rename(x = rdplot_mean_x, y = rdplot_mean_y) %>%
    arrange(x)
  
  datapoly <- rdout$vars_poly %>%
    rename(x = rdplot_x, y = rdplot_y) %>%
    arrange(x)
  
  ## 2) Choose a symmetric bandwidth h (or keep your existing h)
  if (!exists("h")) {
    bw <- rdrobust(Y, D[,idx], p = 3)$bws  # c(h_left, h_right, ...)
    h  <- min(bw[1], bw[2], na.rm = TRUE) # symmetric choice
  }
  hnumber <- if (!exists("hnumber")) 1 else hnumber  # label like h_1 if you use it
  
  ## 3) Mark bins inside/outside the band
  databins <- databins %>%
    mutate(in_band = abs(x) <= h)
  
  ## 4) Colors/linetypes to match your scheme
  col_points_inside    <- "purple"
  col_points_outside   <- "grey"
  col_boundary_band    <- "purple"
  linety_boundary_band <- "solid"
  linewd_boundary_band <- 0.8
  
  ## 5) Build the plot (conditional means + smooth + band markers)
  plot <- ggplot() +
    # conditional means (colored by band)
    geom_point(data = databins,
               aes(x = x, y = y, color = in_band),
               alpha = 0.9, size = 1) +
    
    # regression curves (left/right, leave a tiny gap at 0 if you like)
    geom_path(data = subset(datapoly, x <= -0.01), aes(x = x, y = y),
              color = "black", linewidth = 1.1) +
    geom_path(data = subset(datapoly, x >=  0.01), aes(x = x, y = y),
              color = "black", linewidth = 1.1) +
    
    # band edges and cutoff
    geom_vline(xintercept =  h, color = col_boundary_band,
               linetype = linety_boundary_band, linewidth = linewd_boundary_band) +
    geom_vline(xintercept = -h, color = col_boundary_band,
               linetype = linety_boundary_band, linewidth = linewd_boundary_band) +
    geom_vline(xintercept =  0, color = "black", linetype = "dashed",
               linewidth = linewd_boundary_band) +
    
    # label the two halves of the band
    annotate("segment", x = -h, y = 0.90, xend =  0, yend = 0.90,
             arrow = arrow(length = unit(0.3, "cm"), ends = "both"),
             color = "black", linewidth = 1.2) +
    annotate("segment", x =  0, y = 0.90, xend =  h, yend = 0.90,
             arrow = arrow(length = unit(0.3, "cm"), ends = "both"),
             color = "black", linewidth = 1.2) +
    annotate("text", x = -h/2, y = 0.952, label = latex2exp::TeX(paste0("$\\hat{h}_{\\textbf{b}_{", hnumber, "}}$")),
             size = 6, fontface = "bold") +
    annotate("text", x =  h/2, y = 0.952, label = latex2exp::TeX(paste0("$\\hat{h}_{\\textbf{b}_{", hnumber, "}}$")),
             size = 6, fontface = "bold") +
    
    # legend mapping for inside/outside
    scale_color_manual(
      values = c(`TRUE` = col_points_inside, `FALSE` = col_points_outside),
      labels = c(`TRUE` = "Inside [-h, h]", `FALSE` = "Outside band"),
      name   = NULL
    ) +
    
    # axes, limits, theme
    coord_cartesian(xlim = c(-100, 100), ylim = c(0, 1)) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      plot.title   = element_text(size = 20, hjust = 0.5),
      text         = element_text(family = "Times New Roman", face = "bold"),
      axis.text.x  = element_text(face = "bold", size = 12),
      axis.text.y  = element_text(face = "bold", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    xlab(latex2exp::TeX(sprintf("Distance to Boundary Point $\\textbf{b}_{%d}$", idx))) +
    ylab("College Enrollment")
  
  print(plot)
  suffix <- c("b", "c", "d")
  
  ggsave(
    paste0("figures/fig2-", suffix[element], ".png"),
    plot,
    width = 6, height = 5
  )
}


