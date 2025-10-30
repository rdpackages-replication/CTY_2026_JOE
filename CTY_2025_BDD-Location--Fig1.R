
# rd2d: illustration file
# Authors: M. D. Cattaneo, R. Titiunik, R. R. Yu

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
  dist_to_kink[i] <- min(sqrt(eval[i,1]^2 + eval[i,2]^2), 75 - eval[i,1], 75 - eval[i,2])
}

####################### Scatter Plot and Boundary ############################## 

bound <- 21
point_1 <- eval[1,]
x_1 <- point_1$x.1
y_1 <- point_1$x.2

point_11 <- eval[11,]
x_11 <- point_11$x.1
y_11 <- point_11$x.2

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
  geom_segment(aes(x = 0, xend = 60, y = 0, yend = 0),
               linetype = "solid", color = "grey", alpha = 1, size = 3) +
  
  # Use scale_color_manual to define your own colors for d=0 and d=1
  scale_color_manual(
    name = NULL,  # Legend title (optional)
    values = c("0" = "#619CFF",  # Color for d=0
               "1" = "#F8766D"), # Color for d=1
    labels = c("0" = "Control", "1" = "Treatment")
  )

# Arrow pointing to point 21
p1.1 <- p1.1 +  geom_segment(
  aes(x = x_11 - 8, y = y_11 - 8, xend = x_11, yend = y_11),
  arrow = arrow(length = unit(0.2, "cm")),
  color = "black",
  size = 1
) +
  
  # Label next to the arrow
  annotate(
    "text",
    x = x_11 - 10,
    y = y_11 - 12,
    label = TeX("$\\textbf{b}_{11}$"),
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

ggsave("Results/fig-scatterX.png", p1.1, width = 6, height = 5)

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

# kink off
result.dist.kinkoff <- rd2d.dist(Y,D,kink = "off")
summary(result.dist.kinkoff,CBuniform = FALSE, subset = c(1,5,10,15,21,25,30,35,40))
tau.hat.rd2d.kinkoff <- result.dist.kinkoff$tau.hat

# kink on
result.dist.kinkon <- rd2d.dist(Y,D,kink = "on")
summary(result.dist.kinkon,CBuniform = FALSE, subset = c(1,5,10,15,21,25,30,35,40))
tau.hat.rd2d.kinkon <- result.dist.kinkon$tau.hat

################################## RD Plots ####################################

library(dplyr)
library(ggplot2)
library(latex2exp)
library(rdrobust)
library(grid)

idx <- 21 # 1, 11, 21
h <- result.dist.kinkoff$results$h0[idx]
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
ggsave(paste0("Results/fig-rdplot-b",idx,".png"), plot, width = 6, height = 5)

