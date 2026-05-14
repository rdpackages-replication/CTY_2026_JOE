################################################################################################
# Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods
# Empirical Application: SPP Data
# Figures
################################################################################################
rm(list=ls(all=TRUE))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active)) return(dirname(normalizePath(active)))
  }
  getwd()
}
setwd(get_script_dir())

suppressPackageStartupMessages({
  library(ggplot2)
  library(rdrobust)
  library(latex2exp)
  library(dplyr)
  library(grid)
})

output_dir <- "output"
figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
old_files <- list.files(figures_dir, pattern = "[.]png$", full.names = TRUE)
if (length(old_files) > 0) unlink(old_files)

get_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  as.integer(value)
}
emp_repp <- get_env_int("RD2D_EMP_REPP", 5000)

# Prefer ragg on Windows because the default PNG device is very slow for the
# large empirical scatter plots.
ggsave <- function(...) {
  args <- list(...)
  if (is.null(args$device) && requireNamespace("ragg", quietly = TRUE)) {
    args$device <- ragg::agg_png
  }
  do.call(ggplot2::ggsave, args)
}

################################## Load Data ###################################

# Load the SPP dataset using the JASA replication data schema.

load_spp_data <- function(path = "spp.csv") {
  raw <- read.csv(path)
  expected <- c(
    "running_saber11",
    "running_sisben",
    "eligible_spp",
    "beneficiary_spp",
    "spadies_any",
    "icfes_educm1"
  )
  missing_cols <- setdiff(expected, names(raw))
  if (length(missing_cols) > 0) {
    stop(
      "spp.csv is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  dat <- raw[, c("running_saber11", "running_sisben", "spadies_any", "eligible_spp")]
  names(dat) <- c("x.1", "x.2", "y", "d")
  dat <- dat[complete.cases(dat), , drop = FALSE]

  expected_assignment <- as.integer(dat$x.1 >= 0 & dat$x.2 >= 0)
  if (!all(dat$d == expected_assignment)) {
    stop("eligible_spp does not match the quadrant assignment rule.", call. = FALSE)
  }

  dat$y <- as.numeric(dat$y)
  dat$d <- as.numeric(dat$d)
  dat
}

data <- load_spp_data("spp.csv")

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

# subsetting so the kink (0, 0) is the 11th retained boundary point
eval <- eval[c(11:31),]
neval <- 21


dir.create("figures", showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(extrafont))

# Load fonts into R session
invisible(utils::capture.output(suppressMessages(loadfonts(device = "pdf"))))

# # Define the piecewise function
# f <- function(x) {
#   ifelse(x <= 3/4,
#          2/pi * x,
#          (x + 3/4)/(pi - acos(3/(4 * x)))
#   )
# }

bound_vec <- c(0:6)

for (bound in bound_vec){

# Define the piecewise function
f <- function(x) {
  out <- numeric(length(x))
  left <- x <= 1/2
  out[left] <- 2/pi * x[left]
  out[!left] <- (x[!left] + 1/2)/(pi - acos(1/(2 * x[!left])))
  out
}

point_x <- c(0.2,0.4,0.5,0.6,0.8,1)
labels <- c(TeX("$r_1$"),TeX("$r_2$"),TeX("$r_3$"),TeX("$r_4$"),TeX("$r_5$"),TeX("$r_6$"))
labels_y <- c(TeX("$\\theta_{1,\\textbf{b}}(r_1)$"),TeX("$\\theta_{1,\\textbf{b}}(r_2)$"),
              TeX("$\\theta_{1,\\textbf{b}}(r_3)$"),TeX("$\\theta_{1,\\textbf{b}}(r_4)$"),
              TeX("$\\theta_{1,\\textbf{b}}(r_5)$"),TeX("$\\theta_{1,\\textbf{b}}(r_6)$"))

if (bound > 0){
  point_x <- c(0,point_x[c(1:bound)])
  labels <- c(TeX("$0$"), labels[c(1:bound)])
  labels_y <- c(TeX("$\\theta_{1,\\textbf{b}}(0)$"), labels_y[c(1:bound)])
} else {
  point_x <- c(0)
  labels <- c(TeX("$0$"))
  labels_y <- c(TeX("$\\theta_{1,\\textbf{b}}(0)$"))
}

point_y <- f(point_x)
point_data <- data.frame(x = point_x, y = point_y)



# Create a sequence of x values from 0 to 1
x_values <- seq(0, tail(point_x, n = 1), length.out = 1000)

# Calculate y values using the piecewise function
y_values <- sapply(x_values, f)

# Create a data frame for plotting
plot_data <- data.frame(x = x_values, y = y_values)


# Plot the function using ggplot2
plot <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_line(color = "dimgrey", linewidth = 1) +
  # labs(
  #   x = TeX("$$"),
  #   y = TeX("$$")
  # ) +
  # labs(
  #   x = NULL,
  #   y = NULL
  # ) +
  geom_point(data = point_data,
             aes(x = x, y = y),
             color = "blue",
             size = 3) +
  geom_segment(data = point_data, aes(x = x, y = 0, xend = x, yend = y),
               linetype = "dashed", color = "lightgrey") +
  geom_segment(data = point_data, aes(x = 0, y = y, xend = x, yend = y),
               linetype = "dashed", color = "lightgrey") +
  # Annotate axis labels at the ends
  # annotate("text", x = 1 , y = 0, label = TeX("$r$"), hjust = 0, size = 5) +
  # annotate("text", x = 0, y = f(1), label = TeX("$\\theta(r)$"), vjust = 0, size = 5) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(family = "Times", hjust = 0.5, size = 20, ,face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(family = "Times", size = 15,face = "bold"),
    axis.text.y = element_text(family = "Times", size = 12,face = "bold"),
    plot.margin = margin(20, 20, 20, 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5,
                               arrow = grid::arrow(type = "closed", length = unit(0.15, "inches"))),
    axis.line.y = element_line(color = "black", linewidth = 0.5,
                               arrow = grid::arrow(type = "closed", length = unit(0.15, "inches")))
  ) +
  coord_cartesian(
    xlim = c(0, 1),  # Adjust X-axis limits here
    ylim = c(0, 0.75)   # Adjust Y-axis limits here
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = point_x,
    labels = labels # Using latex2exp for the LaTeX-style label
  ) +
  scale_y_continuous(limits = c(0, 0.75),
                     breaks = point_y,
                     labels = labels_y)

if (interactive()) print(plot)

# Save the plot as a PDF
# ggsave(sprintf("figures/conditional_mean_in_r_p%d.png", bound), plot, width = 6, height = 5)
if (bound == 6){
  ggsave("figures/fig1-b.png", width = 6, height = 5)
}
}


################################################################################
# Cartesian View
################################################################################

library(ggplot2)

bound_vec <- c(0:6)

for (bound in bound_vec){

# Circle parameters
cx <- 0.5      # Center x
cy <- 0.0      # Center y

r_vec <- c(0.2, 0.4, 0.5, 0.6, 0.8, 1)
labels <- c(TeX("$r_1$"),TeX("$r_2$"),TeX("$r_3$"),TeX("$r_4$"),TeX("$r_5$"),TeX("$r_6$"))
theta_vec <- c(7 * pi/8, 5 * pi/8, 4 * pi/8, 3 * pi/8, 2 * pi/8, 1 * pi/8)
ex_vec <- cx + r_vec * cos(theta_vec)
ey_vec <- cy + r_vec * sin(theta_vec)
theta_end_vec <- c(pi, pi, pi, pi - acos(0.5/0.6), pi - acos(0.5/0.8), pi - acos(0.5/1))
x_adjust_vec <- c(0.04,0.04,0.05,0.06,0.1,0.06)
y_adjust_vec <- c(0.03,0.03,0.03,0.03,0.03,-0.03)

plot <- ggplot() +
  # 1) Draw x- and y-axis lines (through the origin)
  # geom_segment(aes(x = 0, xend = 0, y = -0.018, yend = 1), linetype = "solid", color = "blue", alpha = 0.5, linewidth = 3) +  # Vertical red line
  # geom_segment(aes(x = -0.018, xend = 1.55, y = 0, yend = 0), linetype = "solid", color = "blue", alpha = 0.5, linewidth = 3) +  # Horizontal red line
  #
  geom_segment(aes(x = 0, xend = 0, y = -0.018, yend = 1), linetype = "solid", color = "grey", alpha = 0.5, linewidth = 3) +  # Vertical red line
  geom_segment(aes(x = -0.018, xend = 1.55, y = 0, yend = 0), linetype = "solid", color = "grey", alpha = 0.5, linewidth = 3) +  # Horizontal red line

  geom_point(aes(x = 0.5, y = 0), color = "black", size = 3) +

  annotate("text",
           x = 0.5,   # Shift slightly right
           y = 0,    # Midpoint in y-direction
           label = TeX("$\\textbf{b}$"),
           color = "black",
           size = 5,
           vjust = 1.7)

if (bound > 0){
  for (i in 1:bound){
    theta <- seq(0, theta_end_vec[i], length.out = 200)
    df_circle <- data.frame(
      x = cx + r_vec[i] * cos(theta),
      y = cy + r_vec[i] * sin(theta))
    # plot <- plot + geom_path(data = df_circle, aes(x = x, y = y),
    #             color = "grey", size = 1)
    plot <- plot + geom_path(data = df_circle, aes(x = x, y = y),
                             color = "blue", linewidth = 1)

    plot <- plot + geom_segment(aes(x = 0.5, y = 0),
                  xend = ex_vec[i], yend = ey_vec[i],
                   arrow = arrow(length = unit(0.2, "cm")),
                   color = "black", linewidth = 1)

    plot <- plot + annotate("text",
                 x = 2/6 * cx + 4/6 * ex_vec[i] + x_adjust_vec[i],   # Shift slightly right
                 y = 2/6 * cy + 4/6 * ey_vec[i] + y_adjust_vec[i],    # Midpoint in y-direction
                 label = labels[i],
                 color = "black",
                 size = 5)
  }
}

plot <- plot +
  coord_fixed(xlim = c(-0.2, 1.7), ylim = c(-0.2, 1.2)) +
  # scale_x_continuous(
  #   limits = c(-0.2, 1.7),
  #   breaks = c(0.75),
  #   labels = c(TeX("$X_1$")) # Using latex2exp for the LaTeX-style label
  # ) +
  # scale_y_continuous(limits = c(-0.2, 1.2),
  #                    breaks = c(0.5),
  #                    labels = c(TeX("$X_2$")))


  # Optional axis labels
  # labs(x = TeX("$X_1$"), y = TeX("$X_2$")) +
  labs(x = TeX("$X_{1i}$"), y = TeX("$X_{2i}$")) +

  # Some minimal styling
  theme_minimal(base_size = 14) +
  # Remove extra grid lines if desired
  theme(
    plot.margin = margin(3, 3, 3, 3),
    panel.background = element_rect(fill = "white", color = NA),  # Set white background
    plot.background = element_rect(fill = "white", color = NA),   # Set white background
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(family = "Times", size = 15,face = "bold"),
    # axis.title = element_blank(),
    axis.ticks = element_blank(),  # Remove tick marks
    axis.text.x = element_blank(),  # Remove x-axis tick labels
    axis.text.y = element_blank(),   # Remove y-axis tick labels
    # axis.text.x = element_text(family = "Times", size = 15,face = "bold"),
    # axis.text.y = element_text(family = "Times", size = 15,face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5),  # Title size and centering
    axis.line.x = element_line(color = "black", linewidth = 0.5,
                               arrow = grid::arrow(type = "closed", length = unit(0.15, "inches"))),
    axis.line.y = element_line(color = "black", linewidth = 0.5,
                               arrow = grid::arrow(type = "closed", length = unit(0.15, "inches")))
  )

if (interactive()) print(plot)
# ggsave(sprintf("figures/conditional_mean_cartesian_p%d.png", bound), plot, width = 6, height = 5)
if (bound == 6){
  ggsave("figures/fig1-a.png", width = 6, height = 5)
}
}


################################################################################
### FIGURE 1: Scatter plot score and RDPLOTs
################################################################################

## Scatter Plot and Boundary

scatter_label_idx <- c(1, 4, 7, 14, 17, 21)
eval_labeled <- eval %>%
  mutate(
    index = row_number(),
    label = sprintf("bold(b)[%d]", index),
    label_x = if_else(x.1 == 0, -3.2, x.1),
    label_y = if_else(x.1 == 0, x.2, -2.4),
    label_hjust = if_else(x.1 == 0, 1, 0.5),
    label_vjust = if_else(x.1 == 0, 0.5, 1)
  ) %>%
  filter(index %in% scatter_label_idx)

kink_point <- eval[11,]
boundary_segments <- data.frame(
  x = c(0, 0),
  xend = c(0, 80),
  y = c(0, 0),
  yend = c(55, 0),
  group = "Boundary",
  stringsAsFactors = FALSE
)

p1.1 <- ggplot() +
  geom_point(
    data = data  %>% sample_frac(0.3),
    aes(x = x.1, y = x.2, color = factor(d), shape = factor(d)),
    alpha = 0.18,
    size = 0.35
  ) +
  geom_segment(
    data = boundary_segments,
    aes(x = x, xend = xend, y = y, yend = yend, color = group),
    linewidth = 1.0,
    linetype = "solid",
    lineend = "butt"
  ) +
  geom_point(
    data = eval,
    aes(x = x.1, y = x.2),
    size = 0.65
  ) +
  geom_text(
    data = eval_labeled,
    aes(
      x = label_x,
      y = label_y,
      label = label,
      hjust = label_hjust,
      vjust = label_vjust
    ),
    color = "black",
    size = 2.7,
    parse = TRUE
  ) +
  geom_segment(
    aes(
      x = kink_point$x.1 - 5,
      y = kink_point$x.2 - 5,
      xend = kink_point$x.1 - 1.4,
      yend = kink_point$x.2 - 1.4
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = kink_point$x.1 - 6.5,
    y = kink_point$x.2 - 7,
    label = TeX("$\\textbf{b}_{11}$"),
    color = "black",
    size = 3.6,
    fontface = "bold",
    hjust = 1
  ) +
  scale_color_manual(
    values = c("0" = "indianred2", "1" = "dodgerblue4", "Boundary" = "grey45"),
    breaks = c("0", "1", "Boundary"),
    name = NULL,
    labels = c("Control", "Treatment", "Boundary")
  ) +
  scale_shape_manual(
    values = c("0" = 15, "1" = 16),
    name = NULL,
    labels = c("Control", "Treatment")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        alpha = c(1, 1, 1),
        size = c(2, 2, 0),
        shape = c(15, 16, NA),
        linetype = c("blank", "blank", "solid"),
        linewidth = c(0, 0, 1)
      )
    ),
    shape = "none"
  ) +
  labs(x = "Saber 11 Score", y = "Sisben Score") +
  coord_cartesian(xlim = c(-80, 100), ylim = c(-40, 60)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title   = element_text(size = 20, hjust = 0.5),
    text = element_text(family = "serif", face = "bold"),
    axis.text.x  = element_text(face = "bold", size = 12),
    axis.text.y  = element_text(face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.23, 0.88),
    legend.justification = c(0.5, 1),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  xlab("Saber 11") +
  ylab("Sisben")

if (interactive()) print(p1.1)

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
assignment <- data$d
b <- eval

distance <- proxy::dist(X, eval, method = "euclidean")  # Use "euclidean" for Euclidean distances
assignment_expanded <- matrix(rep(2 * assignment - 1, times = ncol(distance)), nrow = nrow(distance), ncol = ncol(distance))
distance <- distance * assignment_expanded


################################## Read Output #################################

read_empapp_output <- function(id) {
  path <- file.path(output_dir, paste0("empapp_", id, ".csv"))
  if (!file.exists(path)) {
    stop(sprintf("Missing output file: %s. Run CTY_2026_JOE--empapp.R first.", path), call. = FALSE)
  }
  tab <- read.csv(path, check.names = FALSE)
  required <- c("h0", "h1", "estimate.p", "ci.lower", "ci.upper", "cb.lower", "cb.upper", "estimate.q", "t.value", "p.value")
  missing <- setdiff(required, names(tab))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing columns: %s", path, paste(missing, collapse = ", ")), call. = FALSE)
  }
  if ("row" %in% names(tab)) {
    tab <- tab[!(as.character(tab$row) %in% c("WBATE", "LBATE")), , drop = FALSE]
  }
  plot_cols <- c("h0", "estimate.p", "ci.lower", "ci.upper", "cb.lower", "cb.upper", "estimate.q", "t.value", "p.value")
  as.matrix(tab[, plot_cols, drop = FALSE])
}
out.kinkoff <- read_empapp_output("smooth_itt")
out.adaptive <- read_empapp_output("adaptive_itt")
out.kinkon <- read_empapp_output("unknown_kink_itt")
out.rdrobust <- read_empapp_output("rdrobust_itt")
result.kinkoff <- list(results = data.frame(h0 = out.kinkoff[, 1]))
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

  legend_order <- c("BATEC", "CI", "CB")

  temp_plot <- ggplot(df, aes(x = indx)) +
    geom_vline(
      xintercept = kink_index,
      color = "grey90",
      linewidth = 0.45
    ) +
    geom_ribbon(
      data = df_ribbon,
      aes(x = indx, ymin = ymin, ymax = ymax, fill = label, color = label),
      alpha = 0.14,
      linewidth = 0
    ) +
    geom_errorbar(
      data = df_errorbar,
      aes(x = indx, ymin = ymin, ymax = ymax, color = label),
      width = 0.12,
      linewidth = 0.35
    ) +
    geom_point(aes(y = y, color = label), size = 1.2) +
    xlab("Cutoffs on the Boundary") +
    ylab("Treatment Effect") +
    scale_color_manual(
      values = c("BATEC" = "black", "CI" = "black", "CB" = "dodgerblue4"),
      name = NULL,
      breaks = legend_order
    ) +
    scale_fill_manual(
      values = c("CB" = "dodgerblue4"),
      name = NULL
    ) +
    guides(
      color = guide_legend(
        order = 1,
        override.aes = list(
          shape = c(16, NA, 22),
          linetype = c("blank", "solid", "blank"),
          fill = c(NA, NA, "dodgerblue4"),
          alpha = c(1, 1, 0.14),
          linewidth = c(0, 0.35, 0)
        )
      ),
      fill = "none"
    )

  temp_plot <- temp_plot +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      plot.title   = element_text(size = 20, hjust = 0.5),
      text         = element_text(family = "serif", face = "bold"),
      axis.text.x  = element_text(face = "bold", size = 15),
      axis.text.y  = element_text(face = "bold", size = 12),
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
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

  if (interactive()) print(temp_plot)

  if (!is.null(save_path)) {
    ggsave(save_path, temp_plot, width = 6, height = 5)
  }

  invisible(temp_plot)
}

################################## Produce Four Plots ##################################

show_idx <- c(1, 4, 7, 11, 14, 17, 21)

p_smooth <- plot_rd2d_method(
  out = out.kinkoff,
  method_name = "Smooth Boundary",
  show_idx = show_idx,
  kink_index = 11,
  kink_x_text = 11,
  kink_y = -0.3,
  ylim = c(-0.3, 0.8),
  save_path = "figures/fig4-smooth.png"
)

p_adaptive <- plot_rd2d_method(
  out = out.adaptive,
  method_name = "Adaptive",
  show_idx = show_idx,
  kink_index = 11,
  kink_x_text = 11,
  kink_y = -0.3,
  ylim = c(-0.3, 0.8),
  save_path = "figures/fig4-adaptive.png"
)

p_kinkon <- plot_rd2d_method(
  out = out.kinkon,
  method_name = "Unknown Kink Location",
  show_idx = show_idx,
  kink_index = 11,
  kink_x_text = 11,
  kink_y = -0.3,
  ylim = c(-0.3, 0.8),
  save_path = "figures/fig4-unknown_kink.png"
)

p_rdrobust <- plot_rd2d_method(
  out = out.rdrobust,
  method_name = "Rdrobust",
  show_idx = show_idx,
  kink_index = 11,
  kink_x_text = 11,
  kink_y = -0.3,
  ylim = c(-0.3, 0.8),
  save_path = "figures/fig4-rdrobust.png"
)

################################################################################
### FIGURE 2: (b)-(d)
################################################################################

for (element in c(1:3)){
  idx <- c(1,11,21)[element]

  h <- result.kinkoff$results$h0[idx]
  hnumber <- idx

  ## 1) Run rdplot and extract bin means + smooth curves
  rdout <- rdplot(Y, distance[,idx], nbins = c(500,500), p = 3, hide = TRUE)

  databins <- rdout$vars_bins %>%
    rename(x = rdplot_mean_x, y = rdplot_mean_y) %>%
    arrange(x)

  datapoly <- rdout$vars_poly %>%
    rename(x = rdplot_x, y = rdplot_y) %>%
    arrange(x)

  ## 2) Choose a symmetric bandwidth h (or keep your existing h)
  if (!exists("h")) {
    bw <- rdrobust(Y, distance[,idx], p = 3)$bws  # c(h_left, h_right, ...)
    h  <- min(bw[1], bw[2], na.rm = TRUE) # symmetric choice
  }
  hnumber <- if (!exists("hnumber")) 1 else hnumber  # label like h_1 if you use it

  ## 3) Mark bins inside/outside the band
  databins <- databins %>%
    mutate(
      in_band = abs(x) <= h,
      band_side = dplyr::case_when(
        !in_band ~ "Outside band",
        x < 0 ~ "Control",
        TRUE ~ "Treatment"
      )
    )

  ## 4) Colors/linetypes to match your scheme
  col_control          <- "indianred2"
  col_treatment        <- "dodgerblue4"
  col_points_outside   <- "grey"
  linety_boundary_band <- "solid"
  linewd_boundary_band <- 0.8

  ## 5) Build the plot (conditional means + smooth + band markers)
  plot <- ggplot() +
    # conditional means (active bins colored by side)
    geom_point(data = databins,
               aes(x = x, y = y, color = band_side),
               alpha = 0.9, size = 1) +

    # regression curves (left/right, leave a tiny gap at 0 if you like)
    geom_path(data = subset(datapoly, x <= -0.01), aes(x = x, y = y),
              color = "black", linewidth = 1.1) +
    geom_path(data = subset(datapoly, x >=  0.01), aes(x = x, y = y),
              color = "black", linewidth = 1.1) +

    # band edges and cutoff
    geom_vline(xintercept =  h, color = col_treatment,
               linetype = linety_boundary_band, linewidth = linewd_boundary_band) +
    geom_vline(xintercept = -h, color = col_control,
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
      values = c(
        "Control" = col_control,
        "Treatment" = col_treatment,
        "Outside band" = col_points_outside
      ),
      name   = NULL
    ) +

    # axes, limits, theme
    coord_cartesian(xlim = c(-100, 100), ylim = c(0, 1)) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      plot.title   = element_text(size = 20, hjust = 0.5),
      text         = element_text(family = "serif", face = "bold"),
      axis.text.x  = element_text(face = "bold", size = 12),
      axis.text.y  = element_text(face = "bold", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    xlab(latex2exp::TeX(sprintf("Distance to Boundary Point $\\textbf{b}_{%d}$", idx))) +
    ylab("College Enrollment")

  if (interactive()) print(plot)
  suffix <- c("b", "c", "d")

  ggsave(
    paste0("figures/fig2-", suffix[element], ".png"),
    plot,
    width = 6, height = 5
  )
}

cat(sprintf("Wrote empirical figure files to %s.\n", figures_dir))
