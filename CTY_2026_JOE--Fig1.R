################################################################################################
# Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods
# Figure 1
# Authors: M. D. Cattaneo, R. Titiunik, R. R. Yu
################################################################################################

library(ggplot2)
library(latex2exp)
library(extrafont)
library(grid) # For custom drawing

# Load fonts into R session
loadfonts(device = "pdf")  # Use "win" for Windows, "pdf" for PDF devices

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
  ifelse(x <= 1/2, 
         2/pi * x, 
         (x + 1/2)/(pi - acos(1/(2 * x)))
  )
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
data <- data.frame(x = x_values, y = y_values)


# Plot the function using ggplot2
plot <- ggplot(data, aes(x = x, y = y)) +
  geom_line(color = "dimgrey", size = 1) +
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
    axis.line.x = element_line(color = "black", size = 0.5, 
                               arrow = grid::arrow(type = "closed", length = unit(0.15, "inches"))),
    axis.line.y = element_line(color = "black", size = 0.5, 
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

print(plot)

# Save the plot as a PDF
# ggsave(sprintf("Results/conditional_mean_in_r_p%d.png",bound), plot, width = 6, height = 5)
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
                             color = "blue", size = 1) 
    print(ex_vec[i])
    plot <- plot + geom_segment(aes(x = 0.5, y = 0),
                  xend = ex_vec[i], yend = ey_vec[i],              
                   arrow = arrow(length = unit(0.2, "cm")),
                   color = "black", size = 1) 
        
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
    plot.title = element_text(size = , hjust = 0.5),  # Title size and centering
    axis.line.x = element_line(color = "black", size = 0.5, 
                               arrow = grid::arrow(type = "closed", length = unit(0.15, "inches"))),
    axis.line.y = element_line(color = "black", size = 0.5, 
                               arrow = grid::arrow(type = "closed", length = unit(0.15, "inches")))
  ) 

print(plot)
# ggsave(sprintf("Results/conditional_mean_cartesian_p%d.png",bound), plot, width = 6, height = 5)
if (bound == 6){
  ggsave("figures/fig1-a.png", width = 6, height = 5)
}
}

