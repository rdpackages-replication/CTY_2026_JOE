## Replication Script: Simulation Figures
## Cattaneo, Titiunik and Yu (2026): Boundary Distance RD

options(scipen = 999, digits = 4)

output_dir <- Sys.getenv("RD2D_OUTPUT_DIR", unset = "output")
figures_dir <- Sys.getenv("RD2D_FIGURES_DIR", unset = "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

methods <- c("kinkoff", "adaptive", "kinkon", "rdrobustadj")
method_labels <- c(
  kinkoff = "Smooth",
  adaptive = "Adaptive",
  kinkon = "Unknown Kink",
  rdrobustadj = "Rdrobust"
)
method_cols <- c(
  kinkoff = "blueviolet",
  adaptive = "grey45",
  kinkon = "blue",
  rdrobustadj = "darkgreen"
)
method_pch <- c(kinkoff = 16, adaptive = 15, kinkon = 4, rdrobustadj = 17)
method_lty <- c(kinkoff = 1, adaptive = 5, kinkon = 2, rdrobustadj = 6)

figure_specs <- data.frame(
  dgp = 1:4,
  tag = c(
    "linear_homoskedastic",
    "quadratic_homoskedastic",
    "linear_heteroskedastic",
    "quadratic_heteroskedastic"
  ),
  row = c("linear", "quadratic", "linear", "quadratic"),
  col = c("homoskedastic", "homoskedastic", "heteroskedastic", "heteroskedastic"),
  stringsAsFactors = FALSE
)

## Input helpers
clean_figures_dir <- function(path) {
  files <- list.files(path, pattern = "^(simuls_)?fig3_.*[.]png$", full.names = TRUE)
  if (length(files)) invisible(file.remove(files))
}

read_output_csv <- function(file) {
  path <- file.path(output_dir, file)
  if (!file.exists(path)) stop(sprintf("Missing output file: %s. Run CTY_2026_JOE--simuls.R first.", path), call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

metadata_value <- function(metadata, name) {
  if (!(name %in% names(metadata))) stop(sprintf("Simulation metadata is missing column '%s'.", name), call. = FALSE)
  metadata[[name]][1]
}

read_raw_file <- function(file) {
  if (!file.exists(file)) stop(sprintf("Missing simulation file: %s", file), call. = FALSE)
  read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
}

read_simulation_rows <- function(m, dgp, method, tag) {
  files <- file.path(
    output_dir,
    sprintf("simuls_raw_rep%04d_dgp%02d_sharp_%s_%s_bwmain.csv", seq_len(m), dgp, method, tag)
  )
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) stop(sprintf("Missing simulation file: %s", missing_files[1]), call. = FALSE)
  do.call(rbind, lapply(files, read_raw_file))
}

open_png <- function(file, width, height, units = "in", res = 300) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = file, width = width, height = height, units = units, res = res)
  } else {
    grDevices::png(file, width = width, height = height, units = units, res = res)
  }
}

## Population curves from calibrated ITT coefficients
design_matrix <- function(x, s) {
  if (s == 1L) {
    return(cbind(
      "(Intercept)" = 1,
      "x.1" = x$x.1,
      "x.2" = x$x.2
    ))
  }
  cbind(
    "(Intercept)" = 1,
    "x.1" = x$x.1,
    "x.2" = x$x.2,
    "I(x.1^2)" = x$x.1^2,
    "I(x.1 * x.2)" = x$x.1 * x$x.2,
    "I(x.2^2)" = x$x.2^2
  )
}

coefficient_vector <- function(parameters, tag, side, s) {
  Xnames <- colnames(design_matrix(data.frame(x.1 = 0, x.2 = 0), s))
  out <- setNames(numeric(length(Xnames)), Xnames)
  rows <- parameters[
    parameters$design == "sharp" &
      parameters$tag == tag &
      parameters$component == "beta_y" &
      parameters$side == side,
  ]
  for (term in Xnames) {
    val <- rows$value[rows$term == term]
    if (length(val) == 1L) out[term] <- val
  }
  out
}

population_curve <- function(parameters, boundary, tag, s) {
  X <- design_matrix(boundary, s)
  beta0 <- coefficient_vector(parameters, tag, 0, s)
  beta1 <- coefficient_vector(parameters, tag, 1, s)
  as.numeric(X %*% (beta1[colnames(X)] - beta0[colnames(X)]))
}

## Monte Carlo summaries
summarize_estimate <- function(rows) {
  rows <- rows[rows$estimand == "sharp" & !(rows$row %in% c("WBATE", "LBATE")), ]
  rows$index <- as.integer(rows$row)
  split_rows <- split(rows, rows$index)
  out <- do.call(rbind, lapply(split_rows, function(z) {
    data.frame(
      index = z$index[1],
      mean = mean(z$estimate.p, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$index), ]
}

plot_figure3_panel <- function(truth, plot_data, file, ylim, yticks) {
  open_png(file, width = 6, height = 5, units = "in", res = 300)
  on.exit(grDevices::dev.off(), add = TRUE)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mar = c(4.5, 4.6, 0.8, 0.6), family = "serif", bty = "n")

  x <- seq_along(truth)
  show_idx <- c(1, 4, 7, 11, 14, 17, 21)
  x_labels <- as.expression(lapply(show_idx, function(i) bquote(bold(b)[.(i)])))

  graphics::plot(
    x, truth,
    type = "l", lwd = 2.4, col = "black", lty = 1,
    ylim = ylim,
    axes = FALSE,
    xlab = "Cutoffs on the Boundary",
    ylab = "Treatment Effect",
    main = "",
    frame.plot = FALSE
  )
  graphics::axis(1, at = show_idx, labels = x_labels, tick = FALSE, line = 0)
  graphics::axis(2, at = yticks, tick = FALSE, line = 0, las = 2)
  graphics::abline(v = 11, col = "lightgrey", lwd = 1.8, lty = "solid")

  for (method in methods) {
    z <- plot_data[[method]]
    graphics::lines(z$index, z$mean, col = method_cols[[method]], lwd = 1.25, lty = method_lty[[method]])
    graphics::points(z$index, z$mean, col = method_cols[[method]], pch = method_pch[[method]], cex = 0.75)
  }

  graphics::legend(
    "topright",
    legend = c(unname(method_labels[methods]), "Population"),
    col = c(unname(method_cols[methods]), "black"),
    pch = c(unname(method_pch[methods]), NA),
    lty = c(unname(method_lty[methods]), 1),
    lwd = c(rep(1.25, length(methods)), 2.4),
    bty = "n",
    cex = 0.82
  )
}

check_generated_figures <- function(expected) {
  missing <- expected[!file.exists(expected)]
  if (length(missing)) stop(sprintf("Missing expected figure: %s", missing[1]), call. = FALSE)
  empty <- expected[file.info(expected)$size <= 0]
  if (length(empty)) stop(sprintf("Empty generated figure: %s", empty[1]), call. = FALSE)
  invisible(expected)
}

## Run
clean_figures_dir(figures_dir)
config <- read_output_csv("simuls_config.csv")
m <- as.integer(metadata_value(config, "m"))
boundary <- read_output_csv("simuls_eval_grid.csv")
parameters <- read_output_csv("simuls_dgp_parameters.csv")

expected <- character(0)
for (i in seq_len(nrow(figure_specs))) {
  spec <- figure_specs[i, ]
  s <- if (identical(spec$row, "linear")) 1L else 2L
  truth <- population_curve(parameters, boundary, spec$tag, s)
  method_rows <- lapply(methods, function(method) {
    summarize_estimate(read_simulation_rows(m, spec$dgp, method, spec$tag))
  })
  names(method_rows) <- methods

  file <- file.path(figures_dir, sprintf("fig3_%s.png", spec$tag))
  ylim <- if (identical(spec$row, "linear")) c(0.34, 0.39) else c(0.30, 0.39)
  yticks <- if (identical(spec$row, "linear")) c(0.34, 0.365, 0.39) else c(0.30, 0.345, 0.39)
  plot_figure3_panel(truth, method_rows, file, ylim, yticks)
  expected <- c(expected, file)
  cat(sprintf("Wrote %s.\n", basename(file)))
}

check_generated_figures(expected)
cat("Simulation Figure 3 files complete.\n")
