## Replication Script: Simulation Tables
## Cattaneo, Titiunik and Yu (2026): Boundary Distance RD

options(scipen = 999, digits = 4)

output_dir <- Sys.getenv("RD2D_OUTPUT_DIR", unset = "output")
tables_dir <- Sys.getenv("RD2D_TABLES_DIR", unset = "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

estimand_specs <- data.frame(
  design = c("fuzzy", "fuzzy", "fuzzy", "sharp", "sharp"),
  output = c("main", "itt", "fs", "main", "main.0"),
  estimand = c("fuzzy", "itt", "fs", "sharp", "sharp0"),
  stringsAsFactors = FALSE
)
estimands <- estimand_specs$estimand
estimand_labels <- c(
  fuzzy = "Fuzzy",
  itt = "Intention-to-treat",
  fs = "First stage",
  sharp = "Sharp",
  sharp0 = "Sharp control side"
)

methods <- c("kinkoff", "adaptive", "kinkon", "rdrobustadj")
method_labels <- c(
  kinkoff = "Smooth",
  adaptive = "Adaptive",
  kinkon = "Unknown kink",
  rdrobustadj = "RD robust adjusted"
)

bandwidth_specs <- data.frame(
  bwparam = c("main", "itt"),
  suffix = c("bwmain", "bwitt"),
  stringsAsFactors = FALSE
)

## Input helpers
clean_tables_dir <- function(path) {
  files <- list.files(path, pattern = "^simuls_.*[.]tex$", full.names = TRUE)
  if (length(files)) invisible(file.remove(files))
}

read_output_csv <- function(file) {
  path <- file.path(output_dir, file)
  if (!file.exists(path)) stop(sprintf("Missing output file: %s. Run CTY_2026_JOE--simuls.R first.", path), call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

metadata_value <- function(metadata, name, default = NULL) {
  if (!(name %in% names(metadata))) {
    if (!is.null(default)) return(default)
    stop(sprintf("Simulation metadata is missing column '%s'.", name), call. = FALSE)
  }
  metadata[[name]][1]
}

read_raw_file <- function(file) {
  if (!file.exists(file)) stop(sprintf("Missing simulation file: %s", file), call. = FALSE)
  dat <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  expected <- c(
    "replication", "dgp", "tag", "method", "bwparam", "suffix", "design", "estimand", "row",
    "h0", "h1", "h0.rbc", "h1.rbc", "estimate.p", "std.err.p",
    "estimate.q", "std.err.q", "t.value", "p.value", "ci.lower", "ci.upper",
    "cb.lower", "cb.upper", "N.Co", "N.Tr"
  )
  missing <- setdiff(expected, names(dat))
  if (length(missing)) stop(sprintf("File %s is missing columns: %s", file, paste(missing, collapse = ", ")), call. = FALSE)
  dat[, expected]
}

read_simulation_rows <- function(m, dgp, design, method, tag, suffix, specs) {
  files <- file.path(output_dir, sprintf("simuls_raw_rep%04d_dgp%02d_%s_%s_%s_%s.csv", seq_len(m), dgp, design, method, tag, suffix))
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) stop(sprintf("Missing simulation file: %s", missing_files[1]), call. = FALSE)
  rows <- do.call(rbind, lapply(files, read_raw_file))
  expected_rows <- nrow(specs) * (21 + 2)
  rows_by_file <- table(rows$replication)
  if (any(rows_by_file != expected_rows)) {
    bad <- names(rows_by_file)[rows_by_file != expected_rows][1]
    stop(sprintf("Replication %s has %d rows; expected %d.", bad, rows_by_file[[bad]], expected_rows), call. = FALSE)
  }
  if (!all(rows$design == design)) stop(sprintf("File %s contains the wrong design.", files[1]), call. = FALSE)
  if (!setequal(unique(rows$estimand), specs$estimand)) stop(sprintf("File %s contains the wrong estimands.", files[1]), call. = FALSE)
  rows
}

## DGP targets
load_spp_data <- function(path = "spp.csv") {
  dat <- read.csv(path)
  vars <- c("running_saber11", "running_sisben", "spadies_any", "eligible_spp", "beneficiary_spp")
  missing <- setdiff(vars, names(dat))
  if (length(missing)) stop("Missing variables in spp.csv: ", paste(missing, collapse = ", "))
  dat <- dat[stats::complete.cases(dat[, vars]), vars]
  names(dat) <- c("x.1", "x.2", "y", "d", "w")
  dat$x.1 <- as.numeric(dat$x.1)
  dat$x.2 <- as.numeric(dat$x.2)
  dat$y <- as.numeric(dat$y)
  dat$d <- as.numeric(dat$d)
  dat$w <- as.numeric(dat$w)
  dat
}

fit_side_model <- function(response, data, degree = c("linear", "quadratic")) {
  degree <- match.arg(degree)
  rhs <- if (degree == "linear") {
    "x.1 + x.2"
  } else {
    "x.1 + x.2 + I(x.1^2) + I(x.2^2) + I(x.1 * x.2)"
  }
  f <- stats::as.formula(paste(response, "~", rhs))
  list(
    control = stats::lm(f, data = data[data$d == 0, ]),
    treated = stats::lm(f, data = data[data$d == 1, ])
  )
}

predict_side <- function(models, newdata, d_vec) {
  out <- numeric(nrow(newdata))
  idx0 <- d_vec == 0
  idx1 <- d_vec == 1
  if (any(idx0)) out[idx0] <- stats::predict(models$control, newdata = newdata[idx0, , drop = FALSE])
  if (any(idx1)) out[idx1] <- stats::predict(models$treated, newdata = newdata[idx1, , drop = FALSE])
  out
}

predict_prob <- function(models, newdata, d_vec) {
  p <- predict_side(models, newdata, d_vec)
  pmin(pmax(p, 1e-4), 1 - 1e-4)
}

build_targets <- function(boundary) {
  spp <- load_spp_data()
  mean_models_y <- list(
    linear = fit_side_model("y", spp, "linear"),
    quadratic = fit_side_model("y", spp, "quadratic")
  )
  mean_models_fuzzy <- list(
    linear = fit_side_model("w", spp, "linear"),
    quadratic = fit_side_model("w", spp, "quadratic")
  )

  dgp_grid <- expand.grid(
    mean = c("linear", "quadratic"),
    variance = c("homoskedastic", "heteroskedastic"),
    stringsAsFactors = FALSE
  )
  dgp_grid$tag <- paste(dgp_grid$mean, dgp_grid$variance, sep = "_")

  out <- vector("list", nrow(dgp_grid) * nrow(estimand_specs))
  k <- 1L
  for (g in seq_len(nrow(dgp_grid))) {
    mean_type <- dgp_grid$mean[g]
    mu0 <- predict_side(mean_models_y[[mean_type]], boundary, rep(0L, nrow(boundary)))
    mu1 <- predict_side(mean_models_y[[mean_type]], boundary, rep(1L, nrow(boundary)))
    pi0 <- predict_prob(mean_models_fuzzy[[mean_type]], boundary, rep(0L, nrow(boundary)))
    pi1 <- predict_prob(mean_models_fuzzy[[mean_type]], boundary, rep(1L, nrow(boundary)))
    tau_itt <- mu1 - mu0
    tau_fs <- pi1 - pi0
    tau_fuzzy <- ifelse(abs(tau_fs) < 1e-8, NA_real_, tau_itt / tau_fs)

    target_map <- list(fuzzy = tau_fuzzy, itt = tau_itt, fs = tau_fs, sharp = tau_itt, sharp0 = mu0)
    for (i in seq_len(nrow(estimand_specs))) {
      estimand <- estimand_specs$estimand[i]
      tau <- target_map[[estimand]]
      out[[k]] <- data.frame(
        design = estimand_specs$design[i],
        dgp = g,
        tag = dgp_grid$tag[g],
        estimand = estimand,
        row = c(as.character(seq_along(tau)), "WBATE", "LBATE"),
        target = c(tau, mean(tau, na.rm = TRUE), max(tau, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out[seq_len(k - 1L)])
}

## Summaries
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x) | !is.finite(x), "", sprintf(paste0("%.", digits, "f"), x))
}

coverage <- function(lower, target, upper) {
  is.finite(lower) & is.finite(upper) & lower <= target & target <= upper
}

bandwidth_summary <- function(z) {
  c(h0 = mean(z$h0, na.rm = TRUE), h1 = mean(z$h1, na.rm = TRUE))
}

effective_sample_summary <- function(z) {
  c(N.Co = mean(z$N.Co, na.rm = TRUE), N.Tr = mean(z$N.Tr, na.rm = TRUE))
}

format_boundary_label <- function(row) {
  if (row %in% c("WBATE", "LBATE", "Uniform")) {
    return(sprintf("$\\mathtt{%s}$", row))
  }
  sprintf("$\\bb_{%d}$", as.integer(row))
}

summarize_point_rows <- function(rows, targets) {
  point_rows <- rows[!(rows$row %in% c("WBATE", "LBATE")), ]
  split_rows <- split(point_rows, point_rows$row)
  out <- do.call(rbind, lapply(split_rows, function(z) {
    row <- z$row[1]
    truth <- targets$target[targets$row == row][1]
    bw <- bandwidth_summary(z)
    ns <- effective_sample_summary(z)
    data.frame(
      row = row,
      h0 = bw["h0"],
      h1 = bw["h1"],
      N.Co = ns["N.Co"],
      N.Tr = ns["N.Tr"],
      Bias = mean(z$estimate.p, na.rm = TRUE) - truth,
      SD = stats::sd(z$estimate.p, na.rm = TRUE),
      RMSE = sqrt(mean((z$estimate.p - truth)^2, na.rm = TRUE)),
      EC = mean(coverage(z$ci.lower, truth, z$ci.upper), na.rm = TRUE),
      IL = mean(z$ci.upper - z$ci.lower, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  out[order(as.integer(out$row)), ]
}

summarize_uniform_row <- function(rows, targets) {
  point_rows <- rows[!(rows$row %in% c("WBATE", "LBATE")), ]
  uniform_coverage <- vapply(split(point_rows, point_rows$replication), function(z) {
    z <- z[order(as.integer(z$row)), ]
    truth <- targets$target[match(z$row, targets$row)]
    all(coverage(z$cb.lower, truth, z$cb.upper))
  }, logical(1))
  data.frame(
    row = "Uniform",
    h0 = NA_real_,
    h1 = NA_real_,
    N.Co = NA_real_,
    N.Tr = NA_real_,
    Bias = NA_real_,
    SD = NA_real_,
    RMSE = NA_real_,
    EC = mean(uniform_coverage, na.rm = TRUE),
    IL = mean(point_rows$cb.upper - point_rows$cb.lower, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summarize_aggregate_rows <- function(rows, targets) {
  aggregate_rows <- rows[rows$row %in% c("WBATE", "LBATE"), ]
  do.call(rbind, lapply(c("WBATE", "LBATE"), function(name) {
    z <- aggregate_rows[aggregate_rows$row == name, ]
    if (!nrow(z)) stop(sprintf("Missing aggregate row: %s", name), call. = FALSE)
    truth <- targets$target[targets$row == name][1]
    data.frame(
      row = name,
      h0 = NA_real_,
      h1 = NA_real_,
      N.Co = NA_real_,
      N.Tr = NA_real_,
      Bias = mean(z$estimate.p, na.rm = TRUE) - truth,
      SD = stats::sd(z$estimate.p, na.rm = TRUE),
      RMSE = sqrt(mean((z$estimate.p - truth)^2, na.rm = TRUE)),
      EC = mean(coverage(z$ci.lower, truth, z$ci.upper), na.rm = TRUE),
      IL = mean(z$ci.upper - z$ci.lower, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

simulation_summary <- function(rows, targets, estimand) {
  rows_now <- rows[rows$estimand == estimand, ]
  targets_now <- targets[targets$estimand == estimand, ]
  rbind(
    summarize_point_rows(rows_now, targets_now),
    summarize_uniform_row(rows_now, targets_now),
    summarize_aggregate_rows(rows_now, targets_now)
  )
}

## Output tables
write_simulation_table <- function(summary_rows, file) {
  lines <- c(
    "\\begin{tabular}{lccccrrrrr}",
    "\\toprule\\toprule",
    " & $h_0$ & $h_1$ & $N_{\\mathrm{Co}}$ & $N_{\\mathrm{Tr}}$ & Bias & SD & RMSE & EC & IL \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(summary_rows))) {
    row <- summary_rows[i, ]
    if (row$row %in% c("Uniform", "WBATE", "LBATE")) lines <- c(lines, "\\midrule")
    lines <- c(lines, paste0(
      format_boundary_label(row$row), " & ",
      fmt_num(row$h0), " & ",
      fmt_num(row$h1), " & ",
      fmt_num(row$N.Co, 0), " & ",
      fmt_num(row$N.Tr, 0), " & ",
      fmt_num(row$Bias), " & ",
      fmt_num(row$SD), " & ",
      fmt_num(row$RMSE), " & ",
      fmt_num(row$EC), " & ",
      fmt_num(row$IL), " \\\\")
    )
  }
  lines <- c(lines, "\\bottomrule\\bottomrule", "\\end{tabular}")
  writeLines(lines, file)
  invisible(lines)
}

write_dgp_target_table <- function(targets, estimand, file) {
  rows <- targets[targets$estimand == estimand & !(targets$row %in% c("WBATE", "LBATE")), ]
  rows$row_int <- as.integer(rows$row)
  tags <- unique(rows$tag)
  lines <- c(
    "\\begin{tabular}{lrrrr}",
    "\\toprule\\toprule",
    paste(c("", sprintf("%s", tags)), collapse = " & "),
    "\\\\",
    "\\midrule"
  )
  for (i in sort(unique(rows$row_int))) {
    vals <- rows$target[rows$row_int == i]
    lines <- c(lines, paste0("$\\bb_{", i, "}$ & ", paste(fmt_num(vals), collapse = " & "), " \\\\"))
  }
  lines <- c(lines, "\\bottomrule\\bottomrule", "\\end{tabular}")
  writeLines(lines, file)
}

fmt_sci <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || abs(x) < 1e-12) return("$0$")
  parts <- strsplit(formatC(x, format = "e", digits = 2), "e", fixed = TRUE)[[1]]
  coef <- parts[1]
  exponent <- as.integer(parts[2])
  if (exponent == 0L) return(paste0("$", coef, "$"))
  paste0("$", coef, " \\times 10^{", exponent, "}$")
}

parameter_value <- function(parameters, design, component, term, mean, variance, side,
                            default = NULL) {
  rows <- parameters[
    parameters$design == design &
      parameters$component == component &
      parameters$term == term &
      parameters$mean == mean &
      parameters$variance == variance &
      parameters$side == side,
  ]
  if (nrow(rows) == 0L && !is.null(default)) return(default)
  if (nrow(rows) != 1L) {
    stop(
      sprintf(
        "Expected one DGP parameter for %s/%s/%s/%s/%s/t=%s; found %d.",
        design, component, term, mean, variance, side, nrow(rows)
      ),
      call. = FALSE
    )
  }
  rows$value
}

write_dgp_mean_table <- function(parameters, component, labels, file) {
  terms <- c("(Intercept)", "x.1", "x.2", "I(x.1^2)", "I(x.1 * x.2)", "I(x.2^2)")
  lines <- c(
    "\\begin{tabular}{@{}rrrrr@{}}",
    "\\hline\\hline",
    "    & \\multicolumn{2}{c}{Linear} & \\multicolumn{2}{c}{Quadratic} \\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
    "    & \\multicolumn{1}{c}{$t=0$} & \\multicolumn{1}{c}{$t=1$} & \\multicolumn{1}{c}{$t=0$} & \\multicolumn{1}{c}{$t=1$} \\\\ \\midrule"
  )

  for (i in seq_along(terms)) {
    vals <- c(
      parameter_value(parameters, "fuzzy", component, terms[i], "linear", "homoskedastic", 0, default = 0),
      parameter_value(parameters, "fuzzy", component, terms[i], "linear", "homoskedastic", 1, default = 0),
      parameter_value(parameters, "fuzzy", component, terms[i], "quadratic", "homoskedastic", 0, default = 0),
      parameter_value(parameters, "fuzzy", component, terms[i], "quadratic", "homoskedastic", 1, default = 0)
    )
    lines <- c(lines, paste0(labels[i], " & ", paste(vapply(vals, fmt_sci, character(1)), collapse = " & "), " \\\\"))
  }

  if (identical(component, "beta_w")) {
    vals <- c(
      parameter_value(parameters, "fuzzy", "lambda", "", "linear", "homoskedastic", 0),
      parameter_value(parameters, "fuzzy", "lambda", "", "linear", "homoskedastic", 1),
      parameter_value(parameters, "fuzzy", "lambda", "", "quadratic", "homoskedastic", 0),
      parameter_value(parameters, "fuzzy", "lambda", "", "quadratic", "homoskedastic", 1)
    )
    lines <- c(lines, paste0("$\\lambda_t$ & ", paste(vapply(vals, fmt_sci, character(1)), collapse = " & "), " \\\\"))
  }

  lines <- c(lines, "\\hline\\hline", "\\end{tabular}")
  writeLines(lines, file)
}

write_dgp_variance_table <- function(parameters, file) {
  terms <- c("(Intercept)", "x.1", "x.2", "I(x.1^2)", "I(x.1 * x.2)", "I(x.2^2)")
  labels <- c(
    "$\\delta_{t,0}$", "$\\delta_{t,11}$", "$\\delta_{t,12}$",
    "$\\delta_{t,21}$", "$\\delta_{t,22}$", "$\\delta_{t,23}$"
  )
  lines <- c(
    "\\begin{tabular}{@{}rrrrr@{}}",
    "\\hline\\hline",
    "    & \\multicolumn{2}{c}{Homoskedastic} & \\multicolumn{2}{c}{Heteroskedastic} \\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
    "    & \\multicolumn{1}{c}{$t=0$} & \\multicolumn{1}{c}{$t=1$} & \\multicolumn{1}{c}{$t=0$} & \\multicolumn{1}{c}{$t=1$} \\\\ \\midrule"
  )

  for (i in seq_along(terms)) {
    homo_vals <- if (identical(terms[i], "(Intercept)")) {
      c(
        parameter_value(parameters, "fuzzy", "delta0", "(Intercept)", "linear", "homoskedastic", 0),
        parameter_value(parameters, "fuzzy", "delta0", "(Intercept)", "linear", "homoskedastic", 1)
      )
    } else {
      c(0, 0)
    }
    hetero_vals <- c(
      parameter_value(parameters, "fuzzy", "delta2", terms[i], "quadratic", "heteroskedastic", 0),
      parameter_value(parameters, "fuzzy", "delta2", terms[i], "quadratic", "heteroskedastic", 1)
    )
    vals <- c(homo_vals, hetero_vals)
    lines <- c(lines, paste0(labels[i], " & ", paste(vapply(vals, fmt_sci, character(1)), collapse = " & "), " \\\\"))
  }

  lines <- c(lines, "\\hline\\hline", "\\end{tabular}")
  writeLines(lines, file)
}

write_dgp_parameter_tables <- function(parameters, path) {
  beta_y_labels <- c(
    "$\\beta_{Y,t,0}$", "$\\beta_{Y,t,11}$", "$\\beta_{Y,t,12}$",
    "$\\beta_{Y,t,21}$", "$\\beta_{Y,t,22}$", "$\\beta_{Y,t,23}$"
  )
  beta_w_labels <- c(
    "$\\beta_{W,t,0}$", "$\\beta_{W,t,11}$", "$\\beta_{W,t,12}$",
    "$\\beta_{W,t,21}$", "$\\beta_{W,t,22}$", "$\\beta_{W,t,23}$"
  )
  files <- file.path(
    path,
    c("simuls_dgp_itt_mean.tex", "simuls_dgp_itt_var.tex", "simuls_dgp_fs_mean.tex")
  )
  write_dgp_mean_table(parameters, "beta_y", beta_y_labels, files[1])
  write_dgp_variance_table(parameters, files[2])
  write_dgp_mean_table(parameters, "beta_w", beta_w_labels, files[3])
  files
}

check_generated_tables <- function(expected) {
  missing <- expected[!file.exists(expected)]
  if (length(missing)) stop(sprintf("Missing expected simulation table: %s", missing[1]), call. = FALSE)
  for (file in expected) {
    txt <- paste(readLines(file, warn = FALSE), collapse = "\n")
    needed <- c("Uniform", "WBATE", "LBATE", "$N_{\\mathrm{Co}}$", "$N_{\\mathrm{Tr}}$")
    for (needle in needed) {
      if (!grepl(needle, txt, fixed = TRUE)) stop(sprintf("Generated table %s is missing '%s'.", file, needle), call. = FALSE)
    }
  }
  invisible(TRUE)
}

## Run
clean_tables_dir(tables_dir)

config <- read_output_csv("simuls_config.csv")
m <- as.integer(metadata_value(config, "m"))
boundary <- read_output_csv("simuls_eval_grid.csv")
targets <- read_output_csv("simuls_dgp_targets.csv")
parameters <- read_output_csv("simuls_dgp_parameters.csv")

expected <- character()
for (design in unique(estimand_specs$design)) {
  specs <- estimand_specs[estimand_specs$design == design, ]
  suffixes <- if (identical(design, "sharp")) "bwmain" else bandwidth_specs$suffix
  for (suffix in suffixes) {
    for (tag in unique(targets$tag[targets$design == design])) {
      dgp <- unique(targets$dgp[targets$design == design & targets$tag == tag])
      if (length(dgp) != 1L) stop(sprintf("Could not identify a unique DGP for %s / %s.", design, tag), call. = FALSE)
      targets_now <- targets[targets$design == design & targets$tag == tag, ]
      for (method in methods) {
        rows <- read_simulation_rows(m, dgp, design, method, tag, suffix, specs)
        for (estimand in specs$estimand) {
          summary_rows <- simulation_summary(rows, targets_now, estimand)
          file <- file.path(tables_dir, sprintf("simuls_%s_%s_%s_%s.tex", estimand, tag, method, suffix))
          write_simulation_table(summary_rows, file)
          expected <- c(expected, file)
          cat(sprintf("Wrote %s (%s, %s, %s).\n", basename(file), method_labels[[method]], estimand_labels[[estimand]], suffix))
        }
      }
    }
  }
}

for (estimand in estimands) {
  write_dgp_target_table(targets, estimand, file.path(tables_dir, sprintf("simuls_dgp_%s.tex", estimand)))
}
parameter_files <- write_dgp_parameter_tables(parameters, tables_dir)
cat(sprintf("Wrote %d DGP parameter table file(s).\n", length(parameter_files)))

check_generated_tables(expected)
cat(sprintf("Simulation LaTeX tables complete. Generated %d table file(s).\n", length(expected)))
