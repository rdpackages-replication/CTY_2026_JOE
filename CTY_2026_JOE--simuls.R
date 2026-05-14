## Replication Script: Simulations
## Cattaneo, Titiunik and Yu (2026): Boundary Distance RD

## Plot options
options(scipen = 999, digits = 4)

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/"))
  }
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
setwd(script_dir)

## Dependencies
library(parallel)
library(foreach)
library(doParallel)
library(rdrobust)

## rd2d package
library(rd2d)

## Directories
output_dir <- Sys.getenv("RD2D_OUTPUT_DIR", unset = "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

clean_output_dir <- function(path) {
  files <- list.files(path, pattern = "^simuls_.*[.]csv$", full.names = TRUE)
  if (length(files)) file.remove(files)
}
invisible(clean_output_dir(output_dir))

write_output <- function(x, file) {
  utils::write.csv(x, file = file, row.names = FALSE, na = "")
}

get_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || identical(value, "")) return(as.integer(default))
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0) {
    stop(sprintf("%s must be a positive integer.", name), call. = FALSE)
  }
  parsed
}

simulation_seed <- function(seed_base, replication, design, dgp) {
  design_offset <- if (identical(design, "sharp")) 500000L else 0L
  seed_base + replication * 1000L + design_offset + dgp
}

rd2d_package_version <- function() {
  desc_path <- system.file("DESCRIPTION", package = "rd2d")
  if (!nzchar(desc_path)) stop("Could not find installed rd2d DESCRIPTION file.", call. = FALSE)
  as.character(read.dcf(desc_path, fields = "Version")[1, 1])
}

## Simulation controls
set.seed(20260510)

n <- get_env_int("RD2D_N", 20000)
m <- get_env_int("RD2D_M", 5000)
repp <- get_env_int("RD2D_REPP", 2000)
seed_base <- get_env_int("RD2D_SEED", 20260510)

num_cores <- parallel::detectCores()
if (is.na(num_cores)) num_cores <- 2L
default_workers <- max(1L, num_cores - 4L)
workers <- get_env_int("RD2D_WORKERS", default_workers)
workers <- min(m, max(1L, min(workers, num_cores)))

bandwidth_specs <- data.frame(
  bwparam = c("main", "itt"),
  suffix = c("bwmain", "bwitt"),
  stringsAsFactors = FALSE
)

estimand_specs <- data.frame(
  design = c("fuzzy", "fuzzy", "fuzzy", "sharp", "sharp"),
  output = c("main", "itt", "fs", "main", "main.0"),
  estimand = c("fuzzy", "itt", "fs", "sharp", "sharp0"),
  stringsAsFactors = FALSE
)
fuzzy_estimand_specs <- estimand_specs[estimand_specs$design == "fuzzy", ]
sharp_estimand_specs <- estimand_specs[estimand_specs$design == "sharp", ]

## Data and boundary evaluation grid
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

  expected_assignment <- as.integer(dat$x.1 >= 0 & dat$x.2 >= 0)
  if (!all(dat$d == expected_assignment)) {
    stop("eligible_spp does not match the quadrant assignment rule.", call. = FALSE)
  }

  dat
}

make_eval_grid <- function() {
  neval <- 40L
  half <- ceiling(neval / 2)
  grid <- rbind(
    data.frame(
      x.1 = rep(0, half),
      x.2 = 40 - (seq_len(half) - 1) * 40 / half
    ),
    data.frame(
      x.1 = (seq_len(neval - half) - 1) * 56 / half,
      x.2 = rep(0, neval - half)
    )
  )
  grid[11:31, , drop = FALSE]
}

spp <- load_spp_data()
boundary <- make_eval_grid()
kink.position <- 11L

## Data generating processes
design_matrix <- function(dat, s) {
  x1 <- dat$x.1
  x2 <- dat$x.2

  if (s == 0) {
    return(cbind("(Intercept)" = 1))
  }
  if (s == 1) {
    return(cbind("(Intercept)" = 1, "x.1" = x1, "x.2" = x2))
  }
  if (s == 2) {
    return(cbind(
      "(Intercept)" = 1,
      "x.1" = x1,
      "x.2" = x2,
      "I(x.1^2)" = x1^2,
      "I(x.1 * x.2)" = x1 * x2,
      "I(x.2^2)" = x2^2
    ))
  }

  stop("Only s = 0, s = 1, and s = 2 are supported.", call. = FALSE)
}

as_finite_coef <- function(fit, names, what) {
  coef <- as.numeric(fit$coefficients)
  names(coef) <- names
  if (any(!is.finite(coef))) stop("Non-finite calibration coefficient in ", what, ".", call. = FALSE)
  coef
}

calibrate_dgp <- function(data, mean_type, variance_type, dgp,
                          design = c("fuzzy", "sharp")) {
  design <- match.arg(design)
  s <- switch(mean_type, linear = 1L, quadratic = 2L)
  v <- switch(variance_type, homoskedastic = 0L, heteroskedastic = 2L)
  Xs <- design_matrix(data, s)
  X2 <- design_matrix(data, 2L)

  out <- list(
    design = design,
    dgp = dgp,
    mean = mean_type,
    variance = variance_type,
    tag = paste(mean_type, variance_type, sep = "_"),
    s = s,
    v = v
  )

  for (side in 0:1) {
    ind <- data$d == side
    X_side <- Xs[ind, , drop = FALSE]
    X2_side <- X2[ind, , drop = FALSE]

    y_fit <- stats::lm.fit(X_side, data$y[ind])
    beta_y <- as_finite_coef(y_fit, colnames(X_side), sprintf("outcome, side %d", side))

    w_fit <- suppressWarnings(stats::glm.fit(X_side, data$w[ind], family = stats::binomial()))
    beta_w <- as_finite_coef(w_fit, colnames(X_side), sprintf("first stage, side %d", side))

    mu_y <- as.numeric(X_side %*% beta_y)
    mu_w <- stats::plogis(as.numeric(X_side %*% beta_w))
    eps_y <- data$y[ind] - mu_y
    eps_w <- data$w[ind] - mu_w

    denom <- sum(eps_w^2)
    if (!is.finite(denom) || denom <= 0) {
      stop(sprintf("Degenerate first-stage residual variance for s=%d, side=%d.", s, side), call. = FALSE)
    }

    lambda <- if (identical(design, "sharp")) 0 else sum(eps_y * eps_w) / denom
    u_y <- eps_y - lambda * eps_w
    shock_var <- mean(u_y^2)
    if (!is.finite(shock_var) || shock_var <= 0) {
      stop(sprintf("Degenerate outcome shock variance for s=%d, side=%d.", s, side), call. = FALSE)
    }

    delta0 <- log(shock_var)
    delta2_fit <- stats::lm.fit(X2_side, log(pmax(u_y^2, 1e-8)))
    delta2 <- as_finite_coef(delta2_fit, colnames(X2_side), sprintf("log variance, side %d", side))

    out[[paste0("beta_y_", side)]] <- beta_y
    out[[paste0("beta_w_", side)]] <- beta_w
    out[[paste0("lambda_", side)]] <- lambda
    out[[paste0("delta0_", side)]] <- delta0
    out[[paste0("delta2_", side)]] <- delta2
  }

  out
}

predict_linear <- function(coef, x, s) {
  as.numeric(design_matrix(x, s) %*% coef)
}

predict_probability <- function(coef, x, s) {
  stats::plogis(predict_linear(coef, x, s))
}

predict_sigma <- function(cal, x, side) {
  if (cal$v == 0L) {
    log_sigma2 <- rep(cal[[paste0("delta0_", side)]], nrow(x))
  } else {
    log_sigma2 <- predict_linear(cal[[paste0("delta2_", side)]], x, 2L)
  }
  sqrt(exp(log_sigma2))
}

make_targets <- function(eval, calibrations) {
  out <- vector("list", length(calibrations) * nrow(estimand_specs))
  k <- 1L

  for (cal in calibrations) {
    mu_y0 <- predict_linear(cal$beta_y_0, eval, cal$s)
    mu_y1 <- predict_linear(cal$beta_y_1, eval, cal$s)
    mu_w0 <- predict_probability(cal$beta_w_0, eval, cal$s)
    mu_w1 <- predict_probability(cal$beta_w_1, eval, cal$s)

    tau_itt <- mu_y1 - mu_y0
    tau_fs <- mu_w1 - mu_w0
    tau_fuzzy <- ifelse(abs(tau_fs) < 1e-8, NA_real_, tau_itt / tau_fs)
    target_map <- list(
      fuzzy = tau_fuzzy,
      itt = tau_itt,
      fs = tau_fs,
      sharp = tau_itt,
      sharp0 = mu_y0
    )

    specs <- estimand_specs[estimand_specs$design == cal$design, ]
    for (estimand in specs$estimand) {
      tau <- target_map[[estimand]]
      out[[k]] <- data.frame(
        design = cal$design,
        dgp = cal$dgp,
        tag = cal$tag,
        mean = cal$mean,
        variance = cal$variance,
        s = cal$s,
        v = cal$v,
        estimand = estimand,
        row = c(as.character(seq_along(tau)), "WBATE", "LBATE"),
        index = c(seq_along(tau), NA_integer_, NA_integer_),
        aggregate = c(rep("", length(tau)), "WBATE", "LBATE"),
        target = c(tau, mean(tau, na.rm = TRUE), max(tau, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  do.call(rbind, out[seq_len(k - 1L)])
}

make_dgp_parameter_table <- function(calibrations) {
  rows <- lapply(calibrations, function(cal) {
    coef_rows <- do.call(rbind, lapply(0:1, function(side) {
      do.call(rbind, lapply(c("beta_y", "beta_w", "delta2"), function(component) {
        values <- cal[[paste0(component, "_", side)]]
        data.frame(
          design = cal$design,
          dgp = cal$dgp,
          tag = cal$tag,
          mean = cal$mean,
          variance = cal$variance,
          s = cal$s,
          v = cal$v,
          side = side,
          component = component,
          term = names(values),
          value = as.numeric(values),
          stringsAsFactors = FALSE
        )
      }))
    }))

    scalar_rows <- do.call(rbind, lapply(0:1, function(side) {
      data.frame(
        design = cal$design,
        dgp = cal$dgp,
        tag = cal$tag,
        mean = cal$mean,
        variance = cal$variance,
        s = cal$s,
        v = cal$v,
        side = side,
        component = c("lambda", "delta0"),
        term = c("", "(Intercept)"),
        value = c(cal[[paste0("lambda_", side)]], cal[[paste0("delta0_", side)]]),
        stringsAsFactors = FALSE
      )
    }))

    rbind(coef_rows, scalar_rows)
  })

  do.call(rbind, rows)
}

dgp_grid <- expand.grid(
  mean = c("linear", "quadratic"),
  variance = c("homoskedastic", "heteroskedastic"),
  stringsAsFactors = FALSE
)
dgp_grid$dgp <- seq_len(nrow(dgp_grid))
dgp_grid$s <- ifelse(dgp_grid$mean == "linear", 1L, 2L)
dgp_grid$v <- ifelse(dgp_grid$variance == "homoskedastic", 0L, 2L)
dgp_grid$tag <- paste(dgp_grid$mean, dgp_grid$variance, sep = "_")

fuzzy_calibrations <- lapply(seq_len(nrow(dgp_grid)), function(g) {
  calibrate_dgp(spp, dgp_grid$mean[g], dgp_grid$variance[g], dgp_grid$dgp[g], design = "fuzzy")
})
sharp_calibrations <- lapply(seq_len(nrow(dgp_grid)), function(g) {
  calibrate_dgp(spp, dgp_grid$mean[g], dgp_grid$variance[g], dgp_grid$dgp[g], design = "sharp")
})
calibrations <- c(fuzzy_calibrations, sharp_calibrations)

draw_scores <- function(n) {
  data.frame(
    x.1 = 100 * stats::rbeta(n, 3, 4) - 25,
    x.2 = 100 * stats::rbeta(n, 3, 4) - 25
  )
}

assignment_rule <- function(x) {
  as.integer(x$x.1 >= 0 & x$x.2 >= 0)
}

simulate_sample <- function(cal, n) {
  x <- draw_scores(n)
  assignment <- assignment_rule(x)

  mu_y0 <- predict_linear(cal$beta_y_0, x, cal$s)
  mu_y1 <- predict_linear(cal$beta_y_1, x, cal$s)

  if (identical(cal$design, "fuzzy")) {
    mu_w0 <- predict_probability(cal$beta_w_0, x, cal$s)
    mu_w1 <- predict_probability(cal$beta_w_1, x, cal$s)
    fuzzy0 <- as.numeric(stats::runif(n) <= mu_w0)
    fuzzy1 <- as.numeric(stats::runif(n) <= mu_w1)
    eps_y0 <- stats::rnorm(n, sd = predict_sigma(cal, x, 0L))
    eps_y1 <- stats::rnorm(n, sd = predict_sigma(cal, x, 1L))
    y0 <- mu_y0 + cal$lambda_0 * (fuzzy0 - mu_w0) + eps_y0
    y1 <- mu_y1 + cal$lambda_1 * (fuzzy1 - mu_w1) + eps_y1
    fuzzy <- ifelse(assignment == 1L, fuzzy1, fuzzy0)
  } else {
    eps_y0 <- stats::rnorm(n, sd = predict_sigma(cal, x, 0L))
    eps_y1 <- stats::rnorm(n, sd = predict_sigma(cal, x, 1L))
    y0 <- mu_y0 + eps_y0
    y1 <- mu_y1 + eps_y1
    fuzzy <- NULL
  }

  list(
    x = x,
    assignment = assignment,
    y = ifelse(assignment == 1L, y1, y0),
    fuzzy = fuzzy
  )
}

## Helper functions
make_distances <- function(x, eval, assignment) {
  distance <- as.matrix(proxy::dist(x, eval, method = "Euclidean"))
  distance <- distance * matrix(2 * assignment - 1, nrow = nrow(distance), ncol = ncol(distance))
  colnames(distance) <- sprintf("b%02d", seq_len(ncol(distance)))
  distance
}

rdrobust_bw <- function(outcome, running) {
  out <- tryCatch(rdrobust::rdbwselect(outcome, running, vce = "hc1"), error = function(e) NULL)
  if (is.null(out) || is.null(out$bws)) return(c(NA_real_, NA_real_))
  bws <- as.numeric(out$bws[1:2])
  if (length(bws) < 2 || any(!is.finite(bws)) || any(bws <= 0)) return(c(NA_real_, NA_real_))
  bws
}

rdrobust_tau <- function(outcome, running) {
  out <- tryCatch(rdrobust::rdrobust(outcome, running, vce = "hc1"), error = function(e) NULL)
  if (is.null(out) || is.null(out$coef)) return(NA_real_)
  tau <- suppressWarnings(as.numeric(out$coef[1]))
  if (!is.finite(tau)) NA_real_ else tau
}

linearized_fuzzy_outcome <- function(y, fuzzy, running) {
  tau_itt <- rdrobust_tau(y, running)
  tau_fs <- rdrobust_tau(fuzzy, running)
  if (!is.finite(tau_itt) || !is.finite(tau_fs) || abs(tau_fs) < 1e-8) return(y)
  y / tau_fs - tau_itt * fuzzy / (tau_fs^2)
}

compute_rdrobust_bws <- function(y, fuzzy = NULL, distance, bwparam = "main") {
  bws <- matrix(NA_real_, nrow = ncol(distance), ncol = 2)
  for (i in seq_len(ncol(distance))) {
    outcome_bw <- if (!is.null(fuzzy) && identical(bwparam, "main")) {
      linearized_fuzzy_outcome(y, fuzzy, distance[, i])
    } else {
      y
    }
    bws[i, ] <- rdrobust_bw(outcome_bw, distance[, i])
    if (any(!is.finite(bws[i, ]))) bws[i, ] <- rdrobust_bw(y, distance[, i])
    if (any(!is.finite(bws[i, ]))) {
      fallback <- stats::median(abs(distance[, i]), na.rm = TRUE)
      bws[i, ] <- c(fallback, fallback)
    }
  }
  colnames(bws) <- c("h0", "h1")
  bws
}

summary_table <- function(result, output) {
  s <- NULL
  utils::capture.output({
    s <- summary(
      result,
      output = output,
      cbands = output,
      WBATE = rep(1, nrow(result[[output]])),
      LBATE = TRUE
    )
  })
  tab <- s$tables[[output]]
  tab <- data.frame(row = rownames(tab), tab, row.names = NULL, check.names = FALSE)
  tab
}

make_output_rows <- function(result, design, method, replication, dgp, tag, bwparam, suffix, specs) {
  pieces <- vector("list", nrow(specs))
  for (i in seq_len(nrow(specs))) {
    output <- specs$output[i]
    estimand <- specs$estimand[i]
    tab <- summary_table(result, output)
    pieces[[i]] <- data.frame(
      replication = replication,
      dgp = dgp,
      tag = tag,
      method = method,
      bwparam = bwparam,
      suffix = suffix,
      design = design,
      estimand = estimand,
      tab,
      row.names = NULL,
      check.names = FALSE
    )
  }
  do.call(rbind, pieces)
}

raw_result_columns <- c(
  "replication", "dgp", "tag", "method", "bwparam", "suffix", "design", "estimand", "row",
  "h0", "h1", "h0.rbc", "h1.rbc", "estimate.p", "std.err.p",
  "estimate.q", "std.err.q", "t.value", "p.value", "ci.lower", "ci.upper",
  "cb.lower", "cb.upper", "N.Co", "N.Tr"
)

validate_simulation_rows <- function(x, design, method, replication, dgp, tag, bwparam, suffix, specs) {
  expected <- raw_result_columns
  missing <- setdiff(expected, names(x))
  if (length(missing)) stop("Missing columns for ", method, ": ", paste(missing, collapse = ", "))
  if (!all(x$replication == replication)) stop("Replication metadata mismatch for ", method)
  if (!all(x$dgp == dgp)) stop("DGP metadata mismatch for ", method)
  if (!all(x$tag == tag)) stop("Tag metadata mismatch for ", method)
  if (!all(x$method == method)) stop("Method metadata mismatch for ", method)
  if (!all(x$bwparam == bwparam)) stop("bwparam metadata mismatch for ", method)
  if (!all(x$suffix == suffix)) stop("suffix metadata mismatch for ", method)
  if (!all(x$design == design)) stop("Design metadata mismatch for ", method)
  if (!setequal(unique(x$estimand), specs$estimand)) stop("Estimand rows missing for ", method)
  for (estimand in specs$estimand) {
    xi <- x[x$estimand == estimand, ]
    if (!all(c("WBATE", "LBATE") %in% xi$row)) stop("Aggregate rows missing for ", method, " / ", estimand)
    point <- xi[!(xi$row %in% c("WBATE", "LBATE")), ]
    finite_cols <- if (identical(estimand, "sharp0")) {
      c("h0", "estimate.p", "std.err.p", "estimate.q", "std.err.q", "t.value", "p.value", "N.Co")
    } else {
      c("h0", "h1", "estimate.p", "std.err.p", "estimate.q", "std.err.q", "t.value", "p.value", "N.Co", "N.Tr")
    }
    bad <- finite_cols[!vapply(point[finite_cols], function(z) all(is.finite(z)), logical(1))]
    if (length(bad)) stop("Non-finite pointwise columns for ", method, " / ", estimand, ": ", paste(bad, collapse = ", "))
  }
  invisible(TRUE)
}

raw_replication_file <- function(replication) {
  file.path(output_dir, sprintf("simuls_raw_rep%04d.csv", replication))
}

raw_result_file <- function(dgp, design, method, tag, suffix) {
  file.path(output_dir, sprintf("simuls_raw_dgp%02d_%s_%s_%s_%s.csv", as.integer(dgp), design, method, tag, suffix))
}

raw_result_key <- function(dgp, design, method, tag, suffix) {
  paste(as.integer(dgp), design, method, tag, suffix, sep = "\t")
}

raw_result_rows <- function(result, design, method, replication, dgp, tag, bwparam, suffix, specs) {
  rows <- make_output_rows(result, design, method, replication, dgp, tag, bwparam, suffix, specs)
  validate_simulation_rows(rows, design, method, replication, dgp, tag, bwparam, suffix, specs)
  rows[, raw_result_columns]
}

make_raw_result_specs <- function() {
  methods <- c("kinkoff", "kinkon", "adaptive", "rdrobustadj")

  fuzzy <- do.call(rbind, lapply(fuzzy_calibrations, function(cal) {
    expand.grid(
      dgp = cal$dgp,
      design = cal$design,
      method = methods,
      tag = cal$tag,
      suffix = bandwidth_specs$suffix,
      stringsAsFactors = FALSE
    )
  }))

  sharp <- do.call(rbind, lapply(sharp_calibrations, function(cal) {
    expand.grid(
      dgp = cal$dgp,
      design = cal$design,
      method = methods,
      tag = cal$tag,
      suffix = "bwmain",
      stringsAsFactors = FALSE
    )
  }))

  unique(rbind(fuzzy, sharp))
}

close_raw_connections <- function(connections) {
  for (con in connections) {
    try(if (isOpen(con)) close(con), silent = TRUE)
  }
}

consolidate_raw_results <- function(replication_files, raw_specs) {
  final_files <- raw_result_file(raw_specs$dgp, raw_specs$design, raw_specs$method, raw_specs$tag, raw_specs$suffix)
  if (length(final_files)) invisible(file.remove(final_files[file.exists(final_files)]))

  keys <- raw_result_key(raw_specs$dgp, raw_specs$design, raw_specs$method, raw_specs$tag, raw_specs$suffix)
  connections <- setNames(lapply(final_files, file, open = "wt"), keys)
  on.exit(close_raw_connections(connections), add = TRUE)

  written <- setNames(rep(FALSE, length(keys)), keys)

  for (file in replication_files) {
    raw <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    missing <- setdiff(raw_result_columns, names(raw))
    if (length(missing)) {
      stop(sprintf("Scratch simulation file %s is missing columns: %s", file, paste(missing, collapse = ", ")), call. = FALSE)
    }
    raw <- raw[, raw_result_columns]

    row_keys <- raw_result_key(raw$dgp, raw$design, raw$method, raw$tag, raw$suffix)
    parts <- split(raw, row_keys, drop = TRUE)

    for (key in names(parts)) {
      if (!(key %in% names(connections))) {
        stop(sprintf("Unexpected simulation output cell in %s: %s", file, key), call. = FALSE)
      }
      utils::write.table(
        parts[[key]],
        file = connections[[key]],
        sep = ",",
        row.names = FALSE,
        col.names = !written[[key]],
        na = "",
        qmethod = "double"
      )
      written[[key]] <- TRUE
    }
  }

  missing_cells <- names(written)[!written]
  if (length(missing_cells)) {
    stop(sprintf("No rows were written for expected simulation cell: %s", missing_cells[1]), call. = FALSE)
  }

  close_raw_connections(connections)
  invisible(file.remove(replication_files))
  invisible(final_files)
}

## Metadata
write_output(data.frame(n = n, m = m, repp = repp, workers = workers), file.path(output_dir, "simuls_config.csv"))
write_output(dgp_grid, file.path(output_dir, "simuls_dgps.csv"))
write_output(make_dgp_parameter_table(calibrations), file.path(output_dir, "simuls_dgp_parameters.csv"))
write_output(make_targets(boundary, calibrations), file.path(output_dir, "simuls_dgp_targets.csv"))
write_output(boundary, file.path(output_dir, "simuls_eval_grid.csv"))
write_output(bandwidth_specs, file.path(output_dir, "simuls_bandwidths.csv"))
write_output(estimand_specs, file.path(output_dir, "simuls_estimands.csv"))
write_output(
  data.frame(
    name = c("rd2d.version", "n", "m", "repp", "seed_base", "workers", "generated_at"),
    value = c(rd2d_package_version(), n, m, repp, seed_base, workers, as.character(Sys.time())),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "simuls_metadata.csv")
)

## Run simulations
message(sprintf("Running %d replications with n = %d, repp = %d, workers = %d", m, n, repp, workers))

cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
doParallel::registerDoParallel(cl)
invisible(parallel::clusterExport(cl, "script_dir", envir = globalenv()))
invisible(parallel::clusterEvalQ(cl, {
  setwd(script_dir)
  library(rdrobust)
  library(rd2d)
}))

start_time <- Sys.time()

sim_status <- foreach(
  j = seq_len(m),
  .combine = rbind,
  .packages = c("rdrobust")
) %dopar% {
  status <- data.frame(replication = integer(), design = character(), dgp = integer(), tag = character(), suffix = character(), stringsAsFactors = FALSE)
  raw_rows <- list()

  for (g in seq_along(fuzzy_calibrations)) {
    cal <- fuzzy_calibrations[[g]]
    tag <- cal$tag
    set.seed(simulation_seed(seed_base, j, cal$design, cal$dgp))
    sim <- simulate_sample(cal, n)
    distance <- make_distances(sim$x, boundary, sim$assignment)

    for (b in seq_len(nrow(bandwidth_specs))) {
      bwparam <- bandwidth_specs$bwparam[b]
      suffix <- bandwidth_specs$suffix[b]
      cov_targets <- c("main", "itt", "fs")

      result.kinkoff <- rd2d.distance(sim$y, distance = distance, b = boundary, fuzzy = sim$fuzzy, bwparam = bwparam, repp = repp, vce = "hc1", params.cov = cov_targets)
      result.kinkon <- rd2d.distance(sim$y, distance = distance, b = boundary, fuzzy = sim$fuzzy, kink.unknown = c(TRUE, FALSE), bwparam = bwparam, repp = repp, vce = "hc1", params.cov = cov_targets)
      result.adaptive <- rd2d.distance(sim$y, distance = distance, b = boundary, fuzzy = sim$fuzzy, kink.position = kink.position, bwparam = bwparam, repp = repp, vce = "hc1", params.cov = cov_targets)
      h.rdrobust <- compute_rdrobust_bws(sim$y, sim$fuzzy, distance, bwparam)
      result.rdrobust <- rd2d.distance(sim$y, distance = distance, b = boundary, fuzzy = sim$fuzzy, h = h.rdrobust, bwparam = bwparam, repp = repp, vce = "hc1", params.cov = cov_targets)

      raw_rows[[length(raw_rows) + 1L]] <- raw_result_rows(result.kinkoff, cal$design, "kinkoff", j, cal$dgp, tag, bwparam, suffix, fuzzy_estimand_specs)
      raw_rows[[length(raw_rows) + 1L]] <- raw_result_rows(result.kinkon, cal$design, "kinkon", j, cal$dgp, tag, bwparam, suffix, fuzzy_estimand_specs)
      raw_rows[[length(raw_rows) + 1L]] <- raw_result_rows(result.adaptive, cal$design, "adaptive", j, cal$dgp, tag, bwparam, suffix, fuzzy_estimand_specs)
      raw_rows[[length(raw_rows) + 1L]] <- raw_result_rows(result.rdrobust, cal$design, "rdrobustadj", j, cal$dgp, tag, bwparam, suffix, fuzzy_estimand_specs)

      status <- rbind(status, data.frame(replication = j, design = cal$design, dgp = cal$dgp, tag = tag, suffix = suffix, stringsAsFactors = FALSE))
    }
  }

  for (g in seq_along(sharp_calibrations)) {
    cal <- sharp_calibrations[[g]]
    tag <- cal$tag
    set.seed(simulation_seed(seed_base, j, cal$design, cal$dgp))
    sim <- simulate_sample(cal, n)
    distance <- make_distances(sim$x, boundary, sim$assignment)
    bwparam <- "main"
    suffix <- "bwmain"
    cov_targets <- c("main", "main.0")

    result.kinkoff <- rd2d.distance(sim$y, distance = distance, b = boundary, params.other = "main.0", params.cov = cov_targets, repp = repp, vce = "hc1")
    result.kinkon <- rd2d.distance(sim$y, distance = distance, b = boundary, kink.unknown = c(TRUE, FALSE), params.other = "main.0", params.cov = cov_targets, repp = repp, vce = "hc1")
    result.adaptive <- rd2d.distance(sim$y, distance = distance, b = boundary, kink.position = kink.position, params.other = "main.0", params.cov = cov_targets, repp = repp, vce = "hc1")
    h.rdrobust <- compute_rdrobust_bws(sim$y, distance = distance)
    result.rdrobust <- rd2d.distance(sim$y, distance = distance, b = boundary, h = h.rdrobust, params.other = "main.0", params.cov = cov_targets, repp = repp, vce = "hc1")

    raw_rows[[length(raw_rows) + 1L]] <- raw_result_rows(result.kinkoff, cal$design, "kinkoff", j, cal$dgp, tag, bwparam, suffix, sharp_estimand_specs)
    raw_rows[[length(raw_rows) + 1L]] <- raw_result_rows(result.kinkon, cal$design, "kinkon", j, cal$dgp, tag, bwparam, suffix, sharp_estimand_specs)
    raw_rows[[length(raw_rows) + 1L]] <- raw_result_rows(result.adaptive, cal$design, "adaptive", j, cal$dgp, tag, bwparam, suffix, sharp_estimand_specs)
    raw_rows[[length(raw_rows) + 1L]] <- raw_result_rows(result.rdrobust, cal$design, "rdrobustadj", j, cal$dgp, tag, bwparam, suffix, sharp_estimand_specs)

    status <- rbind(status, data.frame(replication = j, design = cal$design, dgp = cal$dgp, tag = tag, suffix = suffix, stringsAsFactors = FALSE))
  }

  write_output(do.call(rbind, raw_rows), raw_replication_file(j))
  status
}

parallel::stopCluster(cl)
on.exit(NULL)

elapsed <- difftime(Sys.time(), start_time, units = "mins")
message(sprintf("Simulation outputs completed in %.2f minutes", as.numeric(elapsed)))
write_output(sim_status, file.path(output_dir, "simuls_status.csv"))

replication_files <- raw_replication_file(seq_len(m))
missing_replication_files <- replication_files[!file.exists(replication_files)]
if (length(missing_replication_files)) {
  stop(sprintf("Missing expected scratch simulation file: %s", missing_replication_files[1]), call. = FALSE)
}

raw_specs <- make_raw_result_specs()
message(sprintf("Consolidating %d scratch replication file(s) into %d raw CSV file(s)", length(replication_files), nrow(raw_specs)))
consolidate_raw_results(replication_files, raw_specs)

expected_raw <- basename(raw_result_file(raw_specs$dgp, raw_specs$design, raw_specs$method, raw_specs$tag, raw_specs$suffix))
missing_raw <- expected_raw[!file.exists(file.path(output_dir, expected_raw))]
if (length(missing_raw)) stop(sprintf("Missing expected raw simulation file: %s", missing_raw[1]), call. = FALSE)

cat("Simulation output complete.\n")
cat(sprintf("Raw simulation files written to: %s\n", output_dir))
cat(sprintf("Generated %d raw CSV file(s).\n", length(expected_raw)))
