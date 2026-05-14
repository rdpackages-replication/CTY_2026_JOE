################################################################################################
# Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods
# Empirical Application: SPP Data
# Replication output
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

suppressPackageStartupMessages(library(rdrobust))

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

get_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  as.integer(value)
}
emp_repp <- get_env_int("RD2D_EMP_REPP", 5000)
library(rd2d)

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

  dat <- raw[, c(
    "running_saber11", "running_sisben", "spadies_any",
    "eligible_spp", "beneficiary_spp"
  )]
  names(dat) <- c("x.1", "x.2", "y", "d", "w")
  dat <- dat[complete.cases(dat), , drop = FALSE]

  expected_assignment <- as.integer(dat$x.1 >= 0 & dat$x.2 >= 0)
  if (!all(dat$d == expected_assignment)) {
    stop("eligible_spp does not match the quadrant assignment rule.", call. = FALSE)
  }

  dat$y <- as.numeric(dat$y)
  dat$d <- as.numeric(dat$d)
  dat$w <- as.numeric(dat$w)
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
kink.position <- 11

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
fuzzy <- data$w
X <- cbind(data$x.1, data$x.2)
assignment <- data$d
b <- eval

distance <- proxy::dist(X, eval, method = "euclidean")
assignment_expanded <- matrix(rep(2 * assignment - 1, times = ncol(distance)), nrow = nrow(distance), ncol = ncol(distance))
distance <- distance * assignment_expanded

fuzzy_params_cov <- c("main", "itt", "fs")

# The sharp fits reproduce the old empirical application. They also provide
# the ITT bandwidths used below, so fuzzy$itt is directly comparable to the
# previous sharp main estimates.

# smooth boundary
result.itt.kinkoff <- rd2d.distance(Y, distance = distance, b = b, repp = emp_repp, vce = "hc1")
result.kinkoff <- rd2d.distance(
  Y, distance = distance, b = b, repp = emp_repp, vce = "hc1",
  fuzzy = fuzzy, bwparam = "itt", params.cov = fuzzy_params_cov
)

# unknown kink location
result.itt.kinkon <- rd2d.distance(
  Y, distance = distance, b = b, kink.unknown = c(TRUE, FALSE),
  repp = emp_repp, vce = "hc1"
)
result.kinkon <- rd2d.distance(
  Y, distance = distance, b = b, kink.unknown = c(TRUE, FALSE),
  repp = emp_repp, vce = "hc1",
  fuzzy = fuzzy, bwparam = "itt", params.cov = fuzzy_params_cov
)

# adaptive
result.itt.adaptive <- rd2d.distance(
  Y, distance = distance, b = b, kink.position = kink.position,
  repp = emp_repp, vce = "hc1"
)
result.adaptive <- rd2d.distance(
  Y, distance = distance, b = b, kink.position = kink.position,
  repp = emp_repp, vce = "hc1",
  fuzzy = fuzzy, bwparam = "itt", params.cov = fuzzy_params_cov
)

# rdrobust
bws <- matrix(0, nrow = neval, ncol = 2)
for (i in 1:neval) {
  out <- rdbwselect(Y, distance[, i], vce = "hc1")
  bws[i, 1] <- out$bws[1]
  bws[i, 2] <- out$bws[2]
}
result.itt.rdrobust <- rd2d.distance(Y, distance = distance, h = bws, b = b, repp = emp_repp, vce = "hc1")
result.rdrobust <- rd2d.distance(
  Y, distance = distance, h = bws, b = b, repp = emp_repp, vce = "hc1",
  fuzzy = fuzzy, bwparam = "itt", params.cov = fuzzy_params_cov
)

################################ Consistency Checks ############################

check_itt_match <- function(sharp, fuzzy, label, tol = 1e-10) {
  cols <- c(
    "estimate.p", "std.err.p", "estimate.q", "std.err.q",
    "t.value", "p.value", "ci.lower", "ci.upper",
    "h0", "h1", "h0.rbc", "h1.rbc", "N.Co", "N.Tr"
  )
  diffs <- vapply(cols, function(col) {
    max(abs(sharp$main[[col]] - fuzzy$itt[[col]]), na.rm = TRUE)
  }, numeric(1))
  if (any(diffs > tol)) {
    stop(
      sprintf(
        "fuzzy itt does not match old sharp main for %s. Max diff: %s = %.3e",
        label, names(which.max(diffs)), max(diffs)
      ),
      call. = FALSE
    )
  }
  if (!is.null(sharp$params.cov$main) && !is.null(fuzzy$params.cov$itt)) {
    cov_diff <- max(abs(sharp$params.cov$main - fuzzy$params.cov$itt), na.rm = TRUE)
    if (cov_diff > tol) {
      stop(
        sprintf(
          "fuzzy itt covariance does not match old sharp main for %s. Max diff: %.3e",
          label, cov_diff
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

sharp_results <- list(
  smooth = result.itt.kinkoff,
  adaptive = result.itt.adaptive,
  unknown_kink = result.itt.kinkon,
  rdrobust = result.itt.rdrobust
)

fuzzy_results <- list(
  smooth = result.kinkoff,
  adaptive = result.adaptive,
  unknown_kink = result.kinkon,
  rdrobust = result.rdrobust
)

for (nm in names(fuzzy_results)) {
  check_itt_match(sharp_results[[nm]], fuzzy_results[[nm]], nm)
}

################################## Output ######################################

as_result_table <- function(result, eval, output = c("main", "itt", "fs")) {
  output <- match.arg(output)
  point_res <- result[[output]]
  n <- nrow(point_res)
  aggregate_labels <- c("WBATE", "LBATE")

  invisible(capture.output(
    summ <- summary(
      result,
      output = output,
      cbands = output,
      WBATE = rep(1, n),
      LBATE = TRUE
    )
  ))

  tab <- summ$tables[[output]]
  if (!is.data.frame(tab)) {
    stop(sprintf("summary.rd2d.distance did not return a %s summary table.", output), call. = FALSE)
  }
  point_rows <- !(rownames(tab) %in% aggregate_labels)
  if (!identical(sum(point_rows), n)) {
    stop("summary.rd2d.distance returned an unexpected number of point rows.", call. = FALSE)
  }

  cbands <- summ$cbands[[output]]
  if (!is.data.frame(cbands) || !identical(nrow(cbands), n)) {
    stop(sprintf("summary.rd2d.distance did not return %s confidence bands.", output), call. = FALSE)
  }

  col_or_na <- function(name) {
    if (name %in% names(tab)) return(tab[[name]])
    rep(NA_real_, nrow(tab))
  }

  data.frame(
    row = rownames(tab),
    b1 = c(eval$x.1, rep(NA_real_, sum(!point_rows))),
    b2 = c(eval$x.2, rep(NA_real_, sum(!point_rows))),
    h0 = col_or_na("h0"),
    h1 = col_or_na("h1"),
    N.Co = col_or_na("N.Co"),
    N.Tr = col_or_na("N.Tr"),
    estimate.p = col_or_na("estimate.p"),
    ci.lower = col_or_na("ci.lower"),
    ci.upper = col_or_na("ci.upper"),
    cb.lower = col_or_na("cb.lower"),
    cb.upper = col_or_na("cb.upper"),
    estimate.q = col_or_na("estimate.q"),
    t.value = col_or_na("t.value"),
    p.value = col_or_na("p.value"),
    check.names = FALSE
  )
}

write_output_csv <- function(x, file) {
  write.csv(x, file.path(output_dir, file), row.names = FALSE, na = "")
}

old_files <- list.files(output_dir, pattern = "^empapp_.*[.]csv$", full.names = TRUE)
if (length(old_files) > 0) unlink(old_files)

outputs <- c(fuzzy = "main", itt = "itt", fs = "fs")
for (nm in names(fuzzy_results)) {
  for (output_name in names(outputs)) {
    write_output_csv(
      as_result_table(fuzzy_results[[nm]], eval, output = outputs[[output_name]]),
      paste0("empapp_", nm, "_", output_name, ".csv")
    )
  }
}

cat(sprintf(
  "Wrote %d empirical output file(s) to %s. Sharp-main/old ITT checks passed.\n",
  length(fuzzy_results) * length(outputs), output_dir
))
