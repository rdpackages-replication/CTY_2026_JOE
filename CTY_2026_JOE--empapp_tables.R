################################################################################################
# Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods
# Empirical Application: SPP Data
# Tables
################################################################################################

rm(list = ls(all = TRUE))

# This script converts the CSV files produced by CTY_2026_JOE--empapp.R into
# manuscript-ready LaTeX tabular fragments. Each bandwidth rule is reported for
# the fuzzy Wald estimand, the ITT/reduced-form outcome, and the first stage.

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/"))
  }
  normalizePath(getwd(), winslash = "/")
}

setwd(get_script_dir())

output_dir <- "output"
tables_dir <- "tables"
aggregate_labels <- c("WBATE", "LBATE")
methods <- c("smooth", "adaptive", "unknown_kink", "rdrobust")
estimands <- c("fuzzy", "itt", "fs")
estimand_headers <- c(
  fuzzy = "$\\tau(\\bb)$",
  itt = "$\\tau_Y(\\bb)$",
  fs = "$\\tau_W(\\bb)$"
)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

read_empapp_output <- function(method, estimand) {
  path <- file.path(output_dir, paste0("empapp_", method, "_", estimand, ".csv"))
  if (!file.exists(path)) {
    stop(sprintf("Missing output file: %s. Run CTY_2026_JOE--empapp.R first.", path), call. = FALSE)
  }
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "row", "h0", "h1", "N.Co", "N.Tr", "estimate.p",
    "ci.lower", "ci.upper", "cb.lower", "cb.upper", "estimate.q", "t.value", "p.value"
  )
  missing <- setdiff(required, names(tab))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing columns: %s", path, paste(missing, collapse = ", ")), call. = FALSE)
  }
  tab$row <- as.character(tab$row)
  tab
}

out_list <- setNames(vector("list", length(methods)), methods)
for (method in methods) {
  out_list[[method]] <- setNames(vector("list", length(estimands)), estimands)
  for (estimand in estimands) {
    out_list[[method]][[estimand]] <- read_empapp_output(method, estimand)
  }
}

old_files <- list.files(tables_dir, pattern = "^(tab-)?empapp_.*[.]tex$", full.names = TRUE)
if (length(old_files) > 0) unlink(old_files)

fmt_num <- function(x, digits = 3) {
  if (is.na(x) || !is.finite(x)) return("")
  sprintf(paste0("%.", digits, "f"), x)
}

fmt_pvalue <- function(x) {
  if (is.na(x) || !is.finite(x)) return("")
  if (x < 0.001) return("0.000")
  sprintf("%.3f", x)
}

fmt_ci <- function(lo, hi) {
  if (is.na(lo) || is.na(hi) || !is.finite(lo) || !is.finite(hi)) return("")
  sprintf("$(%s,\\, %s)$", fmt_num(lo), fmt_num(hi))
}

fmt_row_label <- function(row_id) {
  if (row_id == "WBATE") return("$\\mathtt{WBATE}$")
  if (row_id == "LBATE") return("$\\mathtt{LBATE}$")
  sprintf("$\\mathbf{b}_{%s}$", row_id)
}

interval_limits <- function(tab) {
  use_cb <- is.finite(tab$cb.lower) & is.finite(tab$cb.upper)
  data.frame(
    lower = ifelse(use_cb, tab$cb.lower, tab$ci.lower),
    upper = ifelse(use_cb, tab$cb.upper, tab$ci.upper)
  )
}

display_rows <- function(tab) {
  point_rows <- tab[!(tab$row %in% aggregate_labels), , drop = FALSE]
  aggregate_rows <- tab[tab$row %in% aggregate_labels, , drop = FALSE]
  aggregate_rows <- aggregate_rows[match(aggregate_labels, aggregate_rows$row), , drop = FALSE]
  if (any(is.na(aggregate_rows$row))) {
    stop("Empirical output is missing WBATE or LBATE rows.", call. = FALSE)
  }

  out <- rbind(point_rows, aggregate_rows)
  ints <- interval_limits(out)

  data.frame(
    label = vapply(out$row, fmt_row_label, character(1)),
    h = out$h0,
    N.Co = out$N.Co,
    N.Tr = out$N.Tr,
    estimate = out$estimate.p,
    pvalue = out$p.value,
    lower = ints$lower,
    upper = ints$upper,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

latex_body <- function(rows) {
  lines <- character(0)
  for (i in seq_len(nrow(rows))) {
    is_aggregate <- grepl("WBATE|LBATE", rows$label[i])
    previous_aggregate <- i > 1 && grepl("WBATE|LBATE", rows$label[i - 1])
    if (is_aggregate && !previous_aggregate) {
      lines <- c(lines, "\\midrule")
    }
    lines <- c(lines, paste0(
      rows$label[i], " & ",
      fmt_num(rows$h[i]), " & ",
      fmt_num(rows$N.Co[i], 0), " & ",
      fmt_num(rows$N.Tr[i], 0), " & ",
      fmt_num(rows$estimate[i]), " & ",
      fmt_pvalue(rows$pvalue[i]), " & ",
      fmt_ci(rows$lower[i], rows$upper[i]), " \\\\"
    ))
  }
  lines
}

write_tabular <- function(rows, file, estimand) {
  lines <- c(
    "\\begin{tabular}{@{}crrrrrc@{}}",
    "\\toprule\\toprule",
    paste0(
      paste(
        "\\multicolumn{1}{c}{$\\bb\\in\\B$}",
        "\\multicolumn{1}{c}{$h$}",
        "\\multicolumn{1}{c}{$N_{\\mathrm{Co}}$}",
        "\\multicolumn{1}{c}{$N_{\\mathrm{Tr}}$}",
        paste0("\\multicolumn{1}{c}{", estimand_headers[[estimand]], "}"),
        "\\multicolumn{1}{c}{p-value}",
        "\\multicolumn{1}{c}{95\\% RBC CI}",
        sep = " & "
      ),
      " \\\\"
    ),
    "\\midrule",
    latex_body(rows),
    "\\bottomrule\\bottomrule",
    "\\end{tabular}"
  )
  writeLines(lines, file)
  invisible(lines)
}

check_generated_tables <- function(expected) {
  missing <- expected[!file.exists(expected)]
  if (length(missing) > 0) {
    stop(sprintf("Missing expected empirical table: %s", missing[1]), call. = FALSE)
  }
  for (file in expected) {
    txt <- paste(readLines(file, warn = FALSE), collapse = "\n")
    if (!grepl("WBATE", txt, fixed = TRUE) || !grepl("LBATE", txt, fixed = TRUE)) {
      stop(sprintf("Generated table is missing WBATE or LBATE: %s", file), call. = FALSE)
    }
    if (!grepl("$N_{\\mathrm{Co}}$", txt, fixed = TRUE) ||
        !grepl("$N_{\\mathrm{Tr}}$", txt, fixed = TRUE)) {
      stop(sprintf("Generated table is missing N.Co or N.Tr headers: %s", file), call. = FALSE)
    }
  }
  invisible(expected)
}

expected <- character(0)
for (method in methods) {
  for (estimand in estimands) {
    file <- file.path(tables_dir, sprintf("empapp_%s_%s.tex", method, estimand))
    write_tabular(display_rows(out_list[[method]][[estimand]]), file, estimand)
    expected <- c(expected, file)
  }
}

check_generated_tables(expected)

cat(sprintf("Wrote %d empirical table file(s) to %s.\n", length(expected), tables_dir))
