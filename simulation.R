source('main_function.R')


library(doSNOW)
library(foreach)

# Setup execution parameters
number_of_replications <- as.integer(Sys.getenv("N_REP", "1000"))
replication_start <- as.integer(Sys.getenv("REP_START", "1"))
output_directory <- Sys.getenv("OUTPUT_DIR", "object_memory_mc_scale_n2000_20260814")
include_fp_sensitivity <- identical(tolower(Sys.getenv("FP_SENSITIVITY", "false")), "true")

rdata_dir <- file.path(output_directory, "RData_results")
dir.create(rdata_dir, showWarnings = FALSE, recursive = TRUE)

replications <- seq.int(replication_start, number_of_replications)

# Setup doSNOW Cluster
num_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

# Progress bar support
pb <- txtProgressBar(max = length(replications), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Export all locally defined functions to the parallel workers
export_funcs <- ls.str(mode = "function")

cat(sprintf("Starting %d replications on %d cores...\n", length(replications), num_cores))

# Parallel execution loop
results <- foreach(rep_id = replications, .options.snow = opts, .export = export_funcs) %dopar% {
  run_single_replication(
    replication = rep_id,
    output_dir = rdata_dir,
    include_fp_sensitivity = include_fp_sensitivity
  )
}

close(pb)
stopCluster(cl)
cat("Parallel simulation complete.\n")

# ----------------------------------------------------------------------------
# Result summary
# ----------------------------------------------------------------------------

run_simulation <- function(
    n_rep_total = 1000L,
    shard_index = 1L,
    shard_count = 1L,
    replication_start = 1L,
    n_values = c(250, 500, 1000, 1500, 2000),
    d_values = c(0, 0.1, 0.2, 0.3, 0.4),
    seed = 2026080811L,
    include_fp_sensitivity = FALSE,
    tunings = make_tuning_set(),
    psd_loading = 0.25,
    wasserstein_base_scale = 1,
    wasserstein_scale_loading = 0.25) {
  first_replication <- as.integer(replication_start) + as.integer(shard_index) - 1L
  replications <- seq.int(first_replication, n_rep_total, by = shard_count)
  n_values <- sort(unique(as.integer(n_values)))
  d_values <- sort(unique(as.numeric(d_values)))
  maximum_n <- max(n_values)
  maximum_lag <- max(vapply(tunings, function(tuning) {
    max(vapply(n_values, function(n) {
      floor(tuning$q * get_bandwidth(n, tuning))
    }, numeric(1)))
  }, numeric(1)))
  fp_specs <- lapply(tunings, function(tuning) {
    setNames(lapply(n_values, function(n) make_fp_spec(n, tuning)), n_values)
  })
  
  designs <- c("PSD/Frobenius", "Wasserstein", "Scalar", "PureLocationWasserstein")
  methods <- if (include_fp_sensitivity) {
    c(
      "Raw", "BC", "FP", "GlobalFP", "FPID", "LocalTotal",
      "PenalizedEqual", "PenalizedCP"
    )
  } else {
    c("Raw", "BC", "FP")
  }
  slopes <- c("Ratio", "LogSlope")
  aggregations <- c("Baseline", "TuneAvg")
  total_rows <- length(replications) * length(d_values) * length(n_values) *
    length(designs) * length(methods) * length(slopes) * length(aggregations)
  output <- data.frame(
    design = character(total_rows),
    n = integer(total_rows),
    d = numeric(total_rows),
    replication = integer(total_rows),
    slope = character(total_rows),
    method = character(total_rows),
    aggregation = character(total_rows),
    estimate = numeric(total_rows),
    projection_hit = logical(total_rows),
    id_identification = numeric(total_rows),
    id_penalty = numeric(total_rows),
    id_unit_penalty = numeric(total_rows),
    stringsAsFactors = FALSE
  )
  output$id_identification[] <- NA_real_
  output$id_penalty[] <- NA_real_
  output$id_unit_penalty[] <- NA_real_
  canonical_d <- c(0.1, 0.2, 0.3, 0.4, 0)
  cursor <- 0L
  for (replication in replications) {
    for (d in d_values) {
      d_seed_index <- match(d, canonical_d)
      if (is.na(d_seed_index)) {
        stop("Each d must belong to {0,.1,.2,.3,.4} for reproducible seeding.")
      }
      set.seed(seed + 100000L * d_seed_index + replication)
      phi <- runif(1L, -0.25, 0.25)
      theta <- runif(1L, -0.25, 0.25)
      gaussian_full <- simulate_arfima_driver(maximum_n, d, phi, theta)
      
      # Generating explicit objects and computing exact metric spaces 
      psd_objects <- generate_psd_matrices(
        gaussian_full, p = 3, loading = psd_loading
      )
      psd_full <- compute_actual_frobenius_distance(psd_objects)
      
      wasserstein_objects <- generate_distributional_objects(
        gaussian_full,
        base_scale = wasserstein_base_scale,
        scale_loading = wasserstein_scale_loading
      )
      wasserstein_full <- compute_actual_wasserstein_distance(wasserstein_objects)
      
      scalar_full <- compute_scalar_distance_matrix(gaussian_full)
      
      pl_wasserstein_objects <- generate_distributional_objects(
        gaussian_full, 
        base_scale = 1, 
        scale_loading = 0
      )
      pl_wasserstein_full <- compute_actual_wasserstein_distance(pl_wasserstein_objects)
      
      for (n in n_values) {
        inputs_by_design <- list(
          `PSD/Frobenius` = prepare_distance_inputs(
            psd_full, n, maximum_lag, tunings[[1L]]$block_counts
          ),
          Wasserstein = prepare_distance_inputs(
            wasserstein_full, n, maximum_lag, tunings[[1L]]$block_counts
          ),
          Scalar = prepare_distance_inputs(
            scalar_full, n, maximum_lag, tunings[[1L]]$block_counts
          ),
          PureLocationWasserstein = prepare_distance_inputs(
            pl_wasserstein_full, n, maximum_lag, tunings[[1L]]$block_counts
          )
        )
        for (design in names(inputs_by_design)) {
          estimates_by_tuning <- array(
            NA_real_,
            dim = c(length(methods), length(slopes), length(tunings)),
            dimnames = list(methods, slopes, names(tunings))
          )
          hits_by_tuning <- array(FALSE, dim = dim(estimates_by_tuning),
                                  dimnames = dimnames(estimates_by_tuning))
          diagnostics_by_tuning <- array(
            NA_real_,
            dim = c(3L, length(slopes), length(tunings)),
            dimnames = list(
              c("identification", "penalty", "unit_penalty"),
              slopes,
              names(tunings)
            )
          )
          for (tuning_name in names(tunings)) {
            family <- estimate_family(
              inputs_by_design[[design]],
              tunings[[tuning_name]],
              fp_specs[[tuning_name]][[as.character(n)]],
              include_fp_sensitivity = include_fp_sensitivity
            )
            estimates_by_tuning[, , tuning_name] <- family$estimate
            hits_by_tuning[, , tuning_name] <- family$projection_hit
            diagnostics_by_tuning[, , tuning_name] <- family$id_diagnostics
          }
          baseline <- estimates_by_tuning[, , "T0"]
          average <- apply(estimates_by_tuning, c(1L, 2L), mean)
          baseline_hits <- hits_by_tuning[, , "T0"]
          average_hits <- apply(hits_by_tuning, c(1L, 2L), any)
          baseline_diagnostics <- diagnostics_by_tuning[, , "T0"]
          average_diagnostics <- apply(
            diagnostics_by_tuning, c(1L, 2L), mean
          )
          for (aggregation in aggregations) {
            values <- if (aggregation == "Baseline") baseline else average
            hit_values <- if (aggregation == "Baseline") {
              baseline_hits
            } else {
              average_hits
            }
            diagnostic_values <- if (aggregation == "Baseline") {
              baseline_diagnostics
            } else {
              average_diagnostics
            }
            for (slope in slopes) for (method in methods) {
              cursor <- cursor + 1L
              is_fpid <- method == "FPID"
              output[cursor, ] <- list(
                design, n, d, replication, slope, method, aggregation,
                values[method, slope], hit_values[method, slope],
                if (is_fpid) diagnostic_values["identification", slope] else NA_real_,
                if (is_fpid) diagnostic_values["penalty", slope] else NA_real_,
                if (is_fpid) diagnostic_values["unit_penalty", slope] else NA_real_
              )
            }
          }
        }
      }
      rm(psd_objects, wasserstein_objects, pl_wasserstein_objects)
      rm(gaussian_full, psd_full, wasserstein_full, scalar_full, pl_wasserstein_full)
    }
    if (replication %% 25L == 0L) {
      cat("Completed replication", replication, "of", n_rep_total, "\n")
    }
  }
  output[seq_len(cursor), ]
}

summarize_simulation <- function(results) {
  key <- interaction(
    results$design, results$n, results$d, results$slope,
    results$method, results$aggregation, drop = TRUE
  )
  groups <- split(results, key)
  finite_mean <- function(x) {
    if (all(!is.finite(x))) return(NA_real_)
    mean(x[is.finite(x)])
  }
  summary <- do.call(rbind, lapply(groups, function(x) {
    error <- x$estimate - x$d
    data.frame(
      design = x$design[1L],
      n = x$n[1L],
      d = x$d[1L],
      slope = x$slope[1L],
      method = x$method[1L],
      aggregation = x$aggregation[1L],
      replications = nrow(x),
      mean = mean(x$estimate),
      bias = mean(error),
      sd = sd(x$estimate),
      mcse = sd(x$estimate) / sqrt(nrow(x)),
      rmse = sqrt(mean(error^2)),
      projection_rate = mean(x$projection_hit),
      mean_id_identification = finite_mean(x$id_identification),
      mean_id_penalty = finite_mean(x$id_penalty),
      unit_penalty_rate = finite_mean(x$id_unit_penalty),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary[order(
    summary$design, summary$slope, summary$n, summary$d,
    summary$method, summary$aggregation
  ), ]
}

aggregate_performance <- function(summary) {
  key <- interaction(
    summary$design, summary$slope, summary$method,
    summary$aggregation, drop = TRUE
  )
  groups <- split(summary, key)
  finite_mean <- function(x) {
    if (all(!is.finite(x))) return(NA_real_)
    mean(x[is.finite(x)])
  }
  answer <- do.call(rbind, lapply(groups, function(x) {
    data.frame(
      design = x$design[1L],
      slope = x$slope[1L],
      method = x$method[1L],
      aggregation = x$aggregation[1L],
      mean_absolute_bias = mean(abs(x$bias)),
      average_rmse = mean(x$rmse),
      maximum_projection_rate = max(x$projection_rate),
      mean_id_identification = finite_mean(x$mean_id_identification),
      mean_id_penalty = finite_mean(x$mean_id_penalty),
      unit_penalty_rate = finite_mean(x$unit_penalty_rate),
      stringsAsFactors = FALSE
    )
  }))
  rownames(answer) <- NULL
  answer[order(
    answer$design, answer$slope, answer$method, answer$aggregation
  ), ]
}

summarize_bias_crossings <- function(summary) {
  key <- interaction(
    summary$design, summary$d, summary$slope, summary$method,
    summary$aggregation, drop = TRUE
  )
  groups <- split(summary, key)
  answer <- do.call(rbind, lapply(groups, function(x) {
    x <- x[order(x$n), ]
    bias <- x$mean - x$d
    crossing_index <- which(
      head(bias, -1L) < 0 & tail(bias, -1L) >= 0
    )
    data.frame(
      design = x$design[1L],
      d = x$d[1L],
      slope = x$slope[1L],
      method = x$method[1L],
      aggregation = x$aggregation[1L],
      smallest_n = x$n[1L],
      bias_at_smallest_n = bias[1L],
      largest_n = x$n[nrow(x)],
      bias_at_largest_n = bias[length(bias)],
      under_to_over_crossing = length(crossing_index) > 0L,
      first_crossing_n = if (length(crossing_index)) {
        x$n[crossing_index[1L] + 1L]
      } else {
        NA_integer_
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(answer) <- NULL
  answer[order(
    answer$design, answer$slope, answer$method, answer$aggregation, answer$d
  ), ]
}

summarize_fp_screen <- function(summary) {
  method_order <- c(
    "BC", "GlobalFP", "LocalTotal", "PenalizedEqual", "PenalizedCP", "FP"
  )
  if (!all(method_order %in% summary$method)) return(NULL)
  baseline <- summary[summary$aggregation == "Baseline", ]
  region_statistics <- function(x) c(
    mean_absolute_bias = mean(abs(x$bias)),
    average_rmse = mean(x$rmse)
  )
  rows <- lapply(method_order, function(method) {
    one <- baseline[baseline$method == method, ]
    data.frame(
      method = method,
      all_mab = region_statistics(one)[1L],
      all_rmse = region_statistics(one)[2L],
      positive_d_mab = region_statistics(one[one$d > 0, ])[1L],
      positive_d_rmse = region_statistics(one[one$d > 0, ])[2L],
      d0_mab = region_statistics(one[one$d == 0, ])[1L],
      d0_rmse = region_statistics(one[one$d == 0, ])[2L],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ----------------------------------------------------------------------------
# LaTeX output
# ----------------------------------------------------------------------------

format_decimal <- function(x, digits = 3L) {
  answer <- sprintf(paste0("%.", digits, "f"), x)
  sub("^-0\\.", "-.", sub("^0\\.", ".", answer))
}

make_detailed_table <- function(summary, design, slope, label) {
  one <- summary[summary$design == design & summary$slope == slope, ]
  if(nrow(one) == 0) return(character(0))
  replication_count <- unique(one$replications)
  if (length(replication_count) != 1L) {
    stop("Detailed-table cells do not have a common replication count.")
  }
  replication_text <- format(replication_count, big.mark = ",", scientific = FALSE)
  replication_noun <- if (replication_count == 1L) {
    "replication"
  } else {
    "replications"
  }
  method_order <- c("Raw", "BC", "FP")
  headings <- c("Raw", "BC", "FP")
  
  design_text <- switch(design,
                        "Wasserstein" = "location--scale Wasserstein",
                        "PSD/Frobenius" = "PSD/Frobenius",
                        "Scalar" = "Scalar/Real-Valued",
                        "PureLocationWasserstein" = "Pure Location Wasserstein",
                        design)
  
  slope_text <- if (slope == "Ratio") "two-bandwidth ratio" else "multi-bandwidth log-slope"
  lines <- c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{Monte Carlo results for the ", design_text,
           " design using the ", slope_text,
           " estimator. Entries from ", replication_text,
           " ", replication_noun, " are means with RMSEs in parentheses.}"),
    paste0("\\label{", label, "}"),
    "\\fontsize{7.3pt}{8.4pt}\\selectfont", "\\setstretch{1.0}",
    "\\setlength{\\tabcolsep}{3.5pt}",
    "\\begin{tabular}{crrrrrr}", "\\toprule",
    "& \\multicolumn{3}{c}{Baseline} & \\multicolumn{3}{c}{Tuning average} \\\\",
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
    paste0("$\\mathrm d$ & ", paste(headings, collapse = " & "), " & ",
           paste(headings, collapse = " & "), " \\\\"), "\\midrule"
  )
  for (n in sort(unique(one$n))) {
    lines <- c(lines, paste0("\\multicolumn{7}{l}{$n=", n, "$} \\\\"))
    for (d in sort(unique(one$d))) {
      cells <- character(0)
      for (aggregation in c("Baseline", "TuneAvg")) {
        cells <- c(cells, vapply(method_order, function(method) {
          row <- one[one$n == n & one$d == d & one$method == method &
                       one$aggregation == aggregation, ]
          paste0(format_decimal(row$mean), " (", format_decimal(row$rmse), ")")
        }, character(1)))
      }
      lines <- c(lines, paste0(format_decimal(d, 2L), " & ",
                               paste(cells, collapse = " & "), " \\\\"))
    }
    if (n != max(one$n)) lines <- c(lines, "\\addlinespace")
  }
  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

write_detailed_tables <- function(summary, file) {
  lines <- c(
    make_detailed_table(summary, "PSD/Frobenius", "Ratio", "tab_mc_psd_ratio"), "",
    make_detailed_table(summary, "PSD/Frobenius", "LogSlope", "tab_mc_psd_ls"), "",
    make_detailed_table(summary, "Wasserstein", "Ratio", "tab_mc_wass_ratio"), "",
    make_detailed_table(summary, "Wasserstein", "LogSlope", "tab_mc_wass_ls"), "",
    make_detailed_table(summary, "Scalar", "Ratio", "tab_mc_scalar_ratio"), "",
    make_detailed_table(summary, "Scalar", "LogSlope", "tab_mc_scalar_ls"), "",
    make_detailed_table(summary, "PureLocationWasserstein", "Ratio", "tab_mc_plwass_ratio"), "",
    make_detailed_table(summary, "PureLocationWasserstein", "LogSlope", "tab_mc_plwass_ls")
  )
  writeLines(lines, file)
}

write_aggregate_table <- function(
    aggregate, file, cell_count, replication_count) {
  lines <- c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{Aggregate Monte Carlo performance over the ",
           cell_count,
           " $(n,\\mathrm d)$ cells in each design, based on ",
           format(replication_count, big.mark = ",", scientific = FALSE),
           if (replication_count == 1L) {
             " replication per cell; "
           } else {
             " replications per cell; "
           },
           "MAB denotes mean absolute bias.}"),
    "\\label{tab_mc_aggregate}", "\\small", "\\setstretch{1.0}",
    "\\begin{tabular}{lllrrrr}", "\\toprule",
    "& & & \\multicolumn{2}{c}{Baseline} & \\multicolumn{2}{c}{Tuning average} \\\\",
    "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
    "Design & Slope & Method & MAB & RMSE & MAB & RMSE \\\\ ", "\\midrule"
  )
  first_design <- TRUE
  for (design in c("PSD/Frobenius", "Wasserstein", "Scalar", "PureLocationWasserstein")) {
    for (slope in c("Ratio", "LogSlope")) {
      for (method in c("Raw", "BC", "FP")) {
        baseline <- aggregate[
          aggregate$design == design & aggregate$slope == slope &
            aggregate$method == method & aggregate$aggregation == "Baseline", ]
        average <- aggregate[
          aggregate$design == design & aggregate$slope == slope &
            aggregate$method == method & aggregate$aggregation == "TuneAvg", ]
        lines <- c(lines, paste(
          if (first_design) design else "",
          if (method == "Raw") {
            if (slope == "Ratio") "Ratio" else "Log slope"
          } else "",
          method,
          format_decimal(baseline$mean_absolute_bias, 4L),
          format_decimal(baseline$average_rmse, 4L),
          format_decimal(average$mean_absolute_bias, 4L),
          paste0(format_decimal(average$average_rmse, 4L), " \\\\"),
          sep = " & "
        ))
        first_design <- FALSE
      }
      if (slope == "Ratio") lines <- c(lines, "\\addlinespace")
    }
    if (design != "PureLocationWasserstein") {
      lines <- c(lines, "\\addlinespace\\midrule\\addlinespace")
      first_design <- TRUE
    }
  }
  writeLines(c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"), file)
}

write_fpid_table <- function(
    aggregate, file, cell_count, replication_count) {
  lines <- c(
    "\\begin{table}[H]", "\\centering",
    paste0(
      "\\caption{Aggregate comparison of FP and identification-adaptive ",
      "FP--ID over the ", cell_count,
      " $(n,\\mathrm d)$ cells in each design, based on ",
      format(replication_count, big.mark = ",", scientific = FALSE), " ",
      if (replication_count == 1L) {
        "replication per cell. MAB denotes mean absolute bias; "
      } else {
        "replications per cell. MAB denotes mean absolute bias; "
      },
      "$\\bar\\lambda_{\\mathrm{ID}}$ and $p_1$ are the average penalty and ",
      "the frequency of its unit-strength branch, respectively.}"
    ),
    "\\label{tab_mc_fpid}", "\\small", "\\setstretch{1.0}",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\begin{tabular}{lllrrrrr}", "\\toprule",
    paste0(
      "Design & Slope & Averaging & \\multicolumn{2}{c}{FP} & ",
      "\\multicolumn{2}{c}{FP--ID} & ",
      "$\\bar\\lambda_{\\mathrm{ID}}$ / $p_1$ \\\\"
    ),
    "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
    " & & & MAB & RMSE & MAB & RMSE & \\\\ ", "\\midrule"
  )
  first_design <- TRUE
  for (design in c("PSD/Frobenius", "Wasserstein", "Scalar", "PureLocationWasserstein")) {
    for (slope in c("Ratio", "LogSlope")) {
      for (aggregation in c("Baseline", "TuneAvg")) {
        fp <- aggregate[
          aggregate$design == design & aggregate$slope == slope &
            aggregate$method == "FP" &
            aggregate$aggregation == aggregation, ]
        fpid <- aggregate[
          aggregate$design == design & aggregate$slope == slope &
            aggregate$method == "FPID" &
            aggregate$aggregation == aggregation, ]
        lines <- c(lines, paste(
          if (first_design) design else "",
          if (aggregation == "Baseline") {
            if (slope == "Ratio") "Ratio" else "Log slope"
          } else "",
          if (aggregation == "Baseline") "No" else "Yes",
          format_decimal(fp$mean_absolute_bias, 4L),
          format_decimal(fp$average_rmse, 4L),
          format_decimal(fpid$mean_absolute_bias, 4L),
          format_decimal(fpid$average_rmse, 4L),
          paste0(
            format_decimal(fpid$mean_id_penalty, 3L), " / ",
            format_decimal(fpid$unit_penalty_rate, 3L), " \\\\"
          ),
          sep = " & "
        ))
        first_design <- FALSE
      }
      if (slope == "Ratio") lines <- c(lines, "\\addlinespace")
    }
    if (design != "PureLocationWasserstein") {
      lines <- c(lines, "\\addlinespace\\midrule\\addlinespace")
      first_design <- TRUE
    }
  }
  writeLines(c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"), file)
}

write_fp_screen_table <- function(
    fp_screen, file, cell_count, replication_count) {
  required <- c(
    "method", "all_mab", "all_rmse", "positive_d_mab",
    "positive_d_rmse", "d0_mab", "d0_rmse"
  )
  if (is.null(fp_screen) || !all(required %in% names(fp_screen))) {
    stop("FP-screen output is missing one or more required columns.")
  }
  method_labels <- c(
    BC = "BC",
    GlobalFP = "Global FP on $[0,.49]$",
    LocalTotal = "Local total",
    PenalizedEqual = "Unit penalty",
    PenalizedCP = "Rate penalty",
    FP = "Local last (reported FP)"
  )
  row_index <- match(names(method_labels), fp_screen$method)
  if (anyNA(row_index)) stop("FP-screen output is missing a required method.")
  one <- fp_screen[row_index, ]
  caption_start <- if (replication_count == 200L) {
    "Two-hundred-replication diagnostic"
  } else {
    paste0(
      format(replication_count, big.mark = ",", scientific = FALSE),
      "-replication diagnostic"
    )
  }
  lines <- c(
    "\\begin{table}[H]", "\\centering",
    paste0(
      "\\caption{", caption_start,
      " for FP specifications, averaged across the four design--slope ",
      "combinations. MAB denotes mean absolute bias.}"
    ),
    "\\label{tab_fp_screen}", "\\scriptsize", "\\setstretch{1.0}",
    "\\setlength{\\tabcolsep}{3.8pt}",
    "\\begin{tabular}{lrrrrrr}", "\\toprule",
    paste0(
      "& \\multicolumn{2}{c}{All ", cell_count, " cells} & ",
      "\\multicolumn{2}{c}{$\\mathrm d>0$} & ",
      "\\multicolumn{2}{c}{$\\mathrm d=0$} \\\\"
    ),
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
    "Method & MAB & RMSE & MAB & RMSE & MAB & RMSE \\\\ ",
    "\\midrule"
  )
  for (i in seq_len(nrow(one))) {
    lines <- c(lines, paste(
      unname(method_labels[one$method[i]]),
      format_decimal(one$all_mab[i], 4L),
      format_decimal(one$all_rmse[i], 4L),
      format_decimal(one$positive_d_mab[i], 4L),
      format_decimal(one$positive_d_rmse[i], 4L),
      format_decimal(one$d0_mab[i], 4L),
      paste0(format_decimal(one$d0_rmse[i], 4L), " \\\\"),
      sep = " & "
    ))
  }
  writeLines(c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"), file)
}

write_manuscript_table_bundle <- function(
    output_directory, include_fp_screen = FALSE) {
  aggregate_file <- file.path(output_directory, "aggregate_table.tex")
  detailed_file <- file.path(output_directory, "detailed_tables.tex")
  required_files <- c(aggregate_file, detailed_file)
  if (!all(file.exists(required_files))) {
    stop("The main LaTeX table files have not been generated.")
  }
  lines <- c(
    "% Auto-generated manuscript table bundle.",
    "% Requires booktabs, float, and setspace in the manuscript preamble.",
    "",
    readLines(aggregate_file, warn = FALSE),
    "",
    readLines(detailed_file, warn = FALSE)
  )
  if (include_fp_screen) {
    fp_screen_file <- file.path(output_directory, "fp_screen_table.tex")
    if (!file.exists(fp_screen_file)) {
      stop("The FP-screen LaTeX table has not been generated.")
    }
    lines <- c(lines, "", readLines(fp_screen_file, warn = FALSE))
  }
  bundle_file <- file.path(output_directory, "manuscript_tables.tex")
  writeLines(lines, bundle_file)
  bundle_file
}

validate_simulation_results <- function(results) {
  required <- c(
    "design", "n", "d", "replication", "slope", "method",
    "aggregation", "estimate", "projection_hit"
  )
  if (!all(required %in% names(results))) {
    stop("Simulation output is missing one or more required columns.")
  }
  if (!nrow(results) || !all(is.finite(results$estimate))) {
    stop("Simulation output contains a missing or non-finite estimate.")
  }
  invisible(TRUE)
}

write_outputs <- function(results, output_directory) {
  validate_simulation_results(results)
  dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
  summary <- summarize_simulation(results)
  aggregate <- aggregate_performance(summary)
  crossings <- summarize_bias_crossings(summary)
  fp_screen <- summarize_fp_screen(summary)
  cell_count <- length(unique(interaction(summary$n, summary$d, drop = TRUE)))
  replication_count <- unique(summary$replications)
  if (length(replication_count) != 1L) {
    stop("Simulation cells do not have a common replication count.")
  }
  write.csv(summary, file.path(output_directory, "simulation_summary.csv"),
            row.names = FALSE)
  write.csv(aggregate, file.path(output_directory, "aggregate_summary.csv"),
            row.names = FALSE)
  write.csv(crossings, file.path(output_directory, "bias_crossing_summary.csv"),
            row.names = FALSE)
  write.csv(summary[summary$d == 0, ],
            file.path(output_directory, "d0_summary.csv"), row.names = FALSE)
  write.csv(summary[summary$n == max(summary$n), ],
            file.path(output_directory, "largest_n_summary.csv"), row.names = FALSE)
  if (!is.null(fp_screen)) {
    write.csv(fp_screen, file.path(output_directory, "fp_screen_summary.csv"),
              row.names = FALSE)
    write_fp_screen_table(
      fp_screen, file.path(output_directory, "fp_screen_table.tex"),
      cell_count, replication_count
    )
  }
  write.csv(compute_population_exponent_check(),
            file.path(output_directory, "population_exponent_check.csv"),
            row.names = FALSE)
  write_detailed_tables(summary, file.path(output_directory, "detailed_tables.tex"))
  write_aggregate_table(aggregate,
                        file.path(output_directory, "aggregate_table.tex"), cell_count,
                        replication_count)
  list(
    summary = summary,
    aggregate = aggregate,
    crossings = crossings,
    fp_screen = fp_screen
  )
}


# ----------------------------------------------------------------------------
# Aggregation and Output Execution
# ----------------------------------------------------------------------------

output_directory <- Sys.getenv("OUTPUT_DIR", "object_memory_mc_scale_n2000_20260814")
rdata_dir <- file.path(output_directory, "RData_results")

# Ensure the simulation output exists
rdata_files <- list.files(rdata_dir, pattern = "^rep_.*\\.RData$", full.names = TRUE)
if (length(rdata_files) == 0) {
  stop(sprintf("No .RData files found in %s. Run '01_run_simulations.R' first.", rdata_dir))
}

cat(sprintf("Found %d replication files. Aggregating...\n", length(rdata_files)))

# Load into an isolated environment to prevent global overwrites
load_rep <- function(file_path) {
  env <- new.env()
  load(file_path, envir = env)
  return(env$output_rep)
}

# Bind into a single master dataframe
simulation_results <- do.call(rbind, lapply(rdata_files, load_rep))
validate_simulation_results(simulation_results)

# Execute identical output pipeline from the original main_0.R code
cat("Generating summaries and CSV tables...\n")
output <- write_outputs(simulation_results, output_directory)

cat("Generating manuscript-ready LaTeX tables...\n")
manuscript_table_file <- write_manuscript_table_bundle(
  output_directory,
  include_fp_screen = !is.null(output$fp_screen)
)

cat("Analysis complete. Check the output directory:\n")
cat(normalizePath(output_directory, winslash = "/", mustWork = TRUE), "\n")
