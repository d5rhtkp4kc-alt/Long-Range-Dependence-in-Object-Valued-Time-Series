source('main_function.R')

##########################################################################
## overlapping subsample variability and LaTeX tables
##########################################################################
subsample_variability <- function(D_full,
                                  tuning_name,
                                  block_fraction = 0.7,
                                  preserve_full_bandwidth = TRUE,
                                  step = 1,
                                  verbose = FALSE) {
  n <- nrow(D_full)
  tuning_full <- make_tuning(tuning_name)
  m_full <- get_bandwidth(n, tuning_full)
  
  b_target <- floor(block_fraction * n)
  if (b_target < 3L || b_target >= n) {
    stop("The initial subsample block length must satisfy 3 <= floor(block_fraction * n) < n.")
  }
  
  b <- b_target
  if (preserve_full_bandwidth) {
    candidate_b <- seq.int(b_target, n - 1L)
    same_bandwidth <- vapply(
      candidate_b,
      function(bb) get_bandwidth(bb, tuning_full) == m_full,
      logical(1)
    )
    
    if (!any(same_bandwidth)) {
      stop(
        "Could not find a subsample size b >= floor(block_fraction * n) ",
        "with the same base Bartlett bandwidth as the full sample."
      )
    }
    b <- candidate_b[which(same_bandwidth)[1L]]
  }
  
  tuning_sub <- make_tuning(tuning_name)
  m_sub <- get_bandwidth(b, tuning_sub)
  
  if (preserve_full_bandwidth && m_sub != m_full) {
    stop("Internal error: the selected subsample bandwidth does not match the full-sample bandwidth.")
  }
  
  if (verbose) {
    message(
      "Subsample setup: n = ", n,
      ", target b = ", b_target,
      ", selected b = ", b,
      ", full-sample m = ", m_full,
      ", subsample m = ", m_sub,
      "."
    )
  }
  
  starts <- seq.int(1L, n - b + 1L, by = step)
  out <- vector("list", length(starts))
  
  for (k in seq_along(starts)) {
    s <- starts[k]
    idx <- s:(s + b - 1L)
    D_sub <- D_full[idx, idx, drop = FALSE]
    
    # Recompute every sample-size-dependent quantity using the selected
    # subsample length b. 
    max_lag_sub <- floor(tuning_sub$q * get_bandwidth(b, tuning_sub))
    max_lag_sub <- min(max_lag_sub, b - 1L)
    
    inputs_sub <- prepare_distance_inputs(
      D_full = D_sub,
      n = b,
      max_lag = max_lag_sub,
      block_counts = tuning_sub$block_counts
    )
    
    fp_spec_sub <- make_fp_spec(b, tuning_sub)
    
    fit_sub <- estimate_family(
      inputs = inputs_sub,
      tuning = tuning_sub,
      fp_spec = fp_spec_sub,
      include_fp_sensitivity = FALSE
    )
    
    est <- fit_sub$estimate
    
    out[[k]] <- data.frame(
      start = s,
      end = s + b - 1L,
      Raw_Ratio = est["Raw", "Ratio"],
      BC_Ratio = est["BC", "Ratio"],
      FP_Ratio = est["FP", "Ratio"],
      Raw_LogSlope = est["Raw", "LogSlope"],
      BC_LogSlope = est["BC", "LogSlope"],
      FP_LogSlope = est["FP", "LogSlope"],
      check.names = FALSE
    )
    
    if (verbose && (k %% 10L == 0L || k == length(starts))) {
      message("Completed ", k, " of ", length(starts), " subsamples.")
    }
  }
  
  estimates <- do.call(rbind, out)
  estimate_cols <- setdiff(names(estimates), c("start", "end"))
  
  summary_table <- do.call(
    rbind,
    lapply(estimate_cols, function(v) {
      x <- estimates[[v]]
      pieces <- strsplit(v, "_", fixed = TRUE)[[1L]]
      method <- pieces[1L]
      estimator <- paste(pieces[-1L], collapse = "_")
      
      data.frame(
        Method = method,
        Estimator = estimator,
        Mean = mean(x, na.rm = TRUE),
        SD = sd(x, na.rm = TRUE),
        Q025 = as.numeric(quantile(x, 0.025, na.rm = TRUE, names = FALSE)),
        Median = median(x, na.rm = TRUE),
        Q975 = as.numeric(quantile(x, 0.975, na.rm = TRUE, names = FALSE)),
        Min = min(x, na.rm = TRUE),
        Max = max(x, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(summary_table) <- NULL
  
  list(
    n = n,
    target_block_size = b_target,
    block_size = b,
    requested_block_fraction = block_fraction,
    achieved_block_fraction = b / n,
    full_bandwidth = m_full,
    subsample_bandwidth = m_sub,
    bandwidth_preserved = (m_sub == m_full),
    n_subsamples = length(starts),
    estimates = estimates,
    summary = summary_table
  )
}

# Build a publication-ready LaTeX table. The 2.5% and 97.5% columns are
# empirical subsample quantiles.
make_latex_subsample_table <- function(full_estimate,
                                       subsample_summary,
                                       caption,
                                       label,
                                       digits = 3) {
  methods <- c("Raw", "BC", "FP")
  estimators <- c("Ratio", "LogSlope")
  
  # Match the numerical style used in the Monte Carlo tables:
  # .198 instead of 0.198, and -.021 instead of -0.021.
  fmt <- function(x) {
    z <- formatC(x, format = "f", digits = digits)
    z <- sub("^0\\.", ".", z)
    z <- sub("^-0\\.", "-.", z)
    z
  }
  
  lines <- c(
    "\\begin{table}[!htb]",
    "\\centering",
    "\\renewcommand{\\arraystretch}{0.6}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\fontsize{7.3pt}{9.4pt}\\selectfont",
    "\\setstretch{1.0}",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\begin{tabular}{@{}llrrrrr@{}}",
    "\\toprule",
    paste0(
      "Method & Estimator & \\shortstack{Full\\\\sample} & \\shortstack{Subsample\\\\mean} & \\shortstack{Subsample\\\\SD} & ",
      "\\shortstack{2.5\\%\\\\quantile} & \\shortstack{97.5\\%\\\\quantile} \\\\"
    ),
    "\\midrule"
  )
  
  for (method in methods) {
    for (estimator in estimators) {
      ss <- subsample_summary[
        subsample_summary$Method == method &
          subsample_summary$Estimator == estimator,
        ,
        drop = FALSE
      ]
      
      if (nrow(ss) != 1L) {
        stop("Could not uniquely match subsample summary row for ",
             method, " / ", estimator, ".")
      }
      
      estimator_label <- if (estimator == "LogSlope") "LogSlope" else "Ratio"
      
      lines <- c(
        lines,
        paste0(
          method, " & ", estimator_label, " & ",
          fmt(full_estimate[method, estimator]), " & ",
          fmt(ss$Mean), " & ",
          fmt(ss$SD), " & ",
          fmt(ss$Q025), " & ",
          fmt(ss$Q975), " \\\\"
        )
      )
    }
    if (method != tail(methods, 1L)) {
      lines <- c(lines, "\\addlinespace")
    }
  }
  
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  
  paste(lines, collapse = "\n")
}

##########################################################################
##########################################################################
## return
##########################################################################
##########################################################################
library(readr)
library(dplyr)
library(ggplot2)
library(transport)

df_intraday <- read_csv("intraday_returns_raw.csv")

dist_ts <- df_intraday %>%
  group_by(asset_name, day) %>%
  summarize(
    empirical_dist = list(log_return),
    obs_count = n(),
    .groups = "drop"
  )

##########################################################################
## Plot: Intraday return distributions

library(tidyverse)
library(plot3D)


df_intraday <- read_csv("intraday_returns_raw.csv") %>%
  mutate(day = as.Date(day)) %>%
  arrange(day)


selected_days <- df_intraday %>%
  distinct(day) %>%
  arrange(day) %>%
  slice_head(n = 258) %>%
  pull(day)

df_plot <- df_intraday %>%
  filter(day %in% selected_days) %>%
  mutate(
    Day = match(day, selected_days),
    Return_Pct = log_return * 100
  )

return_lim <- quantile(
  df_plot$Return_Pct,
  probs = c(0.005, 0.995),
  na.rm = TRUE
)

n_grid <- 250

return_grid <- seq(
  from = return_lim[1],
  to   = return_lim[2],
  length.out = n_grid
)

density_list <- df_plot %>%
  group_by(Day) %>%
  summarise(
    Density = list({
      
      x <- Return_Pct[is.finite(Return_Pct)]
      
      dd <- density(
        x,
        from = return_lim[1],
        to   = return_lim[2],
        n    = n_grid,
        na.rm = TRUE
      )
      
      dd$y
    }),
    .groups = "drop"
  )


density_mat <- do.call(
  rbind,
  density_list$Density
)

rownames(density_mat) <- density_list$Day

z_max <- max(density_mat, na.rm = TRUE)


plot_return_3D <- function(
    matrix_input,
    return_grid,
    title.input = "",
    theta.input = 145,
    phi.input = 25,
    expand.input = 0.7,
    line.adj = -2
) {
  
  ##########################################################
  ## Axes
  ##########################################################
  
  # x-axis: log return
  x_grid <- return_grid
  
  # y-axis: trading day
  y_grid <- seq_len(nrow(matrix_input))
  
  # plot3D requires:
  # rows of z = x values
  # columns of z = y values
  z_final <- t(matrix_input)
  
  
  ##########################################################
  ## Colors
  ##########################################################
  
  my_cols <- colorRampPalette(
    c(
      "#E0F7FA",
      "#80DEEA",
      "#FFF59D",
      "#FFEB3B"
    )
  )(100)
  
  
  ##########################################################
  ## Plot margins and font
  ##########################################################
  
  par(
    mar = c(2.2, 2.0, 1.0, 0.5),
    family = "serif"
  )
  
  
  ##########################################################
  ## 3D surface
  ##########################################################
  
  persp3D(
    x = x_grid,
    y = y_grid,
    z = z_final,
    
    col = my_cols,
    colkey = FALSE,
    
    facets = TRUE,
    border = "#4682B4",
    lwd = 0.10,
    alpha = 0.75,
    
    theta = theta.input,
    phi = phi.input,
    expand = expand.input,
    
    bty = "b2",
    ticktype = "detailed",
    
    # Actual density scale
    zlim = c(0, max(z_final, na.rm = TRUE)),
    clim = c(0, max(z_final, na.rm = TRUE)),
    
    xlab = "Log Return (%)",
    ylab = "Day",
    zlab = "Density",
    
    cex.axis = 0.8,
    cex.lab = 1.0
  )
  
  
  if (title.input != "") {
    mtext(
      title.input,
      side = 3,
      line = line.adj,
      font = 1,
      cex = 1
    )
  }
}

plot_return_3D(
  matrix_input = density_mat,
  return_grid = return_grid,
  title.input = ""
)

##########################################################################
## memory parameter
dist_list <- dist_ts$empirical_dist
n <- length(dist_list)

D_full <- matrix(0, nrow = n, ncol = n)
rownames(D_full) <- as.character(dist_ts$day)
colnames(D_full) <- as.character(dist_ts$day)

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    w_dist <- wasserstein1d(dist_list[[i]], dist_list[[j]], p = 2)
    D_full[i, j] <- w_dist
    D_full[j, i] <- w_dist
  }
}


tuning <- make_tuning("Distribution_Tuning") 
max_lag <- floor(tuning$q * get_bandwidth(n, tuning))

inputs <- prepare_distance_inputs(
  D_full = D_full, 
  n = n, 
  max_lag = max_lag, 
  block_counts = tuning$block_counts
)

# 4. Create the fixed-point specification and estimate
fp_spec <- make_fp_spec(n, tuning)
results <- estimate_family(
  inputs = inputs, 
  tuning = tuning, 
  fp_spec = fp_spec, 
  include_fp_sensitivity = FALSE
)


print(results$estimate)

return_full_estimate <- results$estimate
return_D_full <- D_full

##########################################################################
## Overlapping subsample variability with bandwidth-matched block length
return_subsample <- subsample_variability(
  D_full = return_D_full,
  tuning_name = "Distribution_Tuning",
  block_fraction = 0.7,
  preserve_full_bandwidth = TRUE,
  step = 1,
  verbose = TRUE
)

cat("\nReturn-distribution application: n =", return_subsample$n,
    ", target b =", return_subsample$target_block_size,
    ", selected b =", return_subsample$block_size,
    ", full m =", return_subsample$full_bandwidth,
    ", subsample m =", return_subsample$subsample_bandwidth,
    ", number of overlapping subsamples =", return_subsample$n_subsamples,
    "\n")
print(return_subsample$summary)

# Generate and save the LaTeX table.
return_latex <- make_latex_subsample_table(
  full_estimate = return_full_estimate,
  subsample_summary = return_subsample$summary,
  caption = paste0(
    "Memory-parameter estimates and overlapping-subsample variability for ",
    "the intraday return-distribution data. Starting from ",
    "$\\lfloor 0.7n\\rfloor$, the subsample block length is chosen as the ",
    "smallest value preserving the full-sample integer Bartlett bandwidth; ",
    "here $b=", return_subsample$block_size, "$ and $m=",
    return_subsample$subsample_bandwidth, "$. The last two columns are empirical ",
    "subsample quantiles."
  ),
  label = "tab:return_subsample",
  digits = 3
)

cat("\n\nLaTeX table: return-distribution application\n")
cat(return_latex, "\n")
writeLines(return_latex, con = "return_subsample_table.tex")



##########################################################################
##########################################################################
## energy composition
##########################################################################
##########################################################################
library(patchwork)
library(ggplot2)

compute_d_sphere = function(y){
  dot_products <- y %*% t(y)
  
  dot_products[dot_products > 1] <- 1
  dot_products[dot_products < -1] <- -1
  
  d <- acos(dot_products)
  diag(d) <- 0
  
  return(d)
}





y <- get(load('energy_final_residual.RData'))

##########################################################################
## plot

df <- read.csv('US_energy.csv')
df[df == "--"] <- 0
period_list <- list()
for (i in 2:dim(df)[2]) {
  current_vec <- as.numeric(df[,i])
  compositional_vec <- rep(0,7)
  #coal
  compositional_vec[1] <- current_vec[1]/sum(current_vec)
  #petroleum
  compositional_vec[2] <- sum(current_vec[2:3])/sum(current_vec)
  #gas
  compositional_vec[3] <- sum(current_vec[4:5])/sum(current_vec)
  #nuclear
  compositional_vec[4] <- current_vec[6]/sum(current_vec)
  #Conventional hydroelectric
  compositional_vec[5] <- current_vec[7]/sum(current_vec)
  #Renewables (wind, geothermal, biomass (total))
  compositional_vec[6] <- sum(current_vec[8:10])/sum(current_vec)
  #solar
  compositional_vec[7] <- sum(current_vec[11:12])/sum(current_vec)
  
  period_list[[i-1]] <- sqrt(compositional_vec)
}


y_matrix <- matrix(0,nrow = length(period_list), ncol = length(period_list[[1]]))
for (i in 1:length(period_list)) {
  y_matrix[i,] <- period_list[[i]]
}

y_matrix <- y_matrix[1:120,] 

month_dates <- seq(as.Date("2001-01-01"), as.Date("2010-12-01"), by = "month")

prepare_plot_data <- function(mat, x_axis_vals) {
  df <- as.data.frame(100 * mat^2)
  colnames(df) <- c('Coal','Petroleum','Gas','Nuclear','Conventional hydroelectric','Renewables','Solar')
  df$x <- x_axis_vals
  # Pivot to long format for ggplot2
  return(tidyr::pivot_longer(df, 
                             cols = -x, 
                             names_to = "Energy_Source", 
                             values_to = "value"))
}

custom_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "11")

create_styled_plot <- function(data, title) {
  
  start_date <- as.Date("2001-01-01")
  end_limit  <- as.Date("2010-12-01")
  month_breaks  <- seq(start_date, end_limit, by = "1 years")
  
  ggplot(data, aes(x = x, y = value, color = Energy_Source, linetype = Energy_Source)) +
    geom_line(linewidth = 0.7) +
    scale_linetype_manual(values = custom_linetypes) +
    scale_x_date(
      breaks = month_breaks,
      date_labels = "%Y %b",
      limits = c(start_date, end_limit)
    ) +
    labs(y = 'Percentage(%)', title = title) +
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text=element_text(size = 12),
          legend.text=element_text(size = 14),
          plot.margin=margin(5,15,5,5),
          panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0,size=14))
}
# Plot (1): original data
p1 <- create_styled_plot(prepare_plot_data(y_matrix, month_dates), "(a) Original time series")

# Plot (2): remove trend
p2 <- create_styled_plot(prepare_plot_data(y, month_dates), "(b) De-trended and de-seasonalized time series")


combined_plot <- p1 + p2 +
  plot_layout(nrow = 2, guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    # Increase the width of the legend lines
    legend.key.width = unit(3, "line"), 
    # Optional: add a bit of space between legend items
    legend.spacing.x = unit(0.5, 'cm')
  )

combined_plot
##########################################################################
## memory parameter
y <- y[1:120,]
n <- nrow(y)
D_full <- compute_d_sphere(y)

tuning <- make_tuning("Spherical_Tuning") 
max_lag <- floor(tuning$q * get_bandwidth(n, tuning))

inputs <- prepare_distance_inputs(
  D_full = D_full, 
  n = n, 
  max_lag = max_lag, 
  block_counts = tuning$block_counts
)

fp_spec <- make_fp_spec(n, tuning)
results <- estimate_family(
  inputs = inputs, 
  tuning = tuning, 
  fp_spec = fp_spec, 
  include_fp_sensitivity = FALSE
)


print(results$estimate)
#.      Ratio  LogSlope
#Raw 0.1977376 0.1912929
#BC  0.2609682 0.2538511
#FP  0.2669306 0.2594035

energy_full_estimate <- results$estimate
energy_D_full <- D_full

##########################################################################
## Overlapping subsample variability with bandwidth-matched block length
energy_subsample <- subsample_variability(
  D_full = energy_D_full,
  tuning_name = "Spherical_Tuning",
  block_fraction = 0.7,
  preserve_full_bandwidth = TRUE,
  step = 1,
  verbose = TRUE
)

cat("\nEnergy application: n =", energy_subsample$n,
    ", target b =", energy_subsample$target_block_size,
    ", selected b =", energy_subsample$block_size,
    ", full m =", energy_subsample$full_bandwidth,
    ", subsample m =", energy_subsample$subsample_bandwidth,
    ", number of overlapping subsamples =", energy_subsample$n_subsamples,
    "\n")
print(energy_subsample$summary)

# Generate and save the LaTeX table.
energy_latex <- make_latex_subsample_table(
  full_estimate = energy_full_estimate,
  subsample_summary = energy_subsample$summary,
  caption = paste0(
    "Memory-parameter estimates and overlapping-subsample variability for ",
    "the U.S. electricity-generation composition data. Starting from ",
    "$\\lfloor 0.7n\\rfloor$, the subsample block length is chosen as the ",
    "smallest value preserving the full-sample integer Bartlett bandwidth; ",
    "here $b=", energy_subsample$block_size, "$ and $m=",
    energy_subsample$subsample_bandwidth, "$. The last two columns are empirical ",
    "subsample quantiles."
  ),
  label = "tab:energy_subsample",
  digits = 3
)

cat("\n\nLaTeX table: energy application\n")
cat(energy_latex, "\n")
writeLines(energy_latex, con = "energy_subsample_table.tex")



