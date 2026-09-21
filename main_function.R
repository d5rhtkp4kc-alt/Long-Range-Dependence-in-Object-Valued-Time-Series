# ----------------------------------------------------------------------------
# 1. General helpers and tuning specifications
# ----------------------------------------------------------------------------

project_interval <- function(x, interval) {
  if (length(interval) != 2L || interval[1L] > interval[2L]) {
    stop("interval must be an increasing two-element vector.")
  }
  pmin(pmax(x, interval[1L]), interval[2L])
}

make_tuning <- function(
    name,
    m_multiplier = 1,
    m_exponent = 1 / 3,
    q = 2,
    eta = 1 / 8,
    ridge_constant = 3 / 8,
    block_counts = c(1, 2, 4, 8, 16),
    pilot_projection = c(-0.10, 0.45),
    correction_active = c(0.05, 0.45),
    output_projection = c(-0.25, 0.75),
    global_fp_domain = c(0, 0.49)) {
  if (!is.character(name) || length(name) != 1L) stop("Invalid name.")
  if (m_multiplier <= 0 || m_exponent <= 0 || m_exponent >= 0.5) {
    stop("Invalid bandwidth specification.")
  }
  if (q <= 1 || eta < 0 || ridge_constant < 0) {
    stop("Invalid q, eta, or ridge constant.")
  }
  block_counts <- sort(unique(as.integer(block_counts)))
  if (length(block_counts) < 2L || block_counts[1L] != 1L) {
    stop("block_counts must contain 1 and at least one larger count.")
  }
  if (pilot_projection[1L] >= pilot_projection[2L] ||
      correction_active[1L] <= 0 || correction_active[2L] >= 0.5 ||
      correction_active[1L] >= correction_active[2L]) {
    stop("Invalid pilot or correction interval.")
  }
  if (length(output_projection) != 2L ||
      output_projection[1L] >= output_projection[2L] ||
      length(global_fp_domain) != 2L ||
      global_fp_domain[1L] >= global_fp_domain[2L] ||
      global_fp_domain[1L] < output_projection[1L] ||
      global_fp_domain[2L] > output_projection[2L]) {
    stop("Invalid output or legacy global-FP interval.")
  }
  list(
    name = name,
    m_multiplier = m_multiplier,
    m_exponent = m_exponent,
    q = q,
    eta = eta,
    ridge_constant = ridge_constant,
    block_counts = block_counts,
    pilot_projection = pilot_projection,
    correction_active = correction_active,
    output_projection = output_projection,
    global_fp_domain = global_fp_domain
  )
}

make_tuning_set <- function() {
  list(
    T0 = make_tuning("T0", m_multiplier = 1, q = 2, eta = 1 / 8),
    Tm_low = make_tuning("Tm_low", m_multiplier = 0.8, q = 2, eta = 1 / 8),
    Tm_high = make_tuning("Tm_high", m_multiplier = 1.25, q = 2, eta = 1 / 8),
    Tq_low = make_tuning("Tq_low", m_multiplier = 1, q = 1.5, eta = 1 / 8),
    Tq_high = make_tuning("Tq_high", m_multiplier = 1, q = 2.5, eta = 1 / 8),
    Teta_low = make_tuning("Teta_low", m_multiplier = 1, q = 2, eta = 1 / 12),
    Teta_high = make_tuning("Teta_high", m_multiplier = 1, q = 2, eta = 1 / 6)
  )
}

get_bandwidth <- function(n, tuning) {
  max(3L, as.integer(floor(
    tuning$m_multiplier * n^tuning$m_exponent + 0.5 + 1e-10
  )))
}

get_ridge_weight <- function(n, m, tuning) {
  tuning$ridge_constant * (m / n)^tuning$eta
}

# ----------------------------------------------------------------------------
# 2. Standardized Gaussian ARFIMA driver
# ----------------------------------------------------------------------------

simulate_gaussian_toeplitz <- function(gamma) {
  gamma <- as.numeric(gamma)
  n <- length(gamma)
  if (n < 2L) stop("gamma must have length at least two.")
  embedding_length <- 2L * (n - 1L)
  first_row <- c(gamma, rev(gamma[2L:(n - 1L)]))
  eigenvalues <- Re(fft(first_row))
  if (min(eigenvalues) < -1e-7) {
    stop("The circulant embedding is not nonnegative definite.")
  }
  eigenvalues <- pmax(eigenvalues, 0)
  noise <- rnorm(embedding_length)
  Re(fft(sqrt(eigenvalues) * fft(noise), inverse = TRUE))[seq_len(n)] /
    embedding_length
}

get_arfima_correlation <- function(n, d, phi, theta, grid_size = 32768L) {
  grid_size <- max(as.integer(grid_size), 2L * as.integer(n))
  lambda <- 2 * pi * (0:(grid_size - 1L)) / grid_size
  sin_half <- abs(2 * sin(lambda / 2))
  arma_part <- (1 + theta^2 + 2 * theta * cos(lambda)) /
    (1 + phi^2 - 2 * phi * cos(lambda))
  if (d == 0) {
    spectrum <- arma_part
  } else {
    spectrum <- numeric(grid_size)
    spectrum[-1L] <- arma_part[-1L] * sin_half[-1L]^(-2 * d)
    arma_zero <- ((1 + theta) / (1 - phi))^2
    spectrum[1L] <- arma_zero * (grid_size / pi)^(2 * d) / (1 - 2 * d)
  }
  covariance <- Re(fft(spectrum, inverse = TRUE)) / grid_size
  covariance <- covariance[seq_len(n)]
  covariance / covariance[1L]
}

simulate_arfima_driver <- function(n, d, phi, theta) {
  simulate_gaussian_toeplitz(get_arfima_correlation(n, d, phi, theta))
}

compute_population_exponent_check <- function(
    d_values = c(0.1, 0.2, 0.3, 0.4),
    lag_range = 100:500,
    psd_loading = 0.25,
    wasserstein_base_scale = 1,
    wasserstein_scale_loading = 0.25,
    quadrature_order = 32L) {
  standard_normal_quadrature <- function(quadrature_size) {
    quadrature_size <- as.integer(quadrature_size)
    if (quadrature_size < 2L) {
      stop("quadrature_order must be at least two.")
    }
    jacobi <- matrix(0, quadrature_size, quadrature_size)
    off_diagonal <- sqrt(seq_len(quadrature_size - 1L) / 2)
    jacobi[cbind(seq_len(quadrature_size - 1L), 2:quadrature_size)] <-
      off_diagonal
    jacobi[cbind(2:quadrature_size, seq_len(quadrature_size - 1L))] <-
      off_diagonal
    eig <- eigen(jacobi, symmetric = TRUE)
    index <- order(eig$values)
    list(
      nodes = sqrt(2) * eig$values[index],
      weights = eig$vectors[1L, index]^2
    )
  }
  standard_exponential_quadrature <- function(quadrature_size) {
    quadrature_size <- as.integer(quadrature_size)
    if (quadrature_size < 2L) {
      stop("quadrature_order must be at least two.")
    }
    jacobi <- diag(2 * seq_len(quadrature_size) - 1)
    off_diagonal <- seq_len(quadrature_size - 1L)
    jacobi[cbind(seq_len(quadrature_size - 1L), 2:quadrature_size)] <-
      off_diagonal
    jacobi[cbind(2:quadrature_size, seq_len(quadrature_size - 1L))] <-
      off_diagonal
    eig <- eigen(jacobi, symmetric = TRUE)
    index <- order(eig$values)
    list(
      nodes = eig$values[index],
      weights = eig$vectors[1L, index]^2
    )
  }
  normal_quadrature <- standard_normal_quadrature(quadrature_order)
  exponential_quadrature <- standard_exponential_quadrature(quadrature_order)
  normal_node <- rep(normal_quadrature$nodes, each = quadrature_order)
  exponential_node <- rep(
    exponential_quadrature$nodes, times = quadrature_order
  )
  joint_weight <- rep(normal_quadrature$weights, each = quadrature_order) *
    rep(exponential_quadrature$weights, times = quadrature_order)
  
  expected_wasserstein_distance <- function(correlation) {
    midpoint <- sqrt((1 + correlation) / 2) * normal_node
    half_difference <- sqrt(
      pmax(1 - correlation, 0) * exponential_node
    )
    scale_ratio <- sinh(
      wasserstein_scale_loading * half_difference
    ) / half_difference
    curvature <- sqrt(
      1 + wasserstein_base_scale^2 *
        exp(2 * wasserstein_scale_loading * midpoint) * scale_ratio^2
    )
    2 * sqrt(pmax(1 - correlation, 0)) / sqrt(pi) *
      sum(joint_weight * curvature)
  }
  
  expected_psd_distance <- function(correlation) {
    second_factor <- 2 * integrate(function(z) {
      sqrt(1 + psd_loading^2 * (1 + correlation) * z^2) * dnorm(z)
    }, lower = 0, upper = Inf, rel.tol = 1e-10)$value
    2 * psd_loading * sqrt(1 - correlation) * sqrt(2 / pi) * second_factor
  }
  
  expected_scalar_distance <- function(correlation) {
    2 * sqrt(pmax(1 - correlation, 0)) / sqrt(pi)
  }
  
  psd_zero <- expected_psd_distance(0)
  wasserstein_zero <- expected_wasserstein_distance(0)
  scalar_zero <- expected_scalar_distance(0)
  derivative_step <- 1e-5
  wasserstein_linear_coefficient <- (
    wasserstein_zero - expected_wasserstein_distance(derivative_step)
  ) / derivative_step
  rows <- lapply(d_values, function(d) {
    correlation <- get_arfima_correlation(max(lag_range) + 1L, d, 0, 0)
    rho <- correlation[lag_range + 1L]
    covariance_scale <- rho
    psd_scale <- psd_zero - vapply(rho, expected_psd_distance, numeric(1))
    wasserstein_scale <- wasserstein_zero - vapply(
      rho, expected_wasserstein_distance, numeric(1)
    )
    scalar_scale <- scalar_zero - vapply(
      rho, expected_scalar_distance, numeric(1)
    )
    effective <- function(scale) {
      slope <- coef(lm(log(scale) ~ log(lag_range)))[2L]
      unname((slope + 1) / 2)
    }
    data.frame(
      d = d,
      gaussian_covariance = effective(covariance_scale),
      psd_frobenius_C = effective(psd_scale),
      wasserstein_location_scale_C = effective(wasserstein_scale),
      scalar_real_valued_C = effective(scalar_scale),
      wasserstein_pure_location_C = effective(scalar_scale),
      wasserstein_linear_coefficient = wasserstein_linear_coefficient,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ----------------------------------------------------------------------------
# 3. Object Generation and Native Distance Computations
# ----------------------------------------------------------------------------

generate_psd_matrices <- function(G, p = 3, loading = 0.25) {
  n <- length(G)
  M0 <- diag(p)
  u <- rep(0, p)
  u[1] <- 1
  v <- rep(0, p)
  v[2] <- 1
  
  X_list <- vector("list", n)
  for (t in seq_len(n)) {
    vec <- u + loading * G[t] * v
    X_list[[t]] <- M0 + vec %*% t(vec)
  }
  return(X_list)
}

compute_actual_frobenius_distance <- function(X_list) {
  n <- length(X_list)
  D <- matrix(0, n, n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      dist_val <- sqrt(sum((X_list[[i]] - X_list[[j]])^2))
      D[i, j] <- dist_val
      D[j, i] <- dist_val
    }
  }
  return(D)
}

generate_distributional_objects <- function(G, base_scale = 1, scale_loading = 0.25, grid_size = 1000) {
  n <- length(G)
  u_grid <- seq(0.5 / grid_size, 1 - 0.5 / grid_size, length.out = grid_size)
  std_quantiles <- qnorm(u_grid)
  # Normalize to exactly match the theoretical Wasserstein distance in population
  std_quantiles <- std_quantiles / sqrt(mean(std_quantiles^2))
  
  dist_list <- vector("list", n)
  for (t in seq_len(n)) {
    mu_t <- G[t]
    sigma_t <- base_scale * exp(scale_loading * mu_t)
    dist_list[[t]] <- mu_t + sigma_t * std_quantiles
  }
  return(dist_list)
}

compute_actual_wasserstein_distance <- function(dist_list) {
  n <- length(dist_list)
  D <- matrix(0, n, n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      dist_val <- sqrt(mean((dist_list[[i]] - dist_list[[j]])^2))
      D[i, j] <- dist_val
      D[j, i] <- dist_val
    }
  }
  return(D)
}

compute_scalar_distance_matrix <- function(G) {
  answer <- abs(outer(G, G, "-"))
  diag(answer) <- 0
  answer
}

prepare_distance_inputs <- function(D_full, n, max_lag, block_counts) {
  index <- seq_len(n)
  D_hat <- sum(D_full[index, index, drop = FALSE]) / (n * (n - 1))
  D_blocks <- vapply(block_counts, function(s) {
    if (s == 1L) return(D_hat)
    ell <- floor(n / s)
    mean(vapply(seq_len(s), function(j) {
      block_index <- ((j - 1L) * ell + 1L):(j * ell)
      sum(D_full[block_index, block_index, drop = FALSE]) /
        (ell * (ell - 1))
    }, numeric(1)))
  }, numeric(1))
  names(D_blocks) <- as.character(block_counts)
  d_hat <- vapply(seq_len(max_lag), function(k) {
    mean(D_full[cbind((k + 1L):n, seq_len(n - k))])
  }, numeric(1))
  list(n = n, D_hat = D_hat, D_blocks = D_blocks, d_hat = d_hat)
}

# ----------------------------------------------------------------------------
# 4. Two-sided punctured Bartlett scales and slope estimators
# ----------------------------------------------------------------------------

get_signed_scale_curve <- function(
    inputs,
    D_estimate,
    tuning,
    ridge_estimate = inputs$D_hat) {
  n <- inputs$n
  m <- get_bandwidth(n, tuning)
  maximum_r <- floor(tuning$q * m)
  if (maximum_r >= n || length(inputs$d_hat) < maximum_r) {
    stop("Insufficient lagged-distance averages.")
  }
  C_hat <- D_estimate - inputs$d_hat
  ridge <- 2 * get_ridge_weight(n, m, tuning) * ridge_estimate
  signed <- rep(NA_real_, maximum_r)
  for (r in 2:maximum_r) {
    k <- seq_len(r - 1L)
    signed[r] <- ridge + 2 * sum((1 - k / r) * C_hat[k])
  }
  signed
}

evaluate_slope <- function(signed_curve, n, tuning, type) {
  type <- match.arg(type, c("Ratio", "LogSlope"))
  m <- get_bandwidth(n, tuning)
  maximum_r <- floor(tuning$q * m)
  if (type == "Ratio") {
    values <- abs(signed_curve[c(m, maximum_r)])
    if (any(!is.finite(values)) || any(values <= 0)) return(0)
    return(log(values[2L] / values[1L]) /
             (2 * log(maximum_r / m)))
  }
  grid <- m:maximum_r
  values <- abs(signed_curve[grid])
  if (any(!is.finite(values)) || any(values <= 0)) return(0)
  x <- log(grid)
  x_centered <- x - mean(x)
  sum(x_centered * log(values)) / (2 * sum(x_centered^2))
}

estimate_from_D <- function(inputs, D_estimate, tuning, type) {
  unprojected <- evaluate_slope(
    get_signed_scale_curve(
      inputs,
      D_estimate,
      tuning,
      ridge_estimate = inputs$D_hat
    ),
    inputs$n,
    tuning,
    type
  )
  projected <- project_interval(unprojected, tuning$output_projection)
  list(
    estimate = projected,
    unprojected = unprojected,
    projection_hit = abs(projected - unprojected) > 1e-12
  )
}

# ----------------------------------------------------------------------------
# 5. Block-difference correction
# ----------------------------------------------------------------------------

finite_power_bias_factor <- function(block_length, d) {
  block_length <- as.integer(block_length)
  if (block_length < 2L || !is.finite(d) || d >= 0.5) return(NA_real_)
  lag <- seq_len(block_length - 1L)
  2 * sum((block_length - lag) * lag^(2 * d - 1)) /
    (block_length * (block_length - 1))
}

estimate_common_bias <- function(inputs, pilot, tuning) {
  pilot_used <- project_interval(pilot, tuning$pilot_projection)
  active <- tuning$correction_active[1L] < pilot_used &&
    pilot_used <= tuning$correction_active[2L]
  if (!is.finite(pilot_used) || !active) {
    return(list(bias = NA_real_, correction = 0, active = FALSE))
  }
  counts <- tuning$block_counts
  lengths <- floor(inputs$n / counts)
  factors <- vapply(lengths, finite_power_bias_factor, numeric(1), d = pilot_used)
  regressors <- factors[-1L] / factors[1L] - 1
  response <- inputs$D_blocks[as.character(counts[-1L])] - inputs$D_hat
  denominator <- sum(regressors^2)
  if (!is.finite(denominator) || denominator <= 1e-14) {
    return(list(bias = NA_real_, correction = 0, active = FALSE))
  }
  bias <- min(sum(regressors * response) / denominator, 0)
  list(bias = bias, correction = -bias, active = TRUE)
}

estimate_update <- function(inputs, pilot, tuning, type) {
  adjustment <- estimate_common_bias(inputs, pilot, tuning)
  result <- estimate_from_D(
    inputs,
    inputs$D_hat + adjustment$correction,
    tuning,
    type
  )
  c(result, adjustment)
}

# ----------------------------------------------------------------------------
# 6. BC-centered fixed-point refinement and sensitivity variants
# ----------------------------------------------------------------------------

make_fp_spec <- function(n, tuning, maximum_grid_step = 0.0025) {
  m <- get_bandwidth(n, tuning)
  maximum_r <- floor(tuning$q * m)
  grid_step <- min(maximum_grid_step, n^(-1 / 2))
  search_domain <- tuning$output_projection
  d_grid <- seq(search_domain[1L], search_domain[2L], by = grid_step)
  if (tail(d_grid, 1L) < search_domain[2L]) {
    d_grid <- c(d_grid, search_domain[2L])
  }
  pilot_grid <- project_interval(d_grid, tuning$pilot_projection)
  lengths <- floor(n / tuning$block_counts)
  factor_matrix <- vapply(pilot_grid, function(d) {
    vapply(lengths, finite_power_bias_factor, numeric(1), d = d)
  }, numeric(length(lengths)))
  if (is.null(dim(factor_matrix))) {
    factor_matrix <- matrix(factor_matrix, nrow = length(lengths))
  }
  regressor_matrix <- sweep(
    factor_matrix[-1L, , drop = FALSE],
    2L,
    factor_matrix[1L, ],
    "/"
  ) - 1
  active <- tuning$correction_active[1L] < pilot_grid &
    pilot_grid <= tuning$correction_active[2L]
  list(
    n = n,
    m = m,
    maximum_r = maximum_r,
    d_grid = d_grid,
    regressor_matrix = regressor_matrix,
    denominator = colSums(regressor_matrix^2),
    active = active,
    tuning = tuning
  )
}

parabolic_grid_minimum <- function(objective, grid) {
  best <- which.min(objective)
  estimate <- grid[best]
  if (best > 1L && best < length(grid)) {
    index <- (best - 1L):(best + 1L)
    x <- grid[index]
    y <- objective[index]
    fit <- try(lm.fit(cbind(1, x, x^2), y), silent = TRUE)
    if (!inherits(fit, "try-error") && is.finite(fit$coefficients[3L]) &&
        fit$coefficients[3L] > 0) {
      vertex <- -fit$coefficients[2L] / (2 * fit$coefficients[3L])
      if (is.finite(vertex) && vertex >= x[1L] && vertex <= x[3L]) {
        estimate <- vertex
      }
    }
  }
  estimate
}

interval_grid_minimum <- function(objective, grid, interval, fallback) {
  interval <- c(max(interval[1L], min(grid)), min(interval[2L], max(grid)))
  if (interval[1L] > interval[2L]) return(fallback)
  index <- which(grid >= interval[1L] & grid <= interval[2L])
  if (!length(index)) return(unname(project_interval(fallback, interval)))
  estimate <- parabolic_grid_minimum(objective[index], grid[index])
  unname(project_interval(estimate, interval))
}

get_fp_map <- function(inputs, spec, type) {
  tuning <- spec$tuning
  response <- inputs$D_blocks[
    as.character(tuning$block_counts[-1L])
  ] - inputs$D_hat
  numerator <- colSums(spec$regressor_matrix * response)
  bias <- numerator / pmax(spec$denominator, 1e-14)
  correction <- pmax(-bias, 0)
  correction[!spec$active] <- 0
  base_curve <- get_signed_scale_curve(
    inputs,
    inputs$D_hat,
    tuning,
    ridge_estimate = inputs$D_hat
  )
  if (type == "Ratio") {
    r <- c(spec$m, spec$maximum_r)
    signed <- outer(correction, r - 1L) +
      matrix(base_curve[r], nrow = length(correction), ncol = 2L, byrow = TRUE)
    map <- log(abs(signed[, 2L]) / abs(signed[, 1L])) /
      (2 * log(spec$maximum_r / spec$m))
  } else {
    r <- spec$m:spec$maximum_r
    signed <- outer(correction, r - 1L) +
      matrix(base_curve[r], nrow = length(correction), ncol = length(r), byrow = TRUE)
    log_scale <- log(abs(signed))
    x <- log(r)
    x_centered <- x - mean(x)
    map <- as.numeric(log_scale %*% x_centered) /
      (2 * sum(x_centered^2))
  }
  map[!is.finite(map)] <- 0
  project_interval(map, tuning$output_projection)
}

estimate_grid_derivative <- function(values, grid, at) {
  if (length(values) != length(grid) || length(grid) < 2L) {
    stop("values and grid must have the same length, at least two.")
  }
  index <- which.min(abs(grid - at))
  if (index == 1L) {
    return((values[2L] - values[1L]) / (grid[2L] - grid[1L]))
  }
  if (index == length(grid)) {
    return((values[index] - values[index - 1L]) /
             (grid[index] - grid[index - 1L]))
  }
  (values[index + 1L] - values[index - 1L]) /
    (grid[index + 1L] - grid[index - 1L])
}

estimate_fp_candidates <- function(
    inputs, spec, type, raw_estimate, first_update, bc_estimate) {
  map <- get_fp_map(inputs, spec, type)
  fp_objective <- (map - spec$d_grid)^2
  fp <- interval_grid_minimum(
    fp_objective, spec$d_grid, spec$tuning$global_fp_domain, bc_estimate
  )
  map_derivative <- estimate_grid_derivative(map, spec$d_grid, fp)
  identification <- abs(1 - map_derivative)
  resolution <- (spec$m - 1) / (2 * spec$n)
  penalty <- resolution / max(identification^2, resolution)
  fpid_objective <- fp_objective + penalty *
    (spec$d_grid - first_update)^2
  fpid <- interval_grid_minimum(
    fpid_objective, spec$d_grid, spec$tuning$global_fp_domain, first_update
  )
  
  output_domain <- spec$tuning$output_projection
  last_radius <- abs(bc_estimate - first_update)
  total_radius <- abs(bc_estimate - raw_estimate)
  local_last <- interval_grid_minimum(
    fp_objective,
    spec$d_grid,
    bc_estimate + c(-last_radius, last_radius),
    bc_estimate
  )
  local_total <- interval_grid_minimum(
    fp_objective,
    spec$d_grid,
    bc_estimate + c(-total_radius, total_radius),
    bc_estimate
  )
  equal_objective <- fp_objective + (spec$d_grid - bc_estimate)^2
  penalized_equal <- interval_grid_minimum(
    equal_objective, spec$d_grid, output_domain, bc_estimate
  )
  bounded_pilot <- project_interval(first_update, c(0, 0.45))
  cp_penalty <- (spec$m / spec$n)^(1 - 2 * bounded_pilot)
  cp_objective <- fp_objective + cp_penalty *
    (spec$d_grid - bc_estimate)^2
  penalized_cp <- interval_grid_minimum(
    cp_objective, spec$d_grid, output_domain, bc_estimate
  )
  list(
    GlobalFP = fp,
    FPID = fpid,
    FP = local_last,
    LocalTotal = local_total,
    PenalizedEqual = penalized_equal,
    PenalizedCP = penalized_cp,
    identification = identification,
    penalty = penalty,
    unit_penalty = as.numeric(identification^2 <= resolution),
    last_radius = last_radius,
    total_radius = total_radius,
    cp_penalty = cp_penalty
  )
}

estimate_family <- function(
    inputs, tuning, fp_spec, include_fp_sensitivity = FALSE) {
  methods <- if (include_fp_sensitivity) {
    c(
      "Raw", "BC", "FP", "GlobalFP", "FPID", "LocalTotal",
      "PenalizedEqual", "PenalizedCP"
    )
  } else {
    c("Raw", "BC", "FP")
  }
  slopes <- c("Ratio", "LogSlope")
  estimates <- matrix(
    NA_real_,
    nrow = length(methods),
    ncol = length(slopes),
    dimnames = list(methods, slopes)
  )
  hits <- matrix(FALSE, nrow = length(methods), ncol = length(slopes),
                 dimnames = list(methods, slopes))
  id_diagnostics <- matrix(
    NA_real_,
    nrow = 3L,
    ncol = length(slopes),
    dimnames = list(c("identification", "penalty", "unit_penalty"), slopes)
  )
  for (type in slopes) {
    raw <- estimate_from_D(inputs, inputs$D_hat, tuning, type)
    iter1 <- estimate_update(inputs, raw$estimate, tuning, type)
    iter2 <- estimate_update(inputs, iter1$estimate, tuning, type)
    fp_pair <- estimate_fp_candidates(
      inputs,
      fp_spec,
      type,
      raw_estimate = raw$estimate,
      first_update = iter1$estimate,
      bc_estimate = iter2$estimate
    )
    all_estimates <- c(
      Raw = unname(raw$estimate),
      BC = unname(iter2$estimate),
      FP = unname(fp_pair$FP),
      GlobalFP = unname(fp_pair$GlobalFP),
      FPID = unname(fp_pair$FPID),
      LocalTotal = unname(fp_pair$LocalTotal),
      PenalizedEqual = unname(fp_pair$PenalizedEqual),
      PenalizedCP = unname(fp_pair$PenalizedCP)
    )
    estimates[, type] <- all_estimates[methods]
    output_interval <- tuning$output_projection
    on_output_boundary <- function(x) {
      x <= output_interval[1L] + 1e-12 ||
        x >= output_interval[2L] - 1e-12
    }
    all_hits <- c(
      Raw = unname(raw$projection_hit),
      BC = unname(iter2$projection_hit),
      FP = on_output_boundary(all_estimates["FP"]),
      GlobalFP = on_output_boundary(all_estimates["GlobalFP"]),
      FPID = on_output_boundary(all_estimates["FPID"]),
      LocalTotal = on_output_boundary(all_estimates["LocalTotal"]),
      PenalizedEqual = on_output_boundary(all_estimates["PenalizedEqual"]),
      PenalizedCP = on_output_boundary(all_estimates["PenalizedCP"])
    )
    hits[, type] <- all_hits[methods]
    id_diagnostics[, type] <- unlist(fp_pair[
      c("identification", "penalty", "unit_penalty")
    ])
  }
  list(
    estimate = estimates,
    projection_hit = hits,
    id_diagnostics = id_diagnostics
  )
}

# ----------------------------------------------------------------------------
# 7. Single Replication Environment
# ----------------------------------------------------------------------------

run_single_replication <- function(
    replication,
    n_values = c(250, 500, 1000, 1500, 2000),
    d_values = c(0, 0.1, 0.2, 0.3, 0.4),
    seed_base = 2026080811L,
    include_fp_sensitivity = FALSE,
    tunings = make_tuning_set(),
    psd_loading = 0.25,
    wasserstein_base_scale = 1,
    wasserstein_scale_loading = 0.25,
    output_dir = "RData_results"
) {
  # Sort and map values just as the main driver did
  n_values <- sort(unique(as.integer(n_values)))
  d_values <- sort(unique(as.numeric(d_values)))
  maximum_n <- max(n_values)
  maximum_lag <- max(vapply(tunings, function(tuning) {
    max(vapply(n_values, function(n) { floor(tuning$q * get_bandwidth(n, tuning)) }, numeric(1)))
  }, numeric(1)))
  
  fp_specs <- lapply(tunings, function(tuning) {
    setNames(lapply(n_values, function(n) make_fp_spec(n, tuning)), n_values)
  })
  
  designs <- c("PSD/Frobenius", "Wasserstein", "Scalar", "PureLocationWasserstein")
  methods <- if (include_fp_sensitivity) {
    c("Raw", "BC", "FP", "GlobalFP", "FPID", "LocalTotal", "PenalizedEqual", "PenalizedCP")
  } else {
    c("Raw", "BC", "FP")
  }
  slopes <- c("Ratio", "LogSlope")
  aggregations <- c("Baseline", "TuneAvg")
  canonical_d <- c(0.1, 0.2, 0.3, 0.4, 0)
  
  # Pre-allocate dataframe for this single replication
  total_rows <- length(d_values) * length(n_values) * length(designs) * 
    length(methods) * length(slopes) * length(aggregations)
  
  output_rep <- data.frame(
    design = character(total_rows), n = integer(total_rows), d = numeric(total_rows),
    replication = integer(total_rows), slope = character(total_rows),
    method = character(total_rows), aggregation = character(total_rows),
    estimate = numeric(total_rows), projection_hit = logical(total_rows),
    id_identification = numeric(total_rows), id_penalty = numeric(total_rows),
    id_unit_penalty = numeric(total_rows), stringsAsFactors = FALSE
  )
  output_rep$id_identification[] <- NA_real_
  output_rep$id_penalty[] <- NA_real_
  output_rep$id_unit_penalty[] <- NA_real_
  
  cursor <- 0L
  
  for (d in d_values) {
    d_seed_index <- match(d, canonical_d)
    if (is.na(d_seed_index)) stop("Each d must belong to {0,.1,.2,.3,.4}.")
    
    # Isolate RNG state for precise reproducibility per replication
    set.seed(seed_base + 100000L * d_seed_index + replication)
    phi <- runif(1L, -0.25, 0.25)
    theta <- runif(1L, -0.25, 0.25)
    
    gaussian_full <- simulate_arfima_driver(maximum_n, d, phi, theta)
    psd_objects <- generate_psd_matrices(gaussian_full, p = 3, loading = psd_loading)
    psd_full <- compute_actual_frobenius_distance(psd_objects)
    wasserstein_objects <- generate_distributional_objects(
      gaussian_full, base_scale = wasserstein_base_scale, scale_loading = wasserstein_scale_loading
    )
    wasserstein_full <- compute_actual_wasserstein_distance(wasserstein_objects)
    scalar_full <- compute_scalar_distance_matrix(gaussian_full)
    pl_wasserstein_objects <- generate_distributional_objects(gaussian_full, base_scale = 1, scale_loading = 0)
    pl_wasserstein_full <- compute_actual_wasserstein_distance(pl_wasserstein_objects)
    
    for (n in n_values) {
      inputs_by_design <- list(
        `PSD/Frobenius` = prepare_distance_inputs(psd_full, n, maximum_lag, tunings[[1L]]$block_counts),
        Wasserstein = prepare_distance_inputs(wasserstein_full, n, maximum_lag, tunings[[1L]]$block_counts),
        Scalar = prepare_distance_inputs(scalar_full, n, maximum_lag, tunings[[1L]]$block_counts),
        PureLocationWasserstein = prepare_distance_inputs(pl_wasserstein_full, n, maximum_lag, tunings[[1L]]$block_counts)
      )
      
      for (design in names(inputs_by_design)) {
        estimates_by_tuning <- array(NA_real_, dim = c(length(methods), length(slopes), length(tunings)), dimnames = list(methods, slopes, names(tunings)))
        hits_by_tuning <- array(FALSE, dim = dim(estimates_by_tuning), dimnames = dimnames(estimates_by_tuning))
        diagnostics_by_tuning <- array(NA_real_, dim = c(3L, length(slopes), length(tunings)), dimnames = list(c("identification", "penalty", "unit_penalty"), slopes, names(tunings)))
        
        for (tuning_name in names(tunings)) {
          family <- estimate_family(inputs_by_design[[design]], tunings[[tuning_name]], fp_specs[[tuning_name]][[as.character(n)]], include_fp_sensitivity = include_fp_sensitivity)
          estimates_by_tuning[, , tuning_name] <- family$estimate
          hits_by_tuning[, , tuning_name] <- family$projection_hit
          diagnostics_by_tuning[, , tuning_name] <- family$id_diagnostics
        }
        
        baseline <- estimates_by_tuning[, , "T0"]
        average <- apply(estimates_by_tuning, c(1L, 2L), mean)
        baseline_hits <- hits_by_tuning[, , "T0"]
        average_hits <- apply(hits_by_tuning, c(1L, 2L), any)
        baseline_diagnostics <- diagnostics_by_tuning[, , "T0"]
        average_diagnostics <- apply(diagnostics_by_tuning, c(1L, 2L), mean)
        
        for (aggregation in aggregations) {
          values <- if (aggregation == "Baseline") baseline else average
          hit_values <- if (aggregation == "Baseline") baseline_hits else average_hits
          diagnostic_values <- if (aggregation == "Baseline") baseline_diagnostics else average_diagnostics
          
          for (slope in slopes) {
            for (method in methods) {
              cursor <- cursor + 1L
              is_fpid <- method == "FPID"
              output_rep[cursor, ] <- list(
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
    }
  }
  
  file_path <- file.path(output_dir, sprintf("rep_%04d.RData", replication))
  save(output_rep, file = file_path)
  return(file_path)
}