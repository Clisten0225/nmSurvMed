# Counterfactual hazard estimation ----------------------------------------

.v1 <- function(t, t1, h, data, z, columns) {
  selected <- .exposure_rows(data, columns["exposure"], z) &
    data[[columns["terminal_time"]]] >= t &
    data[[columns["intermediate_event"]]] == 1 &
    data[[columns["terminal_event"]]] == 1
  values <- data[[columns["intermediate_time"]]][selected]
  if (!length(values)) return(rep(0, length(t1)))
  rowSums(vapply(
    values,
    function(value) .gaussian_kernel((t1 - value) / h) / (nrow(data) * h),
    numeric(length(t1))
  ))
}

.dlambda1_denominator <- function(t, t1, h, data, z, columns) {
  selected <- .exposure_rows(data, columns["exposure"], z) &
    data[[columns["terminal_time"]]] >= t &
    data[[columns["intermediate_event"]]] == 1
  values <- data[[columns["intermediate_time"]]][selected]
  if (!length(values)) return(rep(0, length(t1)))
  rowSums(vapply(
    values,
    function(value) .gaussian_kernel((t1 - value) / h) / (nrow(data) * h),
    numeric(length(t1))
  ))
}

.dlambda1 <- function(t, t1, h, data, z, dt, columns) {
  current <- .v1(t, t1, h, data, z, columns)
  next_value <- .v1(t + dt, t1, h, data, z, columns)
  current[t <= t1] <- 0
  next_value[t <= t1] <- 0
  .safe_ratio(
    current - next_value,
    .dlambda1_denominator(t, t1, h, data, z, columns)
  )
}

.omega1 <- function(t, t1, h, data, z, columns) {
  selected <- .exposure_rows(data, columns["exposure"], z) &
    data[[columns["terminal_time"]]] >= t &
    data[[columns["intermediate_event"]]] == 1
  values <- data[[columns["intermediate_time"]]][selected]
  numerator <- if (!length(values)) {
    rep(0, length(t1))
  } else {
    rowSums(vapply(
      values,
      function(value) .gaussian_kernel((t1 - value) / h) / (nrow(data) * h),
      numeric(length(t1))
    ))
  }
  denominator <- sum(
    .exposure_rows(data, columns["exposure"], z) &
      data[[columns["terminal_time"]]] >= t
  ) / nrow(data)
  .safe_ratio(numerator, rep(denominator, length(numerator)))
}

.dlambda0 <- function(t, data, z, dt, columns) {
  at_risk <- sum(
    .exposure_rows(data, columns["exposure"], z) &
      data[[columns["intermediate_time"]]] >= t
  ) / nrow(data)
  events <- function(at) {
    sum(
      .exposure_rows(data, columns["exposure"], z) &
        data[[columns["terminal_time"]]] >= at &
        data[[columns["intermediate_event"]]] == 0 &
        data[[columns["terminal_event"]]] == 1
    ) / nrow(data)
  }
  .safe_ratio(events(t) - events(t + dt), at_risk)
}

.omega0 <- function(t, data, z, columns) {
  numerator <- sum(
    .exposure_rows(data, columns["exposure"], z) &
      data[[columns["intermediate_time"]]] >= t
  ) / nrow(data)
  denominator <- sum(
    .exposure_rows(data, columns["exposure"], z) &
      data[[columns["terminal_time"]]] >= t
  ) / nrow(data)
  .safe_ratio(numerator, denominator)
}

.estimate_cf_hazard_impl <- function(data, z_a, z_b, bandwidth, time, columns) {
  dt <- time[2L] - time[1L]
  dlambda0 <- vapply(
    time, .dlambda0, numeric(1L),
    data = data, z = z_a, dt = dt, columns = columns
  )
  omega0 <- vapply(
    time, .omega0, numeric(1L),
    data = data, z = z_b, columns = columns
  )
  part0 <- dlambda0 * omega0

  dlambda1 <- vapply(
    time, .dlambda1, numeric(length(time)),
    t1 = time, h = bandwidth, data = data, z = z_a, dt = dt,
    columns = columns
  )
  omega1 <- vapply(
    time, .omega1, numeric(length(time)),
    t1 = time, h = bandwidth, data = data, z = z_b, columns = columns
  )
  part1 <- dlambda1 * omega1
  part1[lower.tri(part1, diag = TRUE)] <- 0
  part1[!is.finite(part1)] <- 0

  as.numeric(part0 + colSums(part1) * dt)
}

#' Estimate a counterfactual hazard
#'
#' Estimates the hazard increment curve under the joint intervention
#' \eqn{(z_a, z_b)} in a non-Markov illness-death model. The transition hazard
#' under `z_a` is combined with the distribution of the intermediate event
#' under `z_b`.
#'
#' @param data A data frame containing intermediate and terminal event times,
#'   event indicators, and exposure.
#' @param z_a Exposure value used for the transition hazards.
#' @param z_b Exposure value used for the intermediate-event distribution.
#' @param bandwidth Positive kernel bandwidth.
#' @param time Equally spaced, strictly increasing evaluation grid.
#' @param columns Optional named character vector mapping the roles
#'   `intermediate_time`, `terminal_time`, `intermediate_event`,
#'   `terminal_event`, and `exposure` to columns in `data`.
#'
#' @return An object of class `nm_effect` with components `time`, `estimate`,
#'   `scale`, `effect`, and `contrast`. Here `estimate` contains hazard
#'   increments, not a cumulative hazard.
#' @export
#'
#' @examples
#' fit <- estimate_cf_hazard(
#'   nmSurvMed:::nm_example_data(), z_a = 1, z_b = 0,
#'   bandwidth = 0.25, time = seq(0, 2, by = 0.1)
#' )
estimate_cf_hazard <- function(data, z_a, z_b, bandwidth, time,
                               columns = NULL) {
  checked <- .validate_nm_data(data, time, columns)
  bandwidth <- .check_bandwidth(bandwidth)
  estimate <- .estimate_cf_hazard_impl(
    checked$data, z_a, z_b, bandwidth, checked$time, checked$columns
  )
  .effect_result(
    estimate, checked$time, "hazard_increment", "counterfactual_hazard",
    c(z_a = z_a, z_b = z_b)
  )
}
