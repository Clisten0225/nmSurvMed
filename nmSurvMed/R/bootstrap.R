# Bootstrap inference ------------------------------------------------------

.stratified_resample <- function(data, columns) {
  group <- paste0(
    data[[columns["intermediate_event"]]], "_",
    data[[columns["terminal_event"]]]
  )
  indices <- unlist(
    lapply(unique(group), function(level) {
      index <- which(group == level)
      sample(index, length(index), replace = TRUE)
    }),
    use.names = FALSE
  )
  data[indices, , drop = FALSE]
}

.bootstrap_statistic <- function(data, effect, bandwidth, time, exposure,
                                 reference, scale, z_a, z_b, columns) {
  sample_data <- .stratified_resample(data, columns)
  switch(
    effect,
    cf_hazard = .estimate_cf_hazard_impl(
      sample_data, z_a, z_b, bandwidth, time, columns
    ),
    direct = estimate_direct_effect(
      sample_data, bandwidth, time, exposure, reference, scale, columns
    )$estimate,
    indirect = estimate_indirect_effect(
      sample_data, bandwidth, time, exposure, reference, scale, columns
    )$estimate
  )
}

.bootstrap_all_effects <- function(data, bandwidth, time, exposure, reference,
                                   B, conf_level, seed, cores, columns) {
  one_replicate <- function(replicate_seed) {
    set.seed(replicate_seed)
    sample_data <- .stratified_resample(data, columns)
    hazard_ee <- .estimate_cf_hazard_impl(
      sample_data, exposure, exposure,
      bandwidth["exposure_exposure"], time, columns
    )
    hazard_er <- .estimate_cf_hazard_impl(
      sample_data, exposure, reference,
      bandwidth["exposure_reference"], time, columns
    )
    hazard_rr <- .estimate_cf_hazard_impl(
      sample_data, reference, reference,
      bandwidth["reference_reference"], time, columns
    )
    list(
      direct_cumulative_hazard = cumsum(hazard_er - hazard_rr),
      indirect_cumulative_hazard = cumsum(hazard_ee - hazard_er),
      direct_survival =
        exp(-cumsum(hazard_er)) - exp(-cumsum(hazard_rr)),
      indirect_survival =
        exp(-cumsum(hazard_ee)) - exp(-cumsum(hazard_er))
    )
  }

  seeds <- seed + seq_len(B)
  results <- .parallel_lapply(seeds, one_replicate, cores)

  replicate_names <- names(results[[1L]])
  replicates <- lapply(replicate_names, function(name) {
    do.call(rbind, lapply(results, `[[`, name))
  })
  names(replicates) <- replicate_names

  alpha <- (1 - conf_level) / 2
  summaries <- lapply(replicates, function(values) {
    list(
      se = apply(values, 2L, stats::sd, na.rm = TRUE),
      lower = apply(
        values, 2L, stats::quantile,
        probs = alpha, na.rm = TRUE, names = FALSE
      ),
      upper = apply(
        values, 2L, stats::quantile,
        probs = 1 - alpha, na.rm = TRUE, names = FALSE
      )
    )
  })

  list(replicates = replicates, summaries = summaries)
}

#' Stratified bootstrap inference
#'
#' Resamples observations with replacement within `(D1,D2)` event strata,
#' re-estimates the requested curve, and constructs pointwise percentile
#' confidence intervals.
#'
#' @inheritParams estimate_direct_effect
#' @param effect Curve to bootstrap: `"direct"`, `"indirect"`, or
#'   `"cf_hazard"`.
#' @param bandwidth One positive bandwidth for `effect = "cf_hazard"`.
#'   For a direct or indirect effect, supply one shared bandwidth or two
#'   bandwidths in the order documented by [estimate_direct_effect()] or
#'   [estimate_indirect_effect()].
#' @param B Number of bootstrap samples.
#' @param conf_level Pointwise confidence level.
#' @param seed Random seed.
#' @param cores Number of parallel worker processes. Forked processes are used
#'   on Linux and macOS; PSOCK workers are used on Windows.
#' @param z_a,z_b Joint-intervention exposure values used only when
#'   `effect = "cf_hazard"`.
#'
#' @return An object of class `nm_bootstrap` containing the original estimate,
#'   bootstrap replicates, bootstrap standard errors, and pointwise percentile
#'   confidence limits.
#' @export
#'
#' @examples
#' boot <- bootstrap_nm(
#'   nmSurvMed:::nm_example_data(), effect = "direct", bandwidth = 0.25,
#'   time = seq(0, 2, by = 0.1), B = 10, seed = 1
#' )
bootstrap_nm <- function(data, effect = c("direct", "indirect", "cf_hazard"),
                         bandwidth, time, B = 500L, conf_level = 0.95,
                         exposure = 1, reference = 0,
                         scale = c("cumulative_hazard", "survival"),
                         z_a = exposure, z_b = reference, seed = 123,
                         cores = 1L, columns = NULL) {
  checked <- .validate_nm_data(data, time, columns)
  effect <- match.arg(effect)
  scale <- match.arg(scale)
  bandwidth <- .check_bandwidth(
    bandwidth, if (effect == "cf_hazard") 1L else 2L
  )
  B <- as.integer(B)
  cores <- as.integer(cores)
  if (length(B) != 1L || is.na(B) || B < 2L) {
    stop("`B` must be an integer of at least 2.", call. = FALSE)
  }
  if (length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be strictly between 0 and 1.", call. = FALSE)
  }
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    stop("`cores` must be a positive integer.", call. = FALSE)
  }
  .check_parallel_platform(cores)

  original <- switch(
    effect,
    cf_hazard = estimate_cf_hazard(
      checked$data, z_a, z_b, bandwidth, checked$time, checked$columns
    ),
    direct = estimate_direct_effect(
      checked$data, bandwidth, checked$time, exposure, reference, scale,
      checked$columns
    ),
    indirect = estimate_indirect_effect(
      checked$data, bandwidth, checked$time, exposure, reference, scale,
      checked$columns
    )
  )
  seeds <- seed + seq_len(B)
  one_replicate <- function(replicate_seed) {
    set.seed(replicate_seed)
    .bootstrap_statistic(
      checked$data, effect, bandwidth, checked$time, exposure, reference,
      scale, z_a, z_b, checked$columns
    )
  }
  replicates <- .parallel_lapply(seeds, one_replicate, cores)
  replicates <- do.call(rbind, replicates)
  alpha <- (1 - conf_level) / 2
  lower <- apply(
    replicates, 2L, stats::quantile, probs = alpha, na.rm = TRUE,
    names = FALSE
  )
  upper <- apply(
    replicates, 2L, stats::quantile, probs = 1 - alpha, na.rm = TRUE,
    names = FALSE
  )

  structure(
    list(
      estimate = original,
      replicates = replicates,
      std_error = apply(replicates, 2L, stats::sd, na.rm = TRUE),
      conf_low = lower,
      conf_high = upper,
      conf_level = conf_level,
      B = B,
      effect = effect
    ),
    class = "nm_bootstrap"
  )
}

#' @export
print.nm_bootstrap <- function(x, ...) {
  cat(
    sprintf(
      "Stratified bootstrap for %s effect (%d replicates, %.1f%% CI)\n",
      x$effect, x$B, 100 * x$conf_level
    )
  )
  invisible(x)
}
