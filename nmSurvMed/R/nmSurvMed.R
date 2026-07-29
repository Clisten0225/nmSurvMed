# High-level analysis interface -------------------------------------------

.paper_bandwidths <- function(n) {
  constants <- c(0.2, 0.5, 0.7)
  powers <- c(2.5, 3, 3.5, 4, 4.5)
  unique(as.vector(outer(
    constants,
    powers,
    function(constant, power) constant * n^(-1 / power)
  )))
}

.nm_diagnostics <- function(data, time, columns, exposure, reference,
                            hazards) {
  exposure_values <- data[[columns["exposure"]]]
  d1 <- data[[columns["intermediate_event"]]]
  d2 <- data[[columns["terminal_event"]]]
  t1 <- data[[columns["intermediate_time"]]]
  t2 <- data[[columns["terminal_time"]]]

  warnings <- character()
  if (any(t1 < 0) || any(t2 < 0)) {
    warnings <- c(warnings, "Negative event or follow-up times were found.")
  }
  if (any(t2 < t1)) {
    warnings <- c(warnings, "Some observations have T2 < T1.")
  }
  if (max(time) > max(t2)) {
    warnings <- c(
      warnings,
      "The time grid extends beyond the largest observed terminal time."
    )
  }
  if (any(!is.finite(unlist(hazards, use.names = FALSE)))) {
    warnings <- c(warnings, "Non-finite counterfactual estimates were found.")
  }

  stratum <- interaction(
    exposure_values, d1, d2, drop = TRUE,
    sep = ":"
  )
  at_risk <- vapply(time, function(at) sum(t2 >= at), integer(1L))
  empty_risk_times <- time[at_risk == 0L]
  if (length(empty_risk_times)) {
    warnings <- c(
      warnings,
      "The risk set is empty at one or more requested time points."
    )
  }

  list(
    sample_size = nrow(data),
    exposure_counts = table(exposure_values, useNA = "ifany"),
    event_counts = c(D1 = sum(d1 == 1), D2 = sum(d2 == 1)),
    stratum_counts = table(stratum),
    empty_risk_times = empty_risk_times,
    nonfinite_estimates = any(
      !is.finite(unlist(hazards, use.names = FALSE))
    ),
    contrast_present = all(c(exposure, reference) %in% exposure_values),
    warnings = unique(warnings)
  )
}

.validate_analysis_arguments <- function(data, bandwidth, exposure, reference,
                                         bootstrap, conf_level, ci_type,
                                         parallel, cores, verbose,
                                         keep_bootstrap) {
  if (identical(exposure, reference)) {
    stop("`exposure` and `reference` must be different.", call. = FALSE)
  }
  if (!all(c(exposure, reference) %in% data)) {
    stop(
      "Both `exposure` and `reference` must occur in the exposure column.",
      call. = FALSE
    )
  }
  if (is.character(bandwidth)) {
    if (length(bandwidth) != 1L ||
        !bandwidth %in% c("paper_recommand", "paper")) {
      stop(
        paste0(
          "Character `bandwidth` must be exactly ",
          "\"paper_recommand\"."
        ),
        call. = FALSE
      )
    }
    if (bandwidth == "paper") {
      warning(
        paste0(
          "`bandwidth = \"paper\"` is deprecated; use ",
          "`bandwidth = \"paper_recommand\"`."
        ),
        call. = FALSE
      )
      bandwidth <- "paper_recommand"
    }
  } else {
    bandwidth <- as.numeric(bandwidth)
    if (!length(bandwidth)) {
      stop("At least one candidate bandwidth is required.", call. = FALSE)
    }
    .check_bandwidth(bandwidth, length(bandwidth))
  }
  bootstrap <- as.integer(bootstrap)
  if (length(bootstrap) != 1L || is.na(bootstrap) || bootstrap < 2L) {
    stop("`bootstrap` must be an integer of at least 2.", call. = FALSE)
  }
  if (length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be strictly between 0 and 1.", call. = FALSE)
  }
  if (length(ci_type) != 1L || ci_type != "percentile") {
    stop("The current version supports only `ci_type = \"percentile\"`.",
         call. = FALSE)
  }
  logical_arguments <- list(
    parallel = parallel,
    verbose = verbose,
    keep_bootstrap = keep_bootstrap
  )
  for (argument in names(logical_arguments)) {
    value <- logical_arguments[[argument]]
    if (length(value) != 1L || is.na(value) || !is.logical(value)) {
      stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    stop("`cores` must be a positive integer.", call. = FALSE)
  }
  .check_parallel_platform(if (parallel) cores else 1L)
  detected <- parallel::detectCores(logical = TRUE)
  if (parallel && is.finite(detected) && cores > detected) {
    stop(
      sprintf("`cores` (%d) exceeds the %d detected logical cores.",
              cores, detected),
      call. = FALSE
    )
  }
  list(
    bandwidth = bandwidth,
    bootstrap = bootstrap,
    cores = if (parallel) cores else 1L
  )
}

#' Fit a complete non-Markov survival mediation analysis
#'
#' Runs the full analysis pipeline: validates and maps the illness-death data,
#' selects three scenario-specific kernel bandwidths by stratified K-fold
#' cross-validation, estimates direct and indirect effect curves, and obtains
#' pointwise percentile confidence intervals from an event-stratified
#' bootstrap. The selected bandwidths are held fixed in every bootstrap
#' sample.
#'
#' @param data A data frame containing the analysis variables.
#' @param T1,D1,T2,D2,Z Single character strings naming, respectively, the
#'   intermediate-event time, intermediate-event indicator, terminal-event
#'   time, terminal-event indicator, and exposure columns in `data`.
#' @param time Equally spaced, strictly increasing evaluation grid.
#' @param bandwidth Either `"paper_recommand"` or a numeric vector of positive
#'   candidate bandwidths. `"paper_recommand"` uses
#'   \eqn{\{C m^{-1/K}: C \in \{0.2,0.5,0.7\},
#'   K \in \{2.5,3,3.5,4,4.5\}\}}, where \eqn{m} is the sample size.
#' @param exposure Exposure level of interest.
#' @param reference Reference exposure level.
#' @param scale Effect scales to return. The default `"both"` returns
#'   survival and cumulative-hazard results from the same hazard estimates and
#'   bootstrap samples. A single scale may be requested for a smaller result
#'   object.
#' @param folds Number of event-stratified cross-validation folds.
#' @param bootstrap Number of event-stratified bootstrap samples.
#' @param conf_level Pointwise confidence level.
#' @param ci_type Confidence interval method. The current version supports
#'   `"percentile"`.
#' @param parallel Logical; use parallel computation when `TRUE`. Forked
#'   processes are used on Linux and macOS; PSOCK workers are used on Windows.
#' @param cores Positive number of worker processes when `parallel = TRUE`.
#'   This cannot exceed the detected number of logical cores.
#' @param seed Random seed for fold construction and bootstrap resampling.
#' @param keep_bootstrap Logical; retain complete bootstrap replicate matrices.
#' @param verbose Logical; print progress messages.
#'
#' @return An object of class `nmSurvMed` containing:
#' \describe{
#'   \item{effects}{For `scale = "both"`, a list containing `survival` and
#'     `cumulative_hazard` data frames of direct and indirect estimates,
#'     standard errors, and pointwise confidence limits. For a single scale,
#'     the corresponding data frame.}
#'   \item{bandwidth}{A user-friendly table giving the \eqn{z_a} and
#'     \eqn{z_b} values for each counterfactual setting and its selected
#'     bandwidth rounded to three decimal places.}
#'   \item{bandwidth_values}{The same selected bandwidths as a named numeric
#'     vector for programmatic use.}
#'   \item{bandwidth_candidates}{Candidate bandwidths that were evaluated.}
#'   \item{bandwidth_scores}{Cross-validation score matrix.}
#'   \item{counterfactual_hazards}{Estimated hazard-increment curves.}
#'   \item{counterfactual_survival}{Estimated counterfactual survival curves.}
#'   \item{bootstrap}{Bootstrap summaries and, optionally, replicates.}
#'   \item{diagnostics}{Sample, event, stratum, risk-set, and warning details.}
#'   \item{settings}{Analysis and computation settings.}
#'   \item{data_roles}{A table mapping Exposure, Mediator, and Outcome roles
#'     to analysis variables and input data columns.}
#'   \item{runtime}{A table of elapsed seconds for each analysis stage,
#'     rounded to three decimal places.}
#' }
#' @export
#'
#' @examples
#' data <- nmSurvMed:::nm_example_data(100)
#' example_time <- seq(
#'   0, ceiling(max(data$T2) * 10) / 10, by = 0.1
#' )
#' fit <- nmSurvMed(
#'   data,
#'   T1 = "T1_hat", D1 = "D1",
#'   T2 = "T2_hat", D2 = "D2", Z = "Z",
#'   time = example_time,
#'   bandwidth = c(0.2, 0.3),
#'   folds = 2, bootstrap = 3,
#'   parallel = FALSE, seed = 1
#' )
nmSurvMed <- function(data, T1, D1, T2, D2, Z, time,
                      bandwidth = "paper_recommand",
                      exposure = 1, reference = 0,
                      scale = c("both", "survival", "cumulative_hazard"),
                      folds = 5L, bootstrap = 500L, conf_level = 0.95,
                      ci_type = "percentile", parallel = FALSE, cores = 1L,
                      seed = 123, keep_bootstrap = FALSE, verbose = TRUE) {
  call <- match.call()
  total_start <- proc.time()[["elapsed"]]
  scale <- match.arg(scale)

  roles <- c(
    intermediate_time = T1,
    intermediate_event = D1,
    terminal_time = T2,
    terminal_event = D2,
    exposure = Z
  )
  if (any(lengths(list(T1, D1, T2, D2, Z)) != 1L) ||
      !all(vapply(list(T1, D1, T2, D2, Z), is.character, logical(1L)))) {
    stop("`T1`, `D1`, `T2`, `D2`, and `Z` must be single column names.",
         call. = FALSE)
  }
  checked <- .validate_nm_data(data, time, roles)
  validated <- .validate_analysis_arguments(
    checked$data[[checked$columns["exposure"]]],
    bandwidth, exposure, reference, bootstrap, conf_level, ci_type,
    parallel, cores, verbose, keep_bootstrap
  )
  candidates <- if (is.character(validated$bandwidth)) {
    .paper_bandwidths(nrow(checked$data))
  } else {
    unique(validated$bandwidth)
  }
  parallel_backend <- .parallel_backend(validated$cores)

  if (verbose) {
    message(sprintf(
      "Parallel backend: %s (%d worker%s)",
      if (parallel_backend == "sequential") {
        "Sequential"
      } else {
        toupper(parallel_backend)
      },
      validated$cores,
      if (validated$cores == 1L) "" else "s"
    ))
    message(
      sprintf(
        paste0(
          "Selecting bandwidths from %d candidate%s using %d-fold CV%s..."
        ),
        length(candidates), if (length(candidates) == 1L) "" else "s", folds
        ,
        if (parallel) sprintf(" on %d cores", validated$cores) else ""
      )
    )
  }
  bandwidth_start <- proc.time()[["elapsed"]]
  bw <- select_bandwidth(
    checked$data,
    candidates = candidates,
    time = checked$time,
    folds = folds,
    exposure = exposure,
    reference = reference,
    seed = seed,
    cores = validated$cores,
    columns = checked$columns
  )
  bandwidth_table <- .bandwidth_table(bw$selected, bw$scenarios)
  bandwidth_seconds <- proc.time()[["elapsed"]] - bandwidth_start

  if (verbose) {
    message("Selected bandwidths:")
    message(paste(
      utils::capture.output(print(bandwidth_table, row.names = FALSE)),
      collapse = "\n"
    ))
    message("Estimating direct and indirect effects...")
  }
  estimation_start <- proc.time()[["elapsed"]]
  hazard_ee <- estimate_cf_hazard(
    checked$data, exposure, exposure,
    bw$selected["exposure_exposure"], checked$time, checked$columns
  )$estimate
  hazard_er <- estimate_cf_hazard(
    checked$data, exposure, reference,
    bw$selected["exposure_reference"], checked$time, checked$columns
  )$estimate
  hazard_rr <- estimate_cf_hazard(
    checked$data, reference, reference,
    bw$selected["reference_reference"], checked$time, checked$columns
  )$estimate
  hazards <- list(
    exposure_exposure = hazard_ee,
    exposure_reference = hazard_er,
    reference_reference = hazard_rr
  )
  survival <- lapply(hazards, function(value) exp(-cumsum(value)))
  point_estimates <- list(
    survival = list(
      direct = .hazard_to_scale(hazard_er, hazard_rr, "survival"),
      indirect = .hazard_to_scale(hazard_ee, hazard_er, "survival")
    ),
    cumulative_hazard = list(
      direct = .hazard_to_scale(
        hazard_er, hazard_rr, "cumulative_hazard"
      ),
      indirect = .hazard_to_scale(
        hazard_ee, hazard_er, "cumulative_hazard"
      )
    )
  )
  estimation_seconds <- proc.time()[["elapsed"]] - estimation_start

  if (verbose) {
    message(
      sprintf(
        "Running %d stratified bootstrap replicates%s...",
        validated$bootstrap,
        if (parallel) sprintf(" on %d cores", validated$cores) else ""
      )
    )
  }
  bootstrap_start <- proc.time()[["elapsed"]]
  all_bootstrap <- .bootstrap_all_effects(
    checked$data,
    bandwidth = bw$selected,
    time = checked$time,
    exposure = exposure,
    reference = reference,
    B = validated$bootstrap,
    conf_level = conf_level,
    seed = seed,
    cores = validated$cores,
    columns = checked$columns
  )
  bootstrap_seconds <- proc.time()[["elapsed"]] - bootstrap_start

  make_effect_table <- function(scale_name) {
    direct_name <- paste0("direct_", scale_name)
    indirect_name <- paste0("indirect_", scale_name)
    data.frame(
      time = checked$time,
      direct = point_estimates[[scale_name]]$direct,
      direct_se = all_bootstrap$summaries[[direct_name]]$se,
      direct_lower = all_bootstrap$summaries[[direct_name]]$lower,
      direct_upper = all_bootstrap$summaries[[direct_name]]$upper,
      indirect = point_estimates[[scale_name]]$indirect,
      indirect_se = all_bootstrap$summaries[[indirect_name]]$se,
      indirect_lower = all_bootstrap$summaries[[indirect_name]]$lower,
      indirect_upper = all_bootstrap$summaries[[indirect_name]]$upper
    )
  }
  all_effects <- list(
    survival = make_effect_table("survival"),
    cumulative_hazard = make_effect_table("cumulative_hazard")
  )
  effects <- if (scale == "both") all_effects else all_effects[[scale]]

  bootstrap_output <- list(
    summaries = all_bootstrap$summaries,
    B = validated$bootstrap,
    conf_level = conf_level,
    ci_type = ci_type,
    bandwidth_reselected = FALSE
  )
  if (keep_bootstrap) {
    bootstrap_output$replicates <- all_bootstrap$replicates
  }

  diagnostics <- .nm_diagnostics(
    checked$data, checked$time, checked$columns,
    exposure, reference, hazards
  )
  if (verbose && length(diagnostics$warnings)) {
    warning(
      paste(diagnostics$warnings, collapse = "\n"),
      call. = FALSE
    )
  }
  total_seconds <- proc.time()[["elapsed"]] - total_start
  data_roles <- data.frame(
    component = c(
      "Exposure", "Mediator", "Mediator", "Outcome", "Outcome"
    ),
    analysis_variable = c("Z", "T1", "D1", "T2", "D2"),
    data_column = c(Z, T1, D1, T2, D2),
    meaning = c(
      "Exposure indicator",
      "Mediator event time",
      "Mediator event status",
      "Outcome event time",
      "Outcome event status"
    ),
    row.names = NULL
  )
  runtime <- data.frame(
    stage = c(
      "Bandwidth selection",
      "Point estimation",
      "Bootstrap",
      "Total"
    ),
    elapsed_seconds = round(
      c(
        bandwidth_seconds,
        estimation_seconds,
        bootstrap_seconds,
        total_seconds
      ),
      3L
    ),
    row.names = NULL
  )
  structure(
    list(
      call = call,
      effects = effects,
      bandwidth = bandwidth_table,
      bandwidth_values = bw$selected,
      bandwidth_candidates = candidates,
      bandwidth_scores = bw$scores,
      counterfactual_hazards = hazards,
      counterfactual_survival = survival,
      bootstrap = bootstrap_output,
      time = checked$time,
      scale = scale,
      data_roles = data_roles,
      contrast = list(exposure = exposure, reference = reference),
      diagnostics = diagnostics,
      settings = list(
        folds = as.integer(folds),
        bootstrap = validated$bootstrap,
        conf_level = conf_level,
        ci_type = ci_type,
        parallel = parallel,
        cores = validated$cores,
        parallel_backend = parallel_backend,
        seed = seed,
        keep_bootstrap = keep_bootstrap,
        bandwidth_rule = if (is.character(bandwidth)) {
          "paper_recommand"
        } else {
          "user"
        }
      ),
      runtime = runtime,
      session_info = list(
        package_version = utils::packageVersion("nmSurvMed"),
        R_version = R.version.string,
        platform = R.version$platform,
        date = Sys.Date()
      )
    ),
    class = "nmSurvMed"
  )
}

#' @export
print.nmSurvMed <- function(x, ...) {
  cat("Non-Markov survival mediation analysis\n")
  cat(sprintf("Scale: %s\n", x$scale))
  cat(sprintf(
    "Contrast: %s versus %s\n",
    as.character(x$contrast$exposure),
    as.character(x$contrast$reference)
  ))
  cat(sprintf(
    "Bootstrap: %d replicates, %.1f%% pointwise percentile intervals\n",
    x$bootstrap$B, 100 * x$bootstrap$conf_level
  ))
  cat("Selected bandwidths:\n")
  print(x$bandwidth, row.names = FALSE)
  if (length(x$diagnostics$warnings)) {
    cat(sprintf(
      "Diagnostics: %d warning%s; inspect `x$diagnostics$warnings`.\n",
      length(x$diagnostics$warnings),
      if (length(x$diagnostics$warnings) == 1L) "" else "s"
    ))
  }
  invisible(x)
}
