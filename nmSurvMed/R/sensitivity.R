# Bandwidth sensitivity analysis ------------------------------------------

.sensitivity_bandwidth_settings <- function(candidates, selected) {
  choices <- lapply(names(selected), function(scenario) {
    unique(c(
      min(candidates),
      unname(selected[scenario]),
      max(candidates)
    ))
  })
  names(choices) <- names(selected)
  combinations <- expand.grid(
    choices,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  combinations <- unique(combinations)
  combinations$setting <- sprintf(
    "Setting %d", seq_len(nrow(combinations))
  )
  combinations[c("setting", names(selected))]
}

.sensitivity_point_estimates <- function(data, bandwidth, time, exposure,
                                         reference, columns) {
  hazard_ee <- .estimate_cf_hazard_impl(
    data, exposure, exposure,
    bandwidth["exposure_exposure"], time, columns
  )
  hazard_er <- .estimate_cf_hazard_impl(
    data, exposure, reference,
    bandwidth["exposure_reference"], time, columns
  )
  hazard_rr <- .estimate_cf_hazard_impl(
    data, reference, reference,
    bandwidth["reference_reference"], time, columns
  )
  list(
    hazards = list(
      exposure_exposure = hazard_ee,
      exposure_reference = hazard_er,
      reference_reference = hazard_rr
    ),
    survival = list(
      direct = exp(-cumsum(hazard_er)) - exp(-cumsum(hazard_rr)),
      indirect = exp(-cumsum(hazard_ee)) - exp(-cumsum(hazard_er))
    ),
    cumulative_hazard = list(
      direct = cumsum(hazard_er - hazard_rr),
      indirect = cumsum(hazard_ee - hazard_er)
    )
  )
}

.sensitivity_effect_table <- function(setting, time, point, bootstrap,
                                      scale_name) {
  direct_name <- paste0("direct_", scale_name)
  indirect_name <- paste0("indirect_", scale_name)
  data.frame(
    setting = setting,
    time = time,
    direct = point[[scale_name]]$direct,
    direct_se = bootstrap$summaries[[direct_name]]$se,
    direct_lower = bootstrap$summaries[[direct_name]]$lower,
    direct_upper = bootstrap$summaries[[direct_name]]$upper,
    indirect = point[[scale_name]]$indirect,
    indirect_se = bootstrap$summaries[[indirect_name]]$se,
    indirect_lower = bootstrap$summaries[[indirect_name]]$lower,
    indirect_upper = bootstrap$summaries[[indirect_name]]$upper,
    row.names = NULL
  )
}

#' Assess sensitivity to the kernel bandwidth
#'
#' For each of the three counterfactual hazard scenarios, constructs the
#' unique set consisting of the smallest candidate bandwidth, its
#' scenario-specific cross-validation bandwidth, and the largest candidate
#' bandwidth. It then evaluates every Cartesian-product combination (at most
#' \eqn{3^3=27}). Duplicate values are removed when an optimal bandwidth
#' equals an endpoint. Bandwidths remain fixed in every bootstrap sample, and
#' the same bootstrap resampling seeds are used for every setting.
#'
#' @inheritParams nmSurvMed
#'
#' @return An object of class `nm_bandwidth_sensitivity` containing the
#' bandwidth settings, effect estimates and pointwise bootstrap confidence
#' intervals, cross-validation details, diagnostics, and runtime.
#' @export
#'
#' @examples
#' data <- nmSurvMed:::nm_example_data(100)
#' time <- seq(0, ceiling(max(data$T2) * 10) / 10, by = 0.1)
#' sensitivity <- bandwidth_sensitivity(
#'   data,
#'   T1 = "T1_hat", D1 = "D1",
#'   T2 = "T2_hat", D2 = "D2", Z = "Z",
#'   time = time, bandwidth = c(0.2, 0.3),
#'   folds = 2, bootstrap = 3, verbose = FALSE
#' )
bandwidth_sensitivity <- function(
    data, T1, D1, T2, D2, Z, time,
    bandwidth = "paper_recommend",
    exposure = 1, reference = 0,
    scale = c("both", "survival", "cumulative_hazard"),
    folds = 5L, bootstrap = 500L, conf_level = 0.95,
    ci_type = "percentile", parallel = FALSE, cores = 1L,
    seed = 123, keep_bootstrap = FALSE, verbose = TRUE) {
  call <- match.call()
  total_start <- proc.time()[["elapsed"]]
  scale <- match.arg(scale)

  if (any(lengths(list(T1, D1, T2, D2, Z)) != 1L) ||
      !all(vapply(list(T1, D1, T2, D2, Z), is.character, logical(1L)))) {
    stop("`T1`, `D1`, `T2`, `D2`, and `Z` must be single column names.",
         call. = FALSE)
  }
  roles <- c(
    intermediate_time = T1,
    intermediate_event = D1,
    terminal_time = T2,
    terminal_event = D2,
    exposure = Z
  )
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
    message("Selecting the optimal scenario-specific bandwidths...")
  }
  selection_start <- proc.time()[["elapsed"]]
  selected <- select_bandwidth(
    checked$data, candidates, checked$time, folds,
    exposure, reference, seed, validated$cores, checked$columns
  )
  selection_seconds <- proc.time()[["elapsed"]] - selection_start
  settings <- .sensitivity_bandwidth_settings(
    candidates, selected$selected
  )

  if (verbose) {
    message(sprintf(
      "Estimating and bootstrapping %d unique bandwidth setting%s...",
      nrow(settings), if (nrow(settings) == 1L) "" else "s"
    ))
  }
  estimation_seconds <- 0
  bootstrap_seconds <- 0
  setting_results <- vector("list", nrow(settings))
  for (index in seq_len(nrow(settings))) {
    setting_bandwidth <- unlist(
      settings[index, -1L, drop = FALSE], use.names = TRUE
    )
    estimation_start <- proc.time()[["elapsed"]]
    point <- .sensitivity_point_estimates(
      checked$data, setting_bandwidth, checked$time,
      exposure, reference, checked$columns
    )
    estimation_seconds <- estimation_seconds +
      proc.time()[["elapsed"]] - estimation_start

    bootstrap_start <- proc.time()[["elapsed"]]
    boot <- .bootstrap_all_effects(
      checked$data, setting_bandwidth, checked$time,
      exposure, reference, validated$bootstrap, conf_level,
      seed, validated$cores, checked$columns
    )
    bootstrap_seconds <- bootstrap_seconds +
      proc.time()[["elapsed"]] - bootstrap_start

    setting_results[[index]] <- list(point = point, bootstrap = boot)
  }

  all_effects <- lapply(
    c("survival", "cumulative_hazard"),
    function(scale_name) {
      do.call(rbind, lapply(seq_len(nrow(settings)), function(index) {
        .sensitivity_effect_table(
          settings$setting[index], checked$time,
          setting_results[[index]]$point,
          setting_results[[index]]$bootstrap,
          scale_name
        )
      }))
    }
  )
  names(all_effects) <- c("survival", "cumulative_hazard")
  effects <- if (scale == "both") all_effects else all_effects[[scale]]

  bootstrap_output <- list(
    B = validated$bootstrap,
    conf_level = conf_level,
    ci_type = ci_type,
    bandwidth_reselected = FALSE,
    shared_resamples_across_settings = TRUE
  )
  if (keep_bootstrap) {
    bootstrap_output$replicates <- lapply(
      setting_results, function(result) result$bootstrap$replicates
    )
    names(bootstrap_output$replicates) <- settings$setting
  }

  hazards <- lapply(setting_results, function(result) result$point$hazards)
  names(hazards) <- settings$setting
  diagnostics <- .nm_diagnostics(
    checked$data, checked$time, checked$columns,
    exposure, reference, hazards
  )
  if (verbose && length(diagnostics$warnings)) {
    warning(paste(diagnostics$warnings, collapse = "\n"), call. = FALSE)
  }

  total_seconds <- proc.time()[["elapsed"]] - total_start
  runtime <- data.frame(
    stage = c(
      "Bandwidth selection", "Point estimation", "Bootstrap", "Total"
    ),
    elapsed_seconds = round(c(
      selection_seconds, estimation_seconds,
      bootstrap_seconds, total_seconds
    ), 3L),
    row.names = NULL
  )
  display_settings <- settings
  display_settings[-1L] <- lapply(
    display_settings[-1L], round, digits = 3L
  )
  names(display_settings) <- c(
    "setting", "z_1_1", "z_1_0", "z_0_0"
  )

  structure(
    list(
      call = call,
      bandwidth_settings = display_settings,
      bandwidth_values = settings,
      bandwidth_candidates = candidates,
      bandwidth_scores = selected$scores,
      effects = effects,
      bootstrap = bootstrap_output,
      counterfactual_hazards = hazards,
      time = checked$time,
      scale = scale,
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
        keep_bootstrap = keep_bootstrap
      ),
      runtime = runtime
    ),
    class = "nm_bandwidth_sensitivity"
  )
}

#' @export
print.nm_bandwidth_sensitivity <- function(x, ...) {
  cat("Bandwidth sensitivity analysis\n")
  cat(sprintf(
    "%d unique setting%s; %d bootstrap replicates\n",
    nrow(x$bandwidth_settings),
    if (nrow(x$bandwidth_settings) == 1L) "" else "s",
    x$bootstrap$B
  ))
  print(x$bandwidth_settings, row.names = FALSE)
  invisible(x)
}

.combine_conclusion_intervals <- function(pointwise, minimum_points) {
  run_id <- cumsum(c(
    TRUE,
    pointwise$conclusion[-1L] !=
      pointwise$conclusion[-nrow(pointwise)]
  ))
  rows <- split(seq_len(nrow(pointwise)), run_id)
  rows <- Filter(function(index) {
    length(index) >= minimum_points
  }, rows)
  if (!length(rows)) {
    return(data.frame(
      conclusion = character(), start_time = numeric(), end_time = numeric(),
      n_time_points = integer(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(rows, function(index) {
    data.frame(
      conclusion = pointwise$conclusion[index[1L]],
      start_time = pointwise$time[index[1L]],
      end_time = pointwise$time[index[length(index)]],
      n_time_points = length(index),
      row.names = NULL
    )
  }))
}

#' Classify conclusions across all bandwidth settings
#'
#' At each time point, confidence intervals from every bandwidth setting are
#' classified as follows:
#' \describe{
#'   \item{Positive}{Every confidence interval is strictly above zero.}
#'   \item{Negative}{Every confidence interval is strictly below zero.}
#'   \item{Not significant}{Every confidence interval contains zero.}
#'   \item{Mixed}{The bandwidth settings do not share one of the preceding
#'     conclusions; for example, some are significant and others are not.}
#' }
#' Consecutive grid points with the same conclusion are combined into time
#' intervals.
#'
#' @param x An object returned by [bandwidth_sensitivity()].
#' @param minimum_points Minimum number of consecutive time-grid points
#'   required for a reported interval.
#'
#' @return An object of class `nm_bandwidth_decision`. `summary` contains one
#'   row per scale and effect and reports only the consistently significant
#'   `Positive` and `Negative` intervals. `intervals` retains all four
#'   conclusions in long format, `pointwise` retains the classification at
#'   every grid point, and `definitions` explains all conclusion labels.
#' @export
bandwidth_decision <- function(x, minimum_points = 1L) {
  if (!inherits(x, "nm_bandwidth_sensitivity")) {
    stop(
      "`x` must be an object returned by `bandwidth_sensitivity()`.",
      call. = FALSE
    )
  }
  minimum_points <- as.integer(minimum_points)
  if (length(minimum_points) != 1L || is.na(minimum_points) ||
      minimum_points < 1L) {
    stop("`minimum_points` must be a positive integer.", call. = FALSE)
  }

  effect_tables <- if (identical(x$scale, "both")) {
    x$effects
  } else {
    stats::setNames(list(x$effects), x$scale)
  }
  pointwise <- list()
  intervals <- list()
  for (scale_name in names(effect_tables)) {
    table <- effect_tables[[scale_name]]
    settings <- unique(table$setting)
    for (effect in c("direct", "indirect")) {
      by_setting <- split(table, table$setting)
      lower <- do.call(cbind, lapply(
        by_setting, `[[`, paste0(effect, "_lower")
      ))
      upper <- do.call(cbind, lapply(
        by_setting, `[[`, paste0(effect, "_upper")
      ))
      valid_ci <- is.finite(lower) & is.finite(upper)
      positive <- rowSums(
        lower > 0 & valid_ci
      ) == length(settings)
      negative <- rowSums(
        upper < 0 & valid_ci
      ) == length(settings)
      not_significant <- rowSums(
        lower <= 0 & upper >= 0 & valid_ci
      ) == length(settings)
      conclusion <- ifelse(
        positive, "Positive",
        ifelse(
          negative, "Negative",
          ifelse(not_significant, "Not significant", "Mixed")
        )
      )
      key <- paste(scale_name, effect, sep = "_")
      pointwise[[key]] <- data.frame(
        scale = scale_name,
        effect = effect,
        time = by_setting[[1L]]$time,
        conclusion = conclusion,
        row.names = NULL
      )
      combined <- .combine_conclusion_intervals(
        pointwise[[key]], minimum_points
      )
      if (nrow(combined)) {
        combined$scale <- scale_name
        combined$effect <- effect
        combined <- combined[
          c("scale", "effect", "conclusion", "start_time", "end_time",
            "n_time_points")
        ]
      }
      intervals[[key]] <- combined
    }
  }
  pointwise <- do.call(rbind, pointwise)
  intervals <- Filter(function(value) nrow(value) > 0L, intervals)
  intervals <- if (length(intervals)) {
    do.call(rbind, intervals)
  } else {
    data.frame(
      scale = character(), effect = character(), conclusion = character(),
      start_time = numeric(), end_time = numeric(),
      n_time_points = integer(), stringsAsFactors = FALSE
    )
  }
  rownames(pointwise) <- NULL
  rownames(intervals) <- NULL
  display_scale <- function(value) {
    ifelse(value == "cumulative_hazard", "Cum. hazard", "Survival")
  }
  display_effect <- function(value) {
    ifelse(value == "direct", "Direct", "Indirect")
  }
  pointwise$scale <- display_scale(pointwise$scale)
  intervals$scale <- display_scale(intervals$scale)
  pointwise$effect <- display_effect(pointwise$effect)
  intervals$effect <- display_effect(intervals$effect)

  combinations <- unique(pointwise[c("scale", "effect")])
  conclusion_levels <- c(
    "Positive", "Negative", "Not significant", "Mixed"
  )
  format_intervals <- function(scale_name, effect_name, conclusion_name) {
    selected <- intervals[
      intervals$scale == scale_name &
        intervals$effect == effect_name &
        intervals$conclusion == conclusion_name, ,
      drop = FALSE
    ]
    if (!nrow(selected)) {
      return("None")
    }
    format_time <- function(value) {
      format(round(value, 10L), trim = TRUE, scientific = FALSE)
    }
    paste0(
      "[", vapply(selected$start_time, format_time, character(1L)),
      ", ", vapply(selected$end_time, format_time, character(1L)), "]",
      collapse = "; "
    )
  }
  summary <- do.call(rbind, lapply(seq_len(nrow(combinations)), function(i) {
    values <- vapply(conclusion_levels, function(conclusion_name) {
      format_intervals(
        combinations$scale[i],
        combinations$effect[i],
        conclusion_name
      )
    }, character(1L))
    data.frame(
      scale = combinations$scale[i],
      effect = combinations$effect[i],
      as.list(values[c("Positive", "Negative")]),
      row.names = NULL,
      check.names = FALSE
    )
  }))
  definitions <- data.frame(
    conclusion = conclusion_levels,
    definition = c(
      "All confidence intervals are above zero.",
      "All confidence intervals are below zero.",
      "All confidence intervals contain zero.",
      "Conclusions differ across bandwidth settings."
    ),
    row.names = NULL
  )

  structure(
    list(
      summary = summary,
      intervals = intervals,
      pointwise = pointwise,
      definitions = definitions,
      minimum_points = minimum_points
    ),
    class = "nm_bandwidth_decision"
  )
}

#' @export
print.nm_bandwidth_decision <- function(x, ...) {
  cat("Bandwidth decision summary\n")
  cat("Consistently significant intervals across all bandwidth settings.\n")
  cat("`None` means no such interval.\n")
  print(x$summary, row.names = FALSE)
  cat(
    "Not-significant and mixed results remain in `$intervals` and ",
    "`$pointwise`.\n",
    sep = ""
  )
  invisible(x)
}
