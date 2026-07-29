# Internal utilities -------------------------------------------------------

.nm_columns <- function(columns = NULL) {
  defaults <- c(
    intermediate_time = "T1_hat",
    terminal_time = "T2_hat",
    intermediate_event = "D1",
    terminal_event = "D2",
    exposure = "Z"
  )

  if (is.null(columns)) {
    return(defaults)
  }
  if (is.null(names(columns)) || any(names(columns) == "")) {
    stop("`columns` must be a named character vector.", call. = FALSE)
  }
  unknown <- setdiff(names(columns), names(defaults))
  if (length(unknown)) {
    stop(
      "Unknown names in `columns`: ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  defaults[names(columns)] <- columns
  defaults
}

.validate_nm_data <- function(data, time, columns = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  columns <- .nm_columns(columns)
  missing_columns <- setdiff(unname(columns), names(data))
  if (length(missing_columns)) {
    stop(
      "Missing required columns: ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  time <- as.numeric(time)
  if (length(time) < 2L || anyNA(time) || any(!is.finite(time)) ||
      is.unsorted(time, strictly = TRUE)) {
    stop("`time` must contain at least two finite, strictly increasing values.",
         call. = FALSE)
  }
  increments <- diff(time)
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(time)))
  if (max(abs(increments - increments[1L])) > tolerance) {
    stop("`time` must be an equally spaced grid.", call. = FALSE)
  }

  event_columns <- unname(columns[c("intermediate_event", "terminal_event")])
  for (column in event_columns) {
    values <- data[[column]]
    if (anyNA(values) || !all(values %in% c(0, 1))) {
      stop("Event column `", column, "` must contain only 0 and 1.",
           call. = FALSE)
    }
  }
  numeric_columns <- unname(columns[c("intermediate_time", "terminal_time")])
  for (column in numeric_columns) {
    values <- data[[column]]
    if (!is.numeric(values) || anyNA(values) || any(!is.finite(values))) {
      stop("Time column `", column, "` must be finite and numeric.",
           call. = FALSE)
    }
  }
  if (nrow(data) == 0L) {
    stop("`data` must contain at least one row.", call. = FALSE)
  }

  list(data = data, time = time, columns = columns)
}

.check_bandwidth <- function(bandwidth, n = 1L) {
  bandwidth <- as.numeric(bandwidth)
  if (length(bandwidth) == 1L && n > 1L) {
    bandwidth <- rep(bandwidth, n)
  }
  if (length(bandwidth) != n || anyNA(bandwidth) ||
      any(!is.finite(bandwidth)) || any(bandwidth <= 0)) {
    stop(
      sprintf("`bandwidth` must contain %d positive finite value(s).", n),
      call. = FALSE
    )
  }
  bandwidth
}

.safe_ratio <- function(numerator, denominator) {
  if (!identical(dim(numerator), dim(denominator)) ||
      length(numerator) != length(denominator)) {
    stop("Numerator and denominator must have matching dimensions.",
         call. = FALSE)
  }
  result <- numeric(length(numerator))
  valid <- is.finite(numerator) & is.finite(denominator) & denominator != 0
  result[valid] <- numerator[valid] / denominator[valid]
  result[!is.finite(result)] <- 0
  dim(result) <- dim(numerator)
  result
}

.exposure_rows <- function(data, exposure, value) {
  data[[exposure]] == value
}

.gaussian_kernel <- function(x) stats::dnorm(x)

.effect_result <- function(estimate, time, scale, effect, contrast) {
  structure(
    list(
      time = time,
      estimate = as.numeric(estimate),
      scale = scale,
      effect = effect,
      contrast = contrast
    ),
    class = "nm_effect"
  )
}

.check_parallel_platform <- function(cores) {
  invisible(TRUE)
}

.parallel_backend <- function(cores = 1L) {
  if (as.integer(cores) <= 1L) {
    return("sequential")
  }
  backend <- getOption(
    "nmSurvMed.parallel_backend",
    if (.Platform$OS.type == "windows") "psock" else "fork"
  )
  match.arg(backend, c("fork", "psock"))
}

.parallel_lapply <- function(X, FUN, cores = 1L,
                             backend = c("auto", "fork", "psock")) {
  backend <- match.arg(backend)
  workers <- min(as.integer(cores), length(X))
  if (workers <= 1L || !length(X)) {
    return(lapply(X, FUN))
  }
  if (backend == "auto") {
    backend <- .parallel_backend(workers)
  }
  if (backend == "fork") {
    return(parallel::mclapply(
      X, FUN, mc.cores = workers, mc.set.seed = FALSE
    ))
  }

  cluster <- parallel::makeCluster(workers, type = "PSOCK")
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  # Export package functions so PSOCK also works with devtools::load_all().
  namespace <- environment(.parallel_lapply)
  namespace_objects <- mget(
    ls(namespace, all.names = TRUE),
    envir = namespace,
    inherits = FALSE
  )
  function_names <- names(namespace_objects)[vapply(
    namespace_objects, is.function, logical(1L)
  )]
  parallel::clusterExport(
    cluster, function_names, envir = namespace
  )
  parallel::parLapply(cluster, X, FUN)
}
