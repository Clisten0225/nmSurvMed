# Bandwidth selection ------------------------------------------------------

.bandwidth_table <- function(selected, scenarios) {
  data.frame(
    z_a = scenarios[, 1L],
    z_b = scenarios[, 2L],
    selected_bandwidth = round(unname(selected), 3L),
    row.names = NULL,
    check.names = FALSE
  )
}

.stratified_folds <- function(data, folds, columns, seed) {
  folds <- as.integer(folds)
  if (length(folds) != 1L || is.na(folds) || folds < 2L ||
      folds > nrow(data)) {
    stop("`folds` must be between 2 and the number of observations.",
         call. = FALSE)
  }
  group <- paste0(
    data[[columns["intermediate_event"]]], "_",
    data[[columns["terminal_event"]]]
  )
  sizes <- table(group)
  if (any(sizes < folds)) {
    stop(
      "Every event stratum must contain at least `folds` observations.",
      call. = FALSE
    )
  }
  set.seed(seed)
  result <- vector("list", folds)
  for (level in unique(group)) {
    indices <- sample(which(group == level))
    group_folds <- split(
      indices,
      cut(seq_along(indices), breaks = folds, labels = FALSE)
    )
    for (k in seq_len(folds)) {
      result[[k]] <- c(result[[k]], group_folds[[k]])
    }
  }
  result
}

.cv_log_likelihood <- function(data, fold_indices, bandwidth, time,
                               z_a, z_b, columns) {
  scores <- numeric(length(fold_indices))
  for (i in seq_along(fold_indices)) {
    validation_indices <- fold_indices[[i]]
    train <- data[-validation_indices, , drop = FALSE]
    validation <- data[validation_indices, , drop = FALSE]
    increments <- .estimate_cf_hazard_impl(
      train, z_a, z_b, bandwidth, time, columns
    )
    survival <- exp(-cumsum(increments))
    probability <- c(0, -diff(survival))
    log_likelihood <- log(probability + 1e-8)

    event_time <- validation[[columns["terminal_time"]]]
    dt <- time[2L] - time[1L]
    # Preserve the indexing rule used in the paper analysis code:
    # terminal times are rounded to one decimal place before conversion to
    # positions on the time grid. Numeric indexing is intentional here.
    index <- round(event_time, 1L) / dt + 1
    validation_log_likelihood <- log_likelihood[index]
    if (anyNA(validation_log_likelihood)) {
      stop(
        "The time grid does not cover all validation terminal times.",
        call. = FALSE
      )
    }
    scores[i] <- sum(validation_log_likelihood)
  }
  sum(scores)
}

#' Select kernel bandwidths by stratified cross-validation
#'
#' Evaluates candidate bandwidths for the three counterfactual scenarios
#' needed by the natural direct and indirect effects: `(exposure, exposure)`,
#' `(exposure, reference)`, and `(reference, reference)`.
#'
#' @inheritParams estimate_direct_effect
#' @param candidates Positive candidate bandwidths.
#' @param folds Number of stratified cross-validation folds.
#' @param seed Random seed used to construct folds.
#' @param cores Number of parallel worker processes. Forked processes are used
#'   on Linux and macOS; PSOCK workers are used on Windows.
#'
#' @return An object of class `nm_bandwidth`. Its `selected` component gives
#'   the maximizing bandwidth for each scenario and `scores` contains all
#'   cross-validation scores.
#' @export
#'
#' @examples
#' example_data <- nmSurvMed:::nm_example_data()
#' example_time <- seq(
#'   0, ceiling(max(example_data$T2) * 10) / 10, by = 0.1
#' )
#' select_bandwidth(
#'   example_data, candidates = c(0.2, 0.3),
#'   time = example_time, folds = 2
#' )
select_bandwidth <- function(data, candidates, time, folds = 5L,
                             exposure = 1, reference = 0, seed = 123,
                             cores = 1L, columns = NULL) {
  checked <- .validate_nm_data(data, time, columns)
  candidates <- unique(as.numeric(candidates))
  .check_bandwidth(candidates, length(candidates))
  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    stop("`cores` must be a positive integer.", call. = FALSE)
  }
  .check_parallel_platform(cores)
  fold_indices <- .stratified_folds(
    checked$data, folds, checked$columns, seed
  )
  scenarios <- rbind(
    exposure_exposure = c(exposure, exposure),
    exposure_reference = c(exposure, reference),
    reference_reference = c(reference, reference)
  )

  task_grid <- expand.grid(
    scenario = seq_len(nrow(scenarios)),
    candidate = seq_along(candidates),
    KEEP.OUT.ATTRS = FALSE
  )
  evaluate_task <- function(task_index) {
    task <- task_grid[task_index, ]
    scenario <- scenarios[task$scenario, ]
    .cv_log_likelihood(
      data = checked$data,
      fold_indices = fold_indices,
      bandwidth = candidates[task$candidate],
      time = checked$time,
      z_a = scenario[1L],
      z_b = scenario[2L],
      columns = checked$columns
    )
  }

  task_indices <- seq_len(nrow(task_grid))
  result <- .parallel_lapply(task_indices, evaluate_task, cores)
  scores <- matrix(
    unlist(result, use.names = FALSE),
    nrow = nrow(scenarios),
    ncol = length(candidates)
  )
  rownames(scores) <- rownames(scenarios)
  colnames(scores) <- format(candidates, trim = TRUE)
  selected <- candidates[max.col(scores, ties.method = "first")]
  names(selected) <- rownames(scenarios)

  structure(
    list(
      selected = selected,
      scores = scores,
      candidates = candidates,
      folds = fold_indices,
      scenarios = scenarios
    ),
    class = "nm_bandwidth"
  )
}

#' @export
print.nm_bandwidth <- function(x, ...) {
  cat("Selected non-Markov survival mediation bandwidths:\n")
  print(.bandwidth_table(x$selected, x$scenarios), row.names = FALSE)
  invisible(x)
}
