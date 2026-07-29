# Direct and indirect effects ---------------------------------------------

.hazard_to_scale <- function(hazard1, hazard0, scale) {
  scale <- match.arg(scale, c("cumulative_hazard", "survival"))
  if (scale == "cumulative_hazard") {
    return(cumsum(hazard1 - hazard0))
  }
  exp(-cumsum(hazard1)) - exp(-cumsum(hazard0))
}

#' Estimate a non-Markov direct effect
#'
#' Compares the counterfactual curves for \eqn{(z_a,z_b)=(1,0)} and
#' \eqn{(0,0)} by default. Two bandwidths may be supplied so each
#' counterfactual curve can use its selected bandwidth.
#'
#' @inheritParams estimate_cf_hazard
#' @param bandwidth One positive bandwidth used for both curves, or two
#'   bandwidths corresponding to `(exposure, reference)` and
#'   `(reference, reference)`.
#' @param exposure Exposure level being compared.
#' @param reference Reference exposure level.
#' @param scale Effect scale: `"cumulative_hazard"` or `"survival"`.
#'
#' @return An `nm_effect` object containing the time-varying direct effect.
#' @export
#'
#' @examples
#' estimate_direct_effect(
#'   nmSurvMed:::nm_example_data(), bandwidth = 0.25,
#'   time = seq(0, 2, by = 0.1)
#' )
estimate_direct_effect <- function(data, bandwidth, time, exposure = 1,
                                   reference = 0,
                                   scale = c("cumulative_hazard", "survival"),
                                   columns = NULL) {
  checked <- .validate_nm_data(data, time, columns)
  bandwidth <- .check_bandwidth(bandwidth, 2L)
  scale <- match.arg(scale)
  treated <- .estimate_cf_hazard_impl(
    checked$data, exposure, reference, bandwidth[1L],
    checked$time, checked$columns
  )
  control <- .estimate_cf_hazard_impl(
    checked$data, reference, reference, bandwidth[2L],
    checked$time, checked$columns
  )
  .effect_result(
    .hazard_to_scale(treated, control, scale), checked$time, scale,
    "direct", c(exposure = exposure, reference = reference)
  )
}

#' Estimate a non-Markov indirect effect
#'
#' Compares the counterfactual curves for \eqn{(z_a,z_b)=(1,1)} and
#' \eqn{(1,0)} by default. Two bandwidths may be supplied so each
#' counterfactual curve can use its selected bandwidth.
#'
#' @inheritParams estimate_direct_effect
#' @param bandwidth One positive bandwidth used for both curves, or two
#'   bandwidths corresponding to `(exposure, exposure)` and
#'   `(exposure, reference)`.
#' @return An `nm_effect` object containing the time-varying indirect effect.
#' @export
#'
#' @examples
#' estimate_indirect_effect(
#'   nmSurvMed:::nm_example_data(), bandwidth = 0.25,
#'   time = seq(0, 2, by = 0.1), scale = "survival"
#' )
estimate_indirect_effect <- function(data, bandwidth, time, exposure = 1,
                                     reference = 0,
                                     scale = c("cumulative_hazard", "survival"),
                                     columns = NULL) {
  checked <- .validate_nm_data(data, time, columns)
  bandwidth <- .check_bandwidth(bandwidth, 2L)
  scale <- match.arg(scale)
  treated_mediator <- .estimate_cf_hazard_impl(
    checked$data, exposure, exposure, bandwidth[1L],
    checked$time, checked$columns
  )
  reference_mediator <- .estimate_cf_hazard_impl(
    checked$data, exposure, reference, bandwidth[2L],
    checked$time, checked$columns
  )
  .effect_result(
    .hazard_to_scale(treated_mediator, reference_mediator, scale),
    checked$time, scale, "indirect",
    c(exposure = exposure, reference = reference)
  )
}
