#' Small simulated illness-death data
#'
#' Generates a compact data set for help examples and documentation checks.
#' This helper is intentionally not exported because it is not part of the
#' estimation API.
#'
#' @param n Number of observations.
#' @param seed Random seed.
#' @return A data frame with the default columns expected by the package.
#' @keywords internal
nm_example_data <- function(n = 80L, seed = 42) {
  set.seed(seed)
  Z <- stats::rbinom(n, 1, 0.5)
  latent_t1 <- stats::rexp(n, exp(0.25 * Z))
  latent_t2_direct <- stats::rexp(n, exp(0.15 * Z))
  D1 <- as.integer(latent_t1 < latent_t2_direct)
  T1_hat <- pmin(latent_t1, latent_t2_direct)
  post_t1 <- stats::rexp(n, exp(0.2 * Z + 0.2 * D1))
  T2_hat <- ifelse(D1 == 1, latent_t1 + post_t1, latent_t2_direct)
  data.frame(
    T1_hat = T1_hat,
    D1 = D1,
    T2_hat = T2_hat,
    D2 = 1L,
    Z = Z
  )
}
