#' Simulate data for psmqte examples
#'
#' Creates a simple observational treatment-effect design with heterogeneous
#' quantile effects and selection on observables.
#'
#' @param n Sample size.
#' @param seed Optional random seed.
#'
#' @return A data frame with observed outcome `Y`, treatment `D`, covariates
#'   `X1` and `X2`, potential outcomes, and the true propensity score.
#' @export
simulate_psmqte_data <- function(n = 1000L, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 10L) {
    stop("n must be at least 10.", call. = FALSE)
  }
  X1 <- stats::rnorm(n)
  X2 <- stats::runif(n, -1, 1)
  p_true <- logit_inv(-0.25 + 0.8 * X1 - 0.6 * X2)
  D <- stats::rbinom(n, 1, p_true)
  e0 <- stats::rnorm(n)
  e1 <- stats::rnorm(n, sd = 1 + 0.25 * abs(X1))
  Y0 <- 0.5 * X1 - 0.3 * X2 + e0
  Y1 <- 0.7 + 0.5 * X1 + 0.4 * X2 + 0.8 * abs(X1) + e1
  Y <- ifelse(D == 1, Y1, Y0)
  data.frame(
    Y = Y,
    D = D,
    X1 = X1,
    X2 = X2,
    Y0 = Y0,
    Y1 = Y1,
    p_true = p_true
  )
}
