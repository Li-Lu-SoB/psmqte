#' Propensity-score matching ATE and ATT
#'
#' Estimates average treatment effects (ATE) and average treatment effects on
#' the treated (ATT) by nearest-neighbor matching on a propensity score. This is
#' included as a companion to `psmqte()` for users who want a full treatment
#' effect workflow.
#'
#' @inheritParams psmqte
#' @param estimands Character vector containing `"ate"` and/or `"att"`.
#'
#' @return An object of class `psmate`.
#' @export
psmate <- function(formula, treatment, data, M = 1L, ps_formula = NULL,
                   ps_score = NULL, estimands = c("ate", "att"), se = TRUE,
                   level = 0.95, eps = 1e-8) {
  estimands <- match.arg(estimands, c("ate", "att"), several.ok = TRUE)
  M <- sort(unique(as.integer(M)))
  if (length(M) == 0L || any(is.na(M)) || any(M < 1L)) {
    stop("M must contain positive integers.", call. = FALSE)
  }
  prep <- prepare_psm_data(formula, treatment, data, ps_formula, ps_score, eps)
  y <- prep$y
  w <- prep$w
  score <- prep$score
  n <- prep$n
  n1 <- prep$n1
  n0 <- prep$n0
  if (max(M) > min(n0, n1)) {
    stop("max(M) cannot exceed the smaller treatment-arm size.", call. = FALSE)
  }

  match_mat <- nearest_match_matrix(score, w, max(M))
  z <- stats::qnorm((1 + level) / 2)
  rows <- lapply(M, function(M_now) {
    K <- match_counts_from_matrix(match_mat, M_now, n)
    Kt <- match_counts_from_matrix(match_mat[w == 1, , drop = FALSE], M_now, n)

    wt1_ate <- w * (1 + K / M_now)
    wt0_ate <- (1 - w) * (1 + K / M_now)
    mu1_ate <- sum(wt1_ate * y) / n
    mu0_ate <- sum(wt0_ate * y) / n
    ate <- mu1_ate - mu0_ate

    mu1_att <- sum(w * y) / n1
    mu0_att <- sum((1 - w) * (Kt / M_now) * y) / n1
    att <- mu1_att - mu0_att

    ate_se <- NA_real_
    att_se <- NA_real_
    if (se) {
      phi_ate <- wt1_ate * (y - mu1_ate) - wt0_ate * (y - mu0_ate)
      ate_se <- sqrt(stats::var(phi_ate) / n)
      pi_hat <- n1 / n
      phi_att <- (w * (y - mu1_att) -
        (1 - w) * (Kt / M_now) * (y - mu0_att)) / pi_hat
      att_se <- sqrt(stats::var(phi_att) / n)
    }

    data.frame(
      M = M_now,
      ate = if ("ate" %in% estimands) ate else NA_real_,
      ate_se = ate_se,
      ate_ci_low = if ("ate" %in% estimands) ate - z * ate_se else NA_real_,
      ate_ci_high = if ("ate" %in% estimands) ate + z * ate_se else NA_real_,
      ate_y1 = if ("ate" %in% estimands) mu1_ate else NA_real_,
      ate_y0 = if ("ate" %in% estimands) mu0_ate else NA_real_,
      att = if ("att" %in% estimands) att else NA_real_,
      att_se = att_se,
      att_ci_low = if ("att" %in% estimands) att - z * att_se else NA_real_,
      att_ci_high = if ("att" %in% estimands) att + z * att_se else NA_real_,
      att_y1 = if ("att" %in% estimands) mu1_att else NA_real_,
      att_y0 = if ("att" %in% estimands) mu0_att else NA_real_,
      n = n,
      n_treated = n1,
      n_control = n0,
      stringsAsFactors = FALSE
    )
  })

  estimates <- do.call(rbind, rows)
  rownames(estimates) <- NULL
  out <- list(
    call = match.call(),
    estimates = estimates,
    diagnostics = list(
      overlap = overlap_diagnostics(score, w),
      reuse = reuse_diagnostics(y, w, match_mat, M)
    ),
    matches = list(match_matrix = match_mat, M = M),
    propensity = list(score = score, fit = prep$ps_fit,
                      formula = prep$ps_formula),
    data_info = list(n = n, n_treated = n1, n_control = n0,
                     row_id = prep$row_id, outcome = prep$outcome,
                     treatment = treatment),
    settings = list(M = M, estimands = estimands, se = se, level = level)
  )
  class(out) <- "psmate"
  out
}
