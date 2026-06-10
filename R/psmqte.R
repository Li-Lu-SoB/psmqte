#' Propensity-score matching QTE and QTT
#'
#' Estimates quantile treatment effects (QTE) and quantile treatment effects on
#' the treated (QTT) by nearest-neighbor matching on a propensity score. The
#' default first step estimates the propensity score by logistic regression using
#' the right-hand side of `formula`.
#'
#' @param formula Formula of the form `outcome ~ covariates`.
#' @param treatment Name of a binary treatment variable in `data`.
#' @param data Data frame.
#' @param taus Quantile indexes in `(0, 1)`.
#' @param M Positive integer vector giving the number of matches.
#' @param ps_formula Optional formula for the propensity-score model. Defaults
#'   to `treatment ~ covariates`.
#' @param ps_score Optional externally supplied propensity score. If supplied,
#'   no propensity model is fit and first-step score adjustment is skipped.
#' @param estimands Character vector containing `"qte"` and/or `"qtt"`.
#' @param se Logical; compute feasible influence-function standard errors.
#' @param level Confidence level for pointwise intervals.
#' @param score_adjust Logical; apply estimated-propensity-score adjustment when
#'   the score is fit internally.
#' @param local_neighbors Number of same-arm nearest neighbors used to estimate
#'   conditional CDFs for standard errors.
#' @param multiplier Number of multiplier draws for uniform bands. Use `0` to
#'   skip uniform bands.
#' @param seed Optional random seed for multiplier bands.
#' @param keep_influence Logical; store influence-function matrices in the
#'   returned object.
#' @param eps Propensity score trimming away from 0 and 1.
#' @param fd_eps Finite-difference step for the QTT target-drift adjustment.
#'
#' @return An object of class `psmqte` with estimates, diagnostics, matches, and
#'   propensity-score information.
#' @export
psmqte <- function(formula, treatment, data, taus = seq(0.1, 0.9, by = 0.1),
                   M = 1L, ps_formula = NULL, ps_score = NULL,
                   estimands = c("qte", "qtt"), se = TRUE, level = 0.95,
                   score_adjust = TRUE, local_neighbors = NULL,
                   multiplier = 0L, seed = NULL, keep_influence = FALSE,
                   eps = 1e-8, fd_eps = 1e-4) {
  estimands <- match.arg(estimands, c("qte", "qtt"), several.ok = TRUE)
  taus <- as.numeric(taus)
  if (any(!is.finite(taus)) || any(taus <= 0 | taus >= 1)) {
    stop("taus must be finite values between 0 and 1.", call. = FALSE)
  }
  M <- sort(unique(as.integer(M)))
  if (length(M) == 0L || any(is.na(M)) || any(M < 1L)) {
    stop("M must contain positive integers.", call. = FALSE)
  }
  if (!isTRUE(se)) {
    multiplier <- 0L
  }
  if (level <= 0 || level >= 1) {
    stop("level must be between 0 and 1.", call. = FALSE)
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
  if (is.null(local_neighbors)) {
    local_neighbors <- max(10L, floor(sqrt(min(n0, n1))))
  }
  local_neighbors <- as.integer(local_neighbors)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.null(ps_score) && score_adjust) {
    warning("score_adjust is skipped when ps_score is supplied.", call. = FALSE)
    score_adjust <- FALSE
  }

  match_mat <- nearest_match_matrix(score, w, max(M))
  estimate_rows <- list()
  influence <- list()

  for (M_now in M) {
    K <- match_counts_from_matrix(match_mat, M_now, n)
    Kt <- match_counts_from_matrix(match_mat[w == 1, , drop = FALSE], M_now, n)
    w_qte_1 <- w * (1 + K / M_now)
    w_qte_0 <- (1 - w) * (1 + K / M_now)
    w_qtt_1 <- w
    w_qtt_0 <- (1 - w) * Kt / M_now

    q1 <- weighted_quantile(y, w_qte_1, taus)
    q0 <- weighted_quantile(y, w_qte_0, taus)
    q1T <- weighted_quantile(y, w_qtt_1, taus)
    q0T <- weighted_quantile(y, w_qtt_0, taus)

    if_qte <- matrix(NA_real_, nrow = n, ncol = length(taus))
    if_qtt <- matrix(NA_real_, nrow = n, ncol = length(taus))
    se_qte <- rep(NA_real_, length(taus))
    se_qtt <- rep(NA_real_, length(taus))

    if (se) {
      if ("qte" %in% estimands) {
        m1_q1 <- local_cdf_same_arm(y, w, score, 1, q1, local_neighbors)
        m0_q0 <- local_cdf_same_arm(y, w, score, 0, q0, local_neighbors)
      }
      if ("qtt" %in% estimands) {
        m0_q0T <- local_cdf_same_arm(y, w, score, 0, q0T, local_neighbors)
      }

      for (j in seq_along(taus)) {
        if ("qte" %in% estimands) {
          F1 <- weighted_cdf_at(y, w_qte_1, q1[j])
          F0 <- weighted_cdf_at(y, w_qte_0, q0[j])
          f1 <- weighted_density_at(y, w_qte_1, q1[j])
          f0 <- weighted_density_at(y, w_qte_0, q0[j])
          B1 <- as.numeric(y <= q1[j])
          B0 <- as.numeric(y <= q0[j])
          psi1 <- m1_q1[, j] - F1 + w * (1 + K / M_now) *
            (B1 - m1_q1[, j])
          psi0 <- m0_q0[, j] - F0 + (1 - w) * (1 + K / M_now) *
            (B0 - m0_q0[, j])
          raw <- -psi1 / f1 + psi0 / f0
          if (score_adjust) {
            raw <- score_adjust_if(raw, prep$x_mat, w, score)$IF_adj
          }
          if_qte[, j] <- raw
          se_qte[j] <- sqrt(max(stats::var(raw), variance_floor(n))) / sqrt(n)
        }

        if ("qtt" %in% estimands) {
          pi_hat <- n1 / n
          F1T <- weighted_cdf_at(y, w_qtt_1, q1T[j])
          F0T <- weighted_cdf_at(y, w_qtt_0, q0T[j])
          f1T <- weighted_density_at(y, w_qtt_1, q1T[j])
          f0T <- weighted_density_at(y, w_qtt_0, q0T[j])
          B1T <- as.numeric(y <= q1T[j])
          B0T <- as.numeric(y <= q0T[j])
          phi1T <- w / pi_hat * (B1T - F1T)
          phi0T <- score / pi_hat * (m0_q0T[, j] - F0T) +
            (1 - w) * (Kt / M_now) / pi_hat * (B0T - m0_q0T[, j])
          rawT <- -phi1T / f1T + phi0T / f0T
          if (score_adjust) {
            drift <- estimate_drift_qtt_fd(
              tau = taus[j], theta_hat = prep$theta_hat, x_mat = prep$x_mat,
              y = y, w = w, phat = score, eps = fd_eps
            )
            rawT <- score_adjust_if(rawT, prep$x_mat, w, score,
                                    drift = drift)$IF_adj
          }
          if_qtt[, j] <- rawT
          se_qtt[j] <- sqrt(max(stats::var(rawT), variance_floor(n))) / sqrt(n)
        }
      }
    }

    z <- stats::qnorm((1 + level) / 2)
    row <- data.frame(
      M = M_now,
      tau = taus,
      qte = if ("qte" %in% estimands) q1 - q0 else NA_real_,
      qte_se = se_qte,
      qte_ci_low = if ("qte" %in% estimands) q1 - q0 - z * se_qte else NA_real_,
      qte_ci_high = if ("qte" %in% estimands) q1 - q0 + z * se_qte else NA_real_,
      qte_q1 = if ("qte" %in% estimands) q1 else NA_real_,
      qte_q0 = if ("qte" %in% estimands) q0 else NA_real_,
      qtt = if ("qtt" %in% estimands) q1T - q0T else NA_real_,
      qtt_se = se_qtt,
      qtt_ci_low = if ("qtt" %in% estimands) q1T - q0T - z * se_qtt else NA_real_,
      qtt_ci_high = if ("qtt" %in% estimands) q1T - q0T + z * se_qtt else NA_real_,
      qtt_q1 = if ("qtt" %in% estimands) q1T else NA_real_,
      qtt_q0 = if ("qtt" %in% estimands) q0T else NA_real_,
      n = n,
      n_treated = n1,
      n_control = n0,
      stringsAsFactors = FALSE
    )

    if (multiplier > 0L && se) {
      row <- add_uniform_bands(row, if_qte, if_qtt, estimands, level,
                               multiplier, n)
    }

    estimate_rows[[as.character(M_now)]] <- row
    if (keep_influence) {
      influence[[as.character(M_now)]] <- list(qte = if_qte, qtt = if_qtt)
    }
  }

  estimates <- do.call(rbind, estimate_rows)
  rownames(estimates) <- NULL
  diagnostics <- list(
    overlap = overlap_diagnostics(score, w),
    reuse = reuse_diagnostics(y, w, match_mat, M)
  )

  out <- list(
    call = match.call(),
    estimates = estimates,
    diagnostics = diagnostics,
    matches = list(match_matrix = match_mat, M = M),
    propensity = list(score = score, fit = prep$ps_fit,
                      formula = prep$ps_formula),
    data_info = list(n = n, n_treated = n1, n_control = n0,
                     row_id = prep$row_id, outcome = prep$outcome,
                     treatment = treatment),
    settings = list(taus = taus, M = M, estimands = estimands, se = se,
                    level = level, score_adjust = score_adjust,
                    local_neighbors = local_neighbors,
                    multiplier = multiplier)
  )
  if (keep_influence) {
    out$influence <- influence
  }
  class(out) <- "psmqte"
  out
}

score_adjust_if <- function(if_raw, x_mat, w, phat, drift = NULL) {
  n <- length(w)
  score_mat <- x_mat * as.numeric(w - phat)
  info <- crossprod(x_mat, x_mat * as.numeric(phat * (1 - phat))) / n
  rho <- colMeans(score_mat * as.numeric(if_raw))
  if (is.null(drift)) {
    coef <- -solve_ridge(info, rho)
  } else {
    coef <- solve_ridge(info, drift - rho)
  }
  list(
    IF_adj = as.numeric(if_raw + score_mat %*% coef),
    rho = rho,
    drift = drift,
    coef = coef
  )
}

estimate_drift_qtt_fd <- function(tau, theta_hat, x_mat, y, w, phat,
                                  eps = 1e-4) {
  k <- length(theta_hat)
  deriv <- numeric(k)
  target <- function(theta) {
    ptheta <- clip01(as.numeric(logit_inv(x_mat %*% theta)))
    wt1 <- w * ptheta / phat
    wt0 <- (1 - w) * ptheta / (1 - phat)
    weighted_quantile(y, wt1, tau) - weighted_quantile(y, wt0, tau)
  }
  for (j in seq_len(k)) {
    step <- rep(0, k)
    step[j] <- eps
    deriv[j] <- (target(theta_hat + step) - target(theta_hat - step)) /
      (2 * eps)
  }
  deriv
}

add_uniform_bands <- function(row, if_qte, if_qtt, estimands, level, B, n) {
  add_one <- function(row, IF, prefix) {
    ok <- colSums(is.finite(IF)) == n
    if (!any(ok)) {
      return(row)
    }
    sds <- sqrt(pmax(apply(IF, 2, stats::var), variance_floor(n)))
    sds[sds <= 1e-12 | !is.finite(sds)] <- NA_real_
    z <- matrix(stats::rnorm(n * B), nrow = n, ncol = B)
    process <- crossprod(IF, z) / sqrt(n)
    stat <- apply(abs(process / sds), 2, max, na.rm = TRUE)
    crit <- as.numeric(stats::quantile(stat, level, names = FALSE, na.rm = TRUE))
    est <- row[[prefix]]
    se <- row[[paste0(prefix, "_se")]]
    row[[paste0(prefix, "_uniform_crit")]] <- crit
    row[[paste0(prefix, "_uniform_low")]] <- est - crit * se
    row[[paste0(prefix, "_uniform_high")]] <- est + crit * se
    row
  }
  if ("qte" %in% estimands) {
    row <- add_one(row, if_qte, "qte")
  }
  if ("qtt" %in% estimands) {
    row <- add_one(row, if_qtt, "qtt")
  }
  row
}
