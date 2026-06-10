#' Standardized mean-difference balance table
#'
#' Computes numeric covariate balance before or after applying separate treated
#' and control weights.
#'
#' @param data Data frame.
#' @param treatment Name of binary treatment variable.
#' @param covariates Character vector of numeric covariates.
#' @param weights_treated Optional weights for the treated mean.
#' @param weights_control Optional weights for the control mean.
#'
#' @return A data frame with weighted means, variances, and standardized mean
#'   differences.
#' @export
psm_balance <- function(data, treatment, covariates,
                        weights_treated = NULL, weights_control = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  data <- as.data.frame(data)
  if (!treatment %in% names(data)) {
    stop("treatment variable not found in data.", call. = FALSE)
  }
  w <- check_binary_treatment(data[[treatment]], treatment)
  if (is.null(weights_treated)) {
    weights_treated <- as.numeric(w == 1)
  }
  if (is.null(weights_control)) {
    weights_control <- as.numeric(w == 0)
  }

  rows <- lapply(covariates, function(v) {
    if (!v %in% names(data)) {
      stop("covariate not found in data: ", v, call. = FALSE)
    }
    x <- data[[v]]
    if (!is.numeric(x)) {
      stop("psm_balance currently supports numeric covariates only: ", v,
           call. = FALSE)
    }
    mt <- weighted_mean(x, weights_treated)
    mc <- weighted_mean(x, weights_control)
    vt <- weighted_var(x, weights_treated)
    vc <- weighted_var(x, weights_control)
    sd_pool <- sqrt((vt + vc) / 2)
    smd <- if (is.finite(sd_pool) && sd_pool > 0) (mt - mc) / sd_pool else NA_real_
    data.frame(
      variable = v,
      mean_treated = mt,
      mean_control = mc,
      var_treated = vt,
      var_control = vc,
      smd = smd,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

overlap_diagnostics <- function(score, w) {
  by_arm <- lapply(0:1, function(arm) {
    s <- score[w == arm]
    data.frame(
      treatment = arm,
      n = length(s),
      min = min(s),
      p05 = as.numeric(stats::quantile(s, 0.05, names = FALSE)),
      median = stats::median(s),
      p95 = as.numeric(stats::quantile(s, 0.95, names = FALSE)),
      max = max(s)
    )
  })
  by_arm <- do.call(rbind, by_arm)
  common <- data.frame(
    lower = max(by_arm$min),
    upper = min(by_arm$max),
    stringsAsFactors = FALSE
  )
  list(by_arm = by_arm, common_support = common)
}

reuse_diagnostics <- function(y, w, match_mat, M_values) {
  n <- length(w)
  out <- lapply(M_values, function(M) {
    K <- match_counts_from_matrix(match_mat, M, n)
    Kt <- match_counts_from_matrix(match_mat[w == 1, , drop = FALSE], M, n)
    w_qte_1 <- w * (1 + K / M)
    w_qte_0 <- (1 - w) * (1 + K / M)
    w_qtt_1 <- w
    w_qtt_0 <- (1 - w) * Kt / M
    data.frame(
      M = M,
      max_reuse_all = max(K),
      max_reuse_treated_target = max(Kt),
      ess_qte_treated = effective_n(w_qte_1),
      ess_qte_control = effective_n(w_qte_0),
      ess_qtt_treated = effective_n(w_qtt_1),
      ess_qtt_control = effective_n(w_qtt_0)
    )
  })
  do.call(rbind, out)
}
