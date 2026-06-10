#' Build propensity-score nearest-neighbor matches
#'
#' Finds the `M` nearest opposite-treatment observations for every unit using a
#' scalar score, usually an estimated propensity score. Matching is with
#' replacement and deterministic tie breaking by row order.
#'
#' @param score Numeric matching score.
#' @param treatment Binary treatment vector coded as 0/1, logical, or
#'   control/treated.
#' @param M Positive integer number of matches.
#'
#' @return A list with a match matrix, reuse counts, score, treatment, and arm
#'   sizes.
#' @export
psm_match <- function(score, treatment, M = 1L) {
  score <- as.numeric(score)
  w <- check_binary_treatment(treatment)
  if (length(score) != length(w)) {
    stop("score and treatment must have the same length.", call. = FALSE)
  }
  keep <- is.finite(score) & !is.na(w)
  if (!all(keep)) {
    score <- score[keep]
    w <- w[keep]
  }
  M <- as.integer(M)
  if (length(M) != 1L || is.na(M) || M < 1L) {
    stop("M must be a single positive integer.", call. = FALSE)
  }
  n0 <- sum(w == 0)
  n1 <- sum(w == 1)
  if (M > min(n0, n1)) {
    stop("M cannot exceed the smaller treatment-arm size.", call. = FALSE)
  }

  match_mat <- nearest_match_matrix(score, w, M)
  list(
    match_matrix = match_mat,
    K_all = match_counts_from_matrix(match_mat, M, length(w)),
    K_treated = match_counts_from_matrix(match_mat[w == 1, , drop = FALSE],
                                         M, length(w)),
    score = score,
    treatment = w,
    M = M,
    n = length(w),
    n1 = n1,
    n0 = n0
  )
}

nearest_match_matrix <- function(score, w, M_max) {
  n <- length(score)
  out <- matrix(NA_integer_, nrow = n, ncol = M_max)

  for (target_arm in 0:1) {
    target_idx <- which(w == target_arm)
    donor_idx <- which(w == 1 - target_arm)
    donor_order <- order(score[donor_idx], donor_idx)
    donor_idx <- donor_idx[donor_order]
    donor_score <- score[donor_idx]
    n_donor <- length(donor_idx)
    if (n_donor < M_max) {
      stop("Not enough opposite-treatment donors for M = ", M_max,
           call. = FALSE)
    }

    for (i in target_idx) {
      pos <- findInterval(score[i], donor_score)
      lo <- max(1L, pos - M_max)
      hi <- min(n_donor, pos + M_max + 1L)
      cand <- donor_idx[lo:hi]
      dist <- abs(score[cand] - score[i])
      out[i, ] <- cand[order(dist, cand)][seq_len(M_max)]
    }
  }

  out
}

match_counts_from_matrix <- function(match_mat, M, n) {
  if (nrow(match_mat) == 0L) {
    return(integer(n))
  }
  tabulate(as.vector(match_mat[, seq_len(M), drop = FALSE]), nbins = n)
}

#' Diverging-M grid used in matching simulations
#'
#' Returns `M = 2^j` for `0 <= j <= floor(log2(n) / 2)`, the grid used in
#' diverging-number-of-matches simulation designs.
#'
#' @param n Sample size.
#' @param max_M Optional upper bound.
#'
#' @return Integer vector of matching counts.
#' @export
m_grid_diverging <- function(n, max_M = NULL) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("n must be a sample size of at least 2.", call. = FALSE)
  }
  grid <- as.integer(2^(0:floor(log2(n) / 2)))
  if (!is.null(max_M)) {
    grid <- grid[grid <= max_M]
  }
  unique(grid)
}
