clip01 <- function(p, eps = 1e-8) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

logit_inv <- function(x) {
  1 / (1 + exp(-x))
}

check_binary_treatment <- function(x, name = "treatment") {
  if (is.logical(x)) {
    return(as.numeric(x))
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    y <- rep(NA_real_, length(x))
    y[x %in% c("0", "control", "Control", "FALSE", "false")] <- 0
    y[x %in% c("1", "treated", "Treated", "TRUE", "true")] <- 1
    if (any(is.na(y) & !is.na(x))) {
      stop(name, " must be coded as 0/1, TRUE/FALSE, or control/treated.",
           call. = FALSE)
    }
    return(y)
  }

  y <- as.numeric(x)
  vals <- sort(unique(y[!is.na(y)]))
  if (!all(vals %in% c(0, 1))) {
    stop(name, " must be binary and coded as 0/1.", call. = FALSE)
  }
  y
}

formula_response_name <- function(formula) {
  vars <- all.vars(formula[[2]])
  if (length(vars) != 1L) {
    stop("The left side of formula must contain exactly one outcome variable.",
         call. = FALSE)
  }
  vars
}

rhs_formula <- function(formula) {
  stats::reformulate(attr(stats::terms(formula), "term.labels"))
}

default_ps_formula <- function(formula, treatment) {
  terms <- attr(stats::terms(formula), "term.labels")
  if (length(terms) == 0L) {
    stop("The formula must include at least one covariate.", call. = FALSE)
  }
  stats::reformulate(terms, response = treatment)
}

prepare_psm_data <- function(formula, treatment, data, ps_formula = NULL,
                             ps_score = NULL, eps = 1e-8) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("formula must look like outcome ~ covariates.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("data must be a data.frame or compatible object.", call. = FALSE)
  }
  data <- as.data.frame(data)
  if (!is.character(treatment) || length(treatment) != 1L) {
    stop("treatment must be the name of a binary treatment variable.",
         call. = FALSE)
  }
  if (!treatment %in% names(data)) {
    stop("treatment variable not found in data: ", treatment, call. = FALSE)
  }
  if (is.null(ps_formula)) {
    ps_formula <- default_ps_formula(formula, treatment)
  }
  if (!inherits(ps_formula, "formula") || length(ps_formula) != 3L) {
    stop("ps_formula must look like treatment ~ covariates.", call. = FALSE)
  }

  y_frame <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(y_frame)

  data_work <- data
  data_work[[treatment]] <- check_binary_treatment(data[[treatment]], treatment)
  ps_frame <- stats::model.frame(ps_formula, data = data_work,
                                 na.action = stats::na.pass)
  w <- stats::model.response(ps_frame)

  keep <- stats::complete.cases(y_frame) & stats::complete.cases(ps_frame)
  if (!is.null(ps_score)) {
    if (length(ps_score) != nrow(data)) {
      stop("ps_score must have one value per row of data.", call. = FALSE)
    }
    keep <- keep & is.finite(ps_score)
  }

  if (!any(keep)) {
    stop("No complete observations remain after removing missing values.",
         call. = FALSE)
  }

  data_keep <- data_work[keep, , drop = FALSE]
  y <- as.numeric(y[keep])
  w <- as.numeric(w[keep])
  row_id <- which(keep)

  if (!all(w %in% c(0, 1))) {
    stop("treatment must be coded as 0/1 after missing values are removed.",
         call. = FALSE)
  }
  n1 <- sum(w == 1)
  n0 <- sum(w == 0)
  if (n1 == 0L || n0 == 0L) {
    stop("Both treatment arms must contain at least one observation.",
         call. = FALSE)
  }

  ps_fit <- NULL
  x_mat <- NULL
  theta_hat <- NULL
  if (is.null(ps_score)) {
    ps_fit <- stats::glm(ps_formula, family = stats::binomial(),
                         data = data_keep)
    score <- clip01(stats::fitted(ps_fit), eps = eps)
    x_mat <- stats::model.matrix(ps_fit)
    theta_hat <- as.numeric(stats::coef(ps_fit))
  } else {
    score <- clip01(ps_score[keep], eps = eps)
  }

  list(
    y = y,
    w = w,
    score = as.numeric(score),
    data = data_keep,
    row_id = row_id,
    formula = formula,
    treatment = treatment,
    ps_formula = ps_formula,
    ps_fit = ps_fit,
    x_mat = x_mat,
    theta_hat = theta_hat,
    outcome = formula_response_name(formula),
    n = length(y),
    n1 = n1,
    n0 = n0
  )
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]
  w <- w[ok]
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    return(NA_real_)
  }
  sum(w * x) / sw
}

weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]
  w <- w[ok]
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    return(NA_real_)
  }
  mu <- sum(w * x) / sw
  sum(w * (x - mu)^2) / sw
}

weighted_quantile <- function(y, w, probs) {
  if (length(y) != length(w)) {
    stop("y and w must have the same length.", call. = FALSE)
  }
  if (any(w < 0, na.rm = TRUE)) {
    stop("Weights must be nonnegative.", call. = FALSE)
  }
  ok <- is.finite(y) & is.finite(w) & w >= 0
  y <- y[ok]
  w <- w[ok]
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(y)
  y <- y[ord]
  w <- w[ord] / sw
  cw <- cumsum(w)
  vapply(probs, function(p) y[which(cw >= p)[1L]], numeric(1))
}

weighted_cdf_at <- function(y, w, q) {
  ok <- is.finite(y) & is.finite(w) & w >= 0
  y <- y[ok]
  w <- w[ok]
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    return(NA_real_)
  }
  sum(w * (y <= q)) / sw
}

weighted_density_at <- function(y, w, q, bw = NULL, floor = 1e-6) {
  ok <- is.finite(y) & is.finite(w) & w > 0
  y <- y[ok]
  w <- w[ok]
  if (length(y) < 5L || sum(w) <= 0) {
    return(floor)
  }
  w <- w / sum(w)
  if (is.null(bw)) {
    sdw <- sqrt(max(weighted_var(y, w), .Machine$double.eps))
    neff <- 1 / sum(w^2)
    bw <- 1.06 * sdw * neff^(-1 / 5)
    if (!is.finite(bw) || bw <= floor) {
      bw <- stats::sd(y) * length(y)^(-1 / 5)
    }
    if (!is.finite(bw) || bw <= floor) {
      bw <- 0.01
    }
  }
  max(sum(w * stats::dnorm((q - y) / bw)) / bw, floor)
}

effective_n <- function(w) {
  w <- as.numeric(w)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    return(NA_real_)
  }
  sw^2 / sum(w^2)
}

local_cdf_same_arm <- function(y, w, score, arm, q_values, L) {
  n <- length(y)
  idx <- which(w == arm)
  if (length(idx) == 0L) {
    stop("No observations in treatment arm ", arm, ".", call. = FALSE)
  }
  score_arm <- score[idx]
  ord <- order(score_arm, idx)
  idx_ord <- idx[ord]
  score_ord <- score[idx_ord]
  n_arm <- length(idx_ord)
  L_eff <- min(as.integer(L), max(1L, n_arm - 1L))

  out <- matrix(NA_real_, nrow = n, ncol = length(q_values))
  for (i in seq_len(n)) {
    pos <- findInterval(score[i], score_ord)
    lo <- max(1L, pos - L_eff - 2L)
    hi <- min(n_arm, pos + L_eff + 3L)
    cand <- idx_ord[lo:hi]
    if (w[i] == arm && length(cand) > 1L) {
      cand <- cand[cand != i]
    }
    if (length(cand) < L_eff) {
      cand <- idx
      if (w[i] == arm && length(cand) > 1L) {
        cand <- cand[cand != i]
      }
    }
    dist <- abs(score[cand] - score[i])
    chosen <- cand[order(dist, cand)][seq_len(min(L_eff, length(cand)))]
    out[i, ] <- vapply(q_values, function(q) mean(y[chosen] <= q), numeric(1))
  }
  out
}

solve_ridge <- function(a, b, ridge = 1e-8) {
  a <- as.matrix(a)
  diag(a) <- diag(a) + ridge
  out <- tryCatch(
    solve(a, b),
    error = function(e) qr.solve(a, b)
  )
  as.numeric(out)
}

variance_floor <- function(n) {
  n^(-2)
}
