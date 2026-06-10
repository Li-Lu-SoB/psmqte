#' @export
print.psmqte <- function(x, ...) {
  cat("psmqte: propensity-score matching QTE/QTT\n")
  cat("Observations:", x$data_info$n,
      " treated:", x$data_info$n_treated,
      " controls:", x$data_info$n_control, "\n")
  cat("M:", paste(x$settings$M, collapse = ", "),
      " taus:", paste(x$settings$taus, collapse = ", "), "\n")
  print(utils::head(x$estimates), row.names = FALSE)
  invisible(x)
}

#' @export
summary.psmqte <- function(object, ...) {
  object$estimates
}

#' @export
coef.psmqte <- function(object, estimand = c("qte", "qtt"), ...) {
  estimand <- match.arg(estimand)
  out <- object$estimates[, c("M", "tau", estimand), drop = FALSE]
  names(out)[3L] <- "estimate"
  out
}

#' @export
plot.psmqte <- function(x, estimand = c("qte", "qtt"), M = NULL, ...) {
  estimand <- match.arg(estimand)
  dat <- x$estimates
  if (!is.null(M)) {
    dat <- dat[dat$M %in% M, , drop = FALSE]
  }
  if (nrow(dat) == 0L) {
    stop("No estimates to plot for the requested M.", call. = FALSE)
  }
  m_values <- sort(unique(dat$M))
  y_mat <- sapply(m_values, function(mm) {
    dat[dat$M == mm, estimand]
  })
  if (is.null(dim(y_mat))) {
    y_mat <- matrix(y_mat, ncol = 1L)
  }
  graphics::matplot(
    sort(unique(dat$tau)), y_mat, type = "b", pch = seq_along(m_values),
    lty = seq_along(m_values), xlab = "Quantile index",
    ylab = toupper(estimand), ...
  )
  graphics::abline(h = 0, lty = 3, col = "gray50")
  graphics::legend("topleft", legend = paste0("M=", m_values),
                   lty = seq_along(m_values), pch = seq_along(m_values),
                   bty = "n")
  invisible(x)
}

#' @export
print.psmate <- function(x, ...) {
  cat("psmate: propensity-score matching ATE/ATT\n")
  cat("Observations:", x$data_info$n,
      " treated:", x$data_info$n_treated,
      " controls:", x$data_info$n_control, "\n")
  cat("M:", paste(x$settings$M, collapse = ", "), "\n")
  print(x$estimates, row.names = FALSE)
  invisible(x)
}

#' @export
summary.psmate <- function(object, ...) {
  object$estimates
}

#' @export
coef.psmate <- function(object, estimand = c("ate", "att"), ...) {
  estimand <- match.arg(estimand)
  out <- object$estimates[, c("M", estimand), drop = FALSE]
  names(out)[2L] <- "estimate"
  out
}
