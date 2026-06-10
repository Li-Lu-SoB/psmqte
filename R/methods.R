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
plot.psmqte <- function(x, estimand = c("qte", "qtt"), M = NULL,
                        ci = TRUE, uniform = FALSE,
                        legend = TRUE,
                        legend_position = c("bottom", "topright", "topleft",
                                            "bottomright", "none"),
                        xlab = "Quantile index", ylab = NULL, main = NULL,
                        ylim = NULL, grid = TRUE, zero_line = TRUE,
                        family = "serif", ...) {
  estimand <- match.arg(estimand)
  legend_position <- match.arg(legend_position)
  dat <- x$estimates
  if (!is.null(M)) {
    dat <- dat[dat$M %in% M, , drop = FALSE]
  }
  if (nrow(dat) == 0L) {
    stop("No estimates to plot for the requested M.", call. = FALSE)
  }

  m_values <- sort(unique(dat$M))
  tau_values <- sort(unique(dat$tau))
  if (is.null(ylab)) {
    ylab <- paste0(toupper(estimand), " estimate")
  }

  interval_cols <- journal_interval_columns(dat, estimand, uniform)
  if (ci && !interval_cols$available) {
    ci <- FALSE
  }

  y_range <- dat[[estimand]]
  if (ci) {
    y_range <- c(y_range, dat[[interval_cols$low]], dat[[interval_cols$high]])
  }
  y_range <- y_range[is.finite(y_range)]
  if (length(y_range) == 0L) {
    stop("No finite estimates are available to plot.", call. = FALSE)
  }
  if (is.null(ylim)) {
    ylim <- range(y_range)
    pad <- diff(ylim) * 0.08
    if (!is.finite(pad) || pad <= 0) {
      pad <- max(0.1, abs(ylim[1]) * 0.08)
    }
    ylim <- ylim + c(-pad, pad)
  }

  line_types <- journal_line_types(length(m_values))
  point_shapes <- journal_point_shapes(length(m_values))
  line_cols <- rep("black", length(m_values))
  ci_cols <- rep("gray55", length(m_values))

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    family = family,
    mar = c(if (legend && legend_position == "bottom") 6.2 else 5.1,
            5.0, 2.0, 1.2),
    mgp = c(2.8, 0.8, 0),
    las = 1,
    tcl = -0.25,
    bty = "l"
  )

  graphics::plot(
    range(tau_values), ylim, type = "n", axes = FALSE,
    xlab = xlab, ylab = ylab, main = main, ...
  )

  if (grid) {
    y_ticks <- graphics::axTicks(2)
    graphics::abline(h = y_ticks, col = "gray90", lwd = 0.45)
    graphics::abline(v = tau_values, col = "gray94", lwd = 0.35)
  }
  if (zero_line && ylim[1] <= 0 && ylim[2] >= 0) {
    graphics::abline(h = 0, lty = 3, col = "gray35", lwd = 0.7)
  }

  graphics::axis(1, at = tau_values, labels = format(tau_values, trim = TRUE),
                 cex.axis = 0.85)
  graphics::axis(2, cex.axis = 0.85)
  graphics::box()

  for (j in seq_along(m_values)) {
    one <- dat[dat$M == m_values[j], , drop = FALSE]
    one <- one[order(one$tau), , drop = FALSE]
    ok <- is.finite(one$tau) & is.finite(one[[estimand]])
    if (!any(ok)) {
      next
    }

    if (ci) {
      low <- one[[interval_cols$low]]
      high <- one[[interval_cols$high]]
      ci_ok <- ok & is.finite(low) & is.finite(high)
      if (any(ci_ok)) {
        journal_error_bars(
          x = one$tau[ci_ok], low = low[ci_ok], high = high[ci_ok],
          col = ci_cols[j], lty = line_types[j]
        )
      }
    }

    graphics::lines(
      one$tau[ok], one[[estimand]][ok],
      type = "b", lty = line_types[j], lwd = 0.8,
      pch = point_shapes[j], cex = 0.75,
      col = line_cols[j]
    )
  }

  if (legend && legend_position != "none") {
    legend_args <- list(
      legend = paste0("M = ", m_values),
      lty = line_types,
      pch = point_shapes,
      lwd = 0.8,
      col = line_cols,
      bty = "n",
      cex = 0.85
    )
    if (legend_position == "bottom") {
      old_xpd <- graphics::par("xpd")
      graphics::par(xpd = NA)
      do.call(graphics::legend, c(
        list(x = "bottom", inset = c(0, -0.28), horiz = TRUE),
        legend_args
      ))
      graphics::par(xpd = old_xpd)
    } else {
      do.call(graphics::legend, c(list(x = legend_position), legend_args))
    }
  }

  invisible(x)
}

journal_line_types <- function(n) {
  rep(c("solid", "dashed", "dotdash", "longdash", "twodash", "dotted"),
      length.out = n)
}

journal_point_shapes <- function(n) {
  rep(c(16, 17, 15, 1, 2, 0, 5, 6), length.out = n)
}

journal_interval_columns <- function(dat, estimand, uniform) {
  if (isTRUE(uniform)) {
    low <- paste0(estimand, "_uniform_low")
    high <- paste0(estimand, "_uniform_high")
  } else {
    low <- paste0(estimand, "_ci_low")
    high <- paste0(estimand, "_ci_high")
  }
  list(
    low = low,
    high = high,
    available = all(c(low, high) %in% names(dat)) &&
      any(is.finite(dat[[low]])) &&
      any(is.finite(dat[[high]]))
  )
}

journal_error_bars <- function(x, low, high, col, lty) {
  graphics::segments(x, low, x, high, col = col, lwd = 0.55, lty = lty)
  tick <- diff(range(x))
  if (!is.finite(tick) || tick <= 0) {
    tick <- 0.01
  } else {
    tick <- tick * 0.012
  }
  graphics::segments(x - tick, low, x + tick, low, col = col, lwd = 0.55)
  graphics::segments(x - tick, high, x + tick, high, col = col, lwd = 0.55)
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
