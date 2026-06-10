\name{plot.psmqte}
\alias{plot.psmqte}
\title{Journal-Style Plot for psmqte Objects}
\usage{
\method{plot}{psmqte}(
  x,
  estimand = c("qte", "qtt"),
  M = NULL,
  ci = TRUE,
  uniform = FALSE,
  legend = TRUE,
  legend_position = c("bottom", "topright", "topleft", "bottomright", "none"),
  xlab = "Quantile index",
  ylab = NULL,
  main = NULL,
  ylim = NULL,
  grid = TRUE,
  zero_line = TRUE,
  family = "serif",
  ...
)
}
\arguments{
\item{x}{An object returned by \code{\link{psmqte}}.}
\item{estimand}{Either \code{"qte"} or \code{"qtt"}.}
\item{M}{Optional subset of matching counts to plot.}
\item{ci}{Logical; draw confidence intervals when they are available.}
\item{uniform}{Logical; use uniform confidence bands when available. If
\code{FALSE}, pointwise intervals are used.}
\item{legend}{Logical; draw a legend.}
\item{legend_position}{Legend position. The default places the legend beneath
the plotting panel.}
\item{xlab, ylab, main}{Axis labels and optional title.}
\item{ylim}{Optional y-axis limits.}
\item{grid}{Logical; draw a light journal-style grid.}
\item{zero_line}{Logical; draw a thin zero reference line when zero is inside
the y-axis range.}
\item{family}{Base graphics font family. The default \code{"serif"} matches the
package's journal-style figures.}
\item{...}{Additional arguments passed to \code{\link[graphics]{plot}}.}
}
\value{
Invisibly returns \code{x}.
}
\description{
Plots QTE or QTT estimates by quantile index using a monochrome, print-friendly
style with distinct line types and point symbols for each value of \code{M}.
}
\examples{
dat <- psmqte_example_data()
fit <- psmqte(
  Y ~ X1 + X2 + X3 + age + size,
  treatment = "D",
  data = dat,
  taus = seq(0.1, 0.9, by = 0.2),
  M = c(1, 2),
  se = FALSE
)
plot(fit, estimand = "qte")
plot(fit, estimand = "qtt", legend_position = "topright")
}
