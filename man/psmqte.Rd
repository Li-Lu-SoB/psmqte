\name{psmqte}
\alias{psmqte}
\title{Propensity-Score Matching QTE and QTT}
\usage{
psmqte(
  formula,
  treatment,
  data,
  taus = seq(0.1, 0.9, by = 0.1),
  M = 1L,
  ps_formula = NULL,
  ps_score = NULL,
  estimands = c("qte", "qtt"),
  se = TRUE,
  level = 0.95,
  score_adjust = TRUE,
  local_neighbors = NULL,
  multiplier = 0L,
  seed = NULL,
  keep_influence = FALSE,
  eps = 1e-08,
  fd_eps = 1e-04
)
}
\arguments{
\item{formula}{Formula of the form \code{outcome ~ covariates}.}
\item{treatment}{Name of a binary treatment variable in \code{data}.}
\item{data}{Data frame.}
\item{taus}{Quantile indexes in \code{(0, 1)}.}
\item{M}{Positive integer vector giving the number of matches.}
\item{ps_formula}{Optional formula for the propensity-score model. Defaults to
\code{treatment ~ covariates}.}
\item{ps_score}{Optional externally supplied propensity score.}
\item{estimands}{Character vector containing \code{"qte"} and/or \code{"qtt"}.}
\item{se}{Logical; compute feasible influence-function standard errors.}
\item{level}{Confidence level for pointwise intervals.}
\item{score_adjust}{Logical; apply estimated-propensity-score adjustment when
the score is fit internally.}
\item{local_neighbors}{Number of same-arm nearest neighbors used to estimate
conditional CDFs for standard errors.}
\item{multiplier}{Number of multiplier draws for uniform bands. Use \code{0}
to skip uniform bands.}
\item{seed}{Optional random seed for multiplier bands.}
\item{keep_influence}{Logical; store influence-function matrices in the
returned object.}
\item{eps}{Propensity score trimming away from 0 and 1.}
\item{fd_eps}{Finite-difference step for the QTT target-drift adjustment.}
}
\value{
An object of class \code{psmqte}. Its \code{estimates} component is a data frame
with one row per \code{M} and quantile index. The object also includes
propensity-score information, matching diagnostics, overlap diagnostics, and
reuse diagnostics.
}
\description{
Estimates quantile treatment effects (QTE) and quantile treatment effects on the
treated (QTT) by nearest-neighbor matching on a propensity score.
}
\details{
The default first step estimates a logistic propensity-score model using the
right-hand side of \code{formula}. Users can instead supply \code{ps_formula} or
an externally computed \code{ps_score}. Matching is with replacement and uses
deterministic tie breaking by row order.
}
\examples{
dat <- simulate_psmqte_data(300, seed = 1)
fit <- psmqte(
  Y ~ X1 + X2,
  treatment = "D",
  data = dat,
  taus = c(0.25, 0.5, 0.75),
  M = c(1, 2),
  se = FALSE
)
summary(fit)
}
