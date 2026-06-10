\name{psmate}
\alias{psmate}
\title{Propensity-Score Matching ATE and ATT}
\usage{
psmate(
  formula,
  treatment,
  data,
  M = 1L,
  ps_formula = NULL,
  ps_score = NULL,
  estimands = c("ate", "att"),
  se = TRUE,
  level = 0.95,
  eps = 1e-08
)
}
\arguments{
\item{formula}{Formula of the form \code{outcome ~ covariates}.}
\item{treatment}{Name of a binary treatment variable in \code{data}.}
\item{data}{Data frame.}
\item{M}{Positive integer vector giving the number of matches.}
\item{ps_formula}{Optional formula for the propensity-score model. Defaults to
\code{treatment ~ covariates}.}
\item{ps_score}{Optional externally supplied propensity score.}
\item{estimands}{Character vector containing \code{"ate"} and/or \code{"att"}.}
\item{se}{Logical; compute large-sample standard errors.}
\item{level}{Confidence level for pointwise intervals.}
\item{eps}{Propensity score trimming away from 0 and 1.}
}
\value{
An object of class \code{psmate}. Its \code{estimates} component contains ATE
and ATT estimates by \code{M}.
}
\description{
Estimates average treatment effects (ATE) and average treatment effects on the
treated (ATT) by nearest-neighbor matching on a propensity score.
}
\examples{
dat <- simulate_psmqte_data(300, seed = 2)
fit <- psmate(Y ~ X1 + X2, treatment = "D", data = dat, M = c(1, 2))
summary(fit)
}
