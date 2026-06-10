\name{psm_balance}
\alias{psm_balance}
\title{Standardized Mean-Difference Balance Table}
\usage{
psm_balance(data, treatment, covariates, weights_treated = NULL, weights_control = NULL)
}
\arguments{
\item{data}{Data frame.}
\item{treatment}{Name of binary treatment variable.}
\item{covariates}{Character vector of numeric covariates.}
\item{weights_treated}{Optional weights for the treated mean.}
\item{weights_control}{Optional weights for the control mean.}
}
\value{
A data frame with weighted means, variances, and standardized mean differences.
}
\description{
Computes numeric covariate balance before or after applying separate treated and
control weights.
}
\examples{
dat <- simulate_psmqte_data(200, seed = 4)
psm_balance(dat, "D", c("X1", "X2"))
}
