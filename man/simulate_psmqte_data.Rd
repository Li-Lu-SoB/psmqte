\name{simulate_psmqte_data}
\alias{simulate_psmqte_data}
\title{Simulate Data for psmqte Examples}
\usage{
simulate_psmqte_data(n = 1000L, seed = NULL)
}
\arguments{
\item{n}{Sample size.}
\item{seed}{Optional random seed.}
}
\value{
A data frame with observed outcome \code{Y}, treatment \code{D}, covariates
\code{X1} and \code{X2}, potential outcomes, and the true propensity score.
}
\description{
Creates a simple observational treatment-effect design with heterogeneous
quantile effects and selection on observables.
}
\examples{
dat <- simulate_psmqte_data(100, seed = 5)
head(dat)
}
