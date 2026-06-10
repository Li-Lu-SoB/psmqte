\name{psmqte_example_data}
\alias{psmqte_example_data}
\title{Load an Example Data Set for psmqte}
\usage{
psmqte_example_data(source = c("demo", "simulated", "nhefs"), n = 1000L, seed = 1L)
}
\arguments{
\item{source}{One of \code{"demo"}, \code{"simulated"}, or \code{"nhefs"}.}
\item{n}{Sample size for the simulated example.}
\item{seed}{Optional random seed for the simulated example.}
}
\value{
A data frame with outcome \code{Y}, treatment \code{D}, and covariates.
}
\description{
Loads an example data set suitable for demonstrating propensity-score matching
QTE/QTT estimation. The default \code{"demo"} source is a fixed synthetic data
set shipped inside this package. The optional \code{"nhefs"} source uses the
public NHEFS example data from the \pkg{causaldata} package without
redistributing those data inside \pkg{psmqte}.
}
\examples{
dat <- psmqte_example_data()
fit <- psmqte(Y ~ X1 + X2 + X3 + age + size,
              treatment = "D", data = dat, se = FALSE)
summary(fit)
}
