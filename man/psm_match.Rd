\name{psm_match}
\alias{psm_match}
\title{Build Propensity-Score Nearest-Neighbor Matches}
\usage{
psm_match(score, treatment, M = 1L)
}
\arguments{
\item{score}{Numeric matching score, usually an estimated propensity score.}
\item{treatment}{Binary treatment vector coded as 0/1, logical, or
control/treated.}
\item{M}{Positive integer number of matches.}
}
\value{
A list containing the match matrix, all-target reuse counts, treated-target
reuse counts, score, treatment vector, and arm sizes.
}
\description{
Finds the \code{M} nearest opposite-treatment observations for every unit using
a scalar score. Matching is with replacement.
}
\examples{
dat <- simulate_psmqte_data(100, seed = 3)
mat <- psm_match(dat$p_true, dat$D, M = 2)
head(mat$match_matrix)
sum(mat$K_all)
}
