\docType{package}
\name{psmqte-package}
\alias{psmqte-package}
\title{Propensity Score Matching for QTE, QTT, ATE, and ATT}
\description{
The \pkg{psmqte} package implements propensity-score matching estimators for
quantile treatment effects, quantile treatment effects on the treated, average
treatment effects, and average treatment effects on the treated.
}
\details{
The main function is \code{\link{psmqte}}, which estimates QTE and QTT curves
using nearest-neighbor matching with replacement on an estimated or supplied
propensity score. The companion function \code{\link{psmate}} estimates ATE and
ATT from the same matching design.
}
\author{
psmqte authors
}
\references{
Abadie, A. and Imbens, G. W. (2016). Matching on the estimated propensity
score. \emph{Econometrica}, 84(2), 781-807.

He, Y. and Han, F. (2024). On propensity score matching with a diverging number
of matches. \emph{Biometrika}, 111(4), 1421-1428.
}
