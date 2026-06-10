\name{m_grid_diverging}
\alias{m_grid_diverging}
\title{Diverging-M Grid Used in Matching Simulations}
\usage{
m_grid_diverging(n, max_M = NULL)
}
\arguments{
\item{n}{Sample size.}
\item{max_M}{Optional upper bound.}
}
\value{
An integer vector of matching counts.
}
\description{
Returns \code{M = 2^j} for \code{0 <= j <= floor(log2(n) / 2)}.
}
\examples{
m_grid_diverging(512)
m_grid_diverging(512, max_M = 8)
}
