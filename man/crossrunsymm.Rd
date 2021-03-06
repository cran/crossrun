% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossrunsymm.r
\name{crossrunsymm}
\alias{crossrunsymm}
\title{Joint Probabilities for Crossings and Runs, Symmetric Case}
\usage{
crossrunsymm(nmax = 100, mult = 2, prec = 120, printn = FALSE)
}
\arguments{
\item{nmax}{; max sequence length.}

\item{mult}{; multiplier for joint probabilities. Default 2.}

\item{prec}{; mpft precision.}

\item{printn}{; logical for including progress output.}
}
\value{
pt, list of joint probabilities, multiplied with \eqn{m^{n-1}}.
In addition cumulative probabilities qt within each row are also included.
}
\description{
Joint probability distribution for the number of crossings
C and the longest run L in a sequence of n independent Bernoulli observations
with success probability p. To enhance precision, results are stored
in mpfr arrays and the probabilities are multiplied by \eqn{m^{n-1}}
for a multiplier m. This is for the symmetric case with success
probability 0.5, in which the multiplied probabilities are
integers for the default value 2 of the multiplier.
}
\examples{
crs10 <- crossrunsymm(nmax=10,printn=TRUE)
}
