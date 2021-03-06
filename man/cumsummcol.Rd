% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary.r
\name{cumsummcol}
\alias{cumsummcol}
\title{Column-Wise Cumulative Sums}
\usage{
cumsummcol(mtrx)
}
\arguments{
\item{mtrx}{mpfr two-dimensional array.}
}
\value{
mpfr array with column-wise cumulative sums,
same dimension as the original array.
}
\description{
Column-wise cumulative sums in
mpfr array.
}
\examples{
nill <- Rmpfr::mpfr(0, 120)
one <- Rmpfr::mpfr(1, 120)
two <- Rmpfr::mpfr(2, 120)
contents <- c(one,nill,nill, one,one,one, two,two,two)
mtrx3 <- Rmpfr::mpfr2array(contents, dim = c(3, 3))
print(mtrx3)
print(cumsummcol(mtrx3))
}
