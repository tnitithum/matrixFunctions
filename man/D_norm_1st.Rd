% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm_deriv.R
\name{D_norm_1st}
\alias{D_norm_1st}
\title{First Order Derivative of Euclidean norm}
\usage{
D_norm_1st(x)
}
\arguments{
\item{x}{numeric vector.}
}
\value{
numeric vector of length x.
}
\description{
First Order Derivative of Euclidean norm of vector \code{x} with respect
to \code{x}.
}
\examples{
library(numDeriv)
func <- function(x) norm(x, type="2")
x <- 1:2
print(numDeriv::grad(func, x))
## [1] 0.4472136 0.8944272
print(D_norm_1st(x))
## [1] 0.4472136 0.8944272
}
\seealso{
\link{D_norm_2nd}, \link{D_norm_3rd}, \link{D_norm_4th}
}
