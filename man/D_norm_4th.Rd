% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm_deriv.R
\name{D_norm_4th}
\alias{D_norm_4th}
\title{Fourth Order Derivative of Euclidean norm}
\usage{
D_norm_4th(x)
}
\arguments{
\item{x}{numeric vector.}
}
\value{
numeric vector of length x.
}
\description{
Fourth Order Derivative of Euclidean norm of vector \code{x} with respect
to \code{x}.
}
\examples{
func <- function(x) norm(x, type="2")
x <- 1:2
print(vector_deriv_4th(func, x))
##               [,1]        [,2]        [,3]        [,4]
## [1,]  5.156024e-07  0.21466240  0.21466240 -0.03577708
## [2,]  2.146624e-01 -0.03577708 -0.03577708 -0.10733140
## [3,]  2.146624e-01 -0.03577708 -0.03577708 -0.10733140
## [4,] -3.577708e-02 -0.10733140 -0.10733140  0.16099882
print(D_norm_4th(x))
##               [,1]        [,2]        [,3]        [,4]
## [1,] -2.081668e-17  0.21466253  0.21466253 -0.03577709
## [2,]  2.146625e-01 -0.03577709 -0.03577709 -0.10733126
## [3,]  2.146625e-01 -0.03577709 -0.03577709 -0.10733126
## [4,] -3.577709e-02 -0.10733126 -0.10733126  0.16099689
}
\seealso{
\link{D_norm_1st}, \link{D_norm_2nd}, \link{D_norm_3rd}
}
