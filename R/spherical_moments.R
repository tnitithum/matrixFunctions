#' Zeroth order spherical moment
#'
#' Spherical integral of \eqn{exp(s'x)} where prime denotes transpose.
#' @param x numeric vector, exponent in the exponential function of the integrand.
#' @param r positive numeric, radius of the spherical integral.
# #' @details
# #' For the 1-dimensional case, it's equivalent to evaluating the integrand at
# #' the endpoints \{-r,r\} which reduces the integral to \eqn{2*cosh(r*|x|)}.
#' @return Scalar.
#' @seealso [spherical_moment_1st], [spherical_moment_2nd], [spherical_moment_3rd], [spherical_moment_4th]
#' @examples
#' spherical_moment_0th(1:3, 2)
#' @export
spherical_moment_0th <- function(x, r=1){
  d <- length(x)
  a <- norm(x, type="2")
  return((2*pi*r)^(d/2) * a^(1-d/2) * besselI(r*a, d/2-1))
}

#' First order spherical moments
#'
#' First order derivative of Sperical integral of \eqn{exp(s'x)} where prime denotes transpose.
#' @inheritParams spherical_moment_0th
#' @return Numeric vector
#' @seealso [spherical_moment_0th], [spherical_moment_2nd], [spherical_moment_3rd], [spherical_moment_4th]
#' @examples
#' spherical_moment_1st(1:3, 2)
#' @export
spherical_moment_1st <- function(x, r=1){
  d <- length(x)
  a <- norm(x, type="2")
  return((2*pi*r)^(d/2) * r * a^(-d/2) * besselI(r*a, d/2) * x)
}

#' Second order spherical moments
#'
#' Second order derivative of Sperical integral of \eqn{exp(s'x)} where prime denotes transpose.
#' @inheritParams spherical_moment_0th
#' @return Numeric matrix
#' @seealso [spherical_moment_0th], [spherical_moment_1st], [spherical_moment_3rd], [spherical_moment_4th]
#' @examples
#' spherical_moment_2nd(1:3, 2)
#' @export
spherical_moment_2nd <- function(x, r=1){
  d <- length(x)
  a <- norm(x, type="2")
  return((2*pi*r)^(d/2) * r * a^(-d/2-1) * (
    a*besselI(r*a, d/2)*diag(d) +
    r*besselI(r*a, d/2+1)*x%*%t(x)
  ))
}

#' Third order spherical moments
#'
#' Third order derivative of Sperical integral of \eqn{exp(s'x)} where prime denotes transpose.
#' @inheritParams spherical_moment_0th
#' @return Numeric matrix
#' @seealso [spherical_moment_0th], [spherical_moment_1st], [spherical_moment_2nd], [spherical_moment_4th]
#' @examples
#' spherical_moment_3rd(1:3, 2)
#' @export
spherical_moment_3rd <- function(x, r=1){
  d <- length(x)
  a <- norm(x, type="2")
  C3 <- r*a*x%x%x%x%t(x)
  D3 <- a^2*(c(diag(d))%*%t(x) + diag(d)%x%x + x%x%diag(d)) +
    -{d+2}*x%x%x%x%t(x)
  return((2*pi*r)^(d/2) * r^2*a^{-d/2-3}*(
    besselI(r*a,d/2)*C3 +
    besselI(r*a,d/2+1)*D3
  ))
}

#' Forth order spherical moments
#'
#' Forth order derivative of Sperical integral of \eqn{exp(s'x)} where prime denotes transpose.
#' @inheritParams spherical_moment_0th
#' @return Numeric matrix
#' @seealso [spherical_moment_0th], [spherical_moment_1st], [spherical_moment_2nd], [spherical_moment_3rd]
#' @examples
#' spherical_moment_4th(1:3, 2)
#' @export
spherical_moment_4th <- function(x, r=1){
  d <- length(x)
  a <- norm(x, type="2")
  A4 <- x%x%x%x%t(x)%x%t(x)
  A2 <- diag(d)%x%x%x%t(x) + x%x%diag(d)%x%t(x) + x%x%x%x%t(c(diag(d)))
  B2 <- t(x)%x%{c(diag(d))%x%t(x) + diag(d)%x%x + x%x%diag(d)}
  B0 <- c(diag(d))%x%t(c(diag(d))) + matrixFunctions::commutation_matrix(d) + diag(d^2)
  C4 <- r*a*{
    -(d+4)*A4 + a^2*{A2 + B2}
  }
  D4 <- -(d+2)*a^2*A2 +
    ((d+4)*(d+2) + r^2*a^2)*A4 +
    a^4*B0 +
    -{d+2}*a^2*B2
  return((2*pi*r)^(d/2) * r^2 * a^(-d/2-5) * (
    besselI(r*a,d/2)*C4 +
    besselI(r*a,d/2+1)*D4
  ))
}
