#' 3rd order vector derivative of a scalar function
#'
#' @description
#' Calculate a numerical approximation to the
#' 3rd order vector derivative of a scalar function at a parameter value.
#' @inheritParams numDeriv::grad
#' @param func scalar function that takes vector input \code{x}.
#' @param method currently only \code{"Richardson"} is implemented for the approximation.
#' @param side an indication of whether one-sided derivatives should be attempted.
#' @param output either \code{"matrix"} or \code{"array"}.
#' @return numeric square matrix or 3-D array depending on \code{output} argument.
#' @seealso [vector_deriv_4th], [matrix_deriv], [numDeriv::grad], [numDeriv::hessian]
#' @examples
#' f <- function(x) norm(x, type = "2")
#' x <- 1:3
#' vector_deriv_3rd(f, x)
#' vector_deriv_3rd(f, x, output = "array")
#' @export
vector_deriv_3rd <- function (func, x, method = "Richardson", side = NULL,
                              method.args = list(), output = "matrix", ...){
  f <- func(x, ...)
  n <- length(x)
  if (is.null(side))
    side <- rep(NA, n)
  else {
    if (n != length(side))
      stop("Non-NULL argument 'side' should have the same length as x")
    if (any(1 != abs(side[!is.na(side)])))
      stop("Non-NULL argument 'side' should have values NA, +1, or -1.")
  }
  if (1 != length(f))
    stop("func must be a scalar valued function.")
  if (method == "Richardson") {
    args <- list(eps = 0.1, d = 0.3, zero.tol = sqrt(.Machine$double.eps/7e-07),
                 r = 4, v = 2, show.details = FALSE)
    args[names(method.args)] <- method.args
    d <- args$d
    r <- args$r
    v <- args$v
    f0 <- func(x, ...)
    show.details <- args$show.details
    a_forward <- array(NA, c(n,n,n,r))
    a_backward <- array(NA, c(n,n,n,r))
    a <- array(NA, c(n,n,n,r))
    ones <- rep(1, n)
    h <- rep(min(abs(d * x) + args$eps * (abs(x) < args$zero.tol)), n)
    pna <- (side == 1) & !is.na(side)
    mna <- (side == -1) & !is.na(side)
    for (s in 1:r) {
      ph <- mh <- h
      ph[pna] <- 2 * ph[pna]
      ph[mna] <- 0
      mh[mna] <- 2 * mh[mna]
      mh[pna] <- 0
      for (i in 1:n) {
        for(j in 1:n) {
          for(k in 1:n) {
            if ((s != 1) && (abs(a[i,j,k,s-1]) < 1e-20)){
              a[i,j,k,s] <- 0
            } else{
              func1 <- function(x, ...) func(x + ph*(i==seq(n)), ...) -  func(x, ...)
              func2 <- function(x, ...) func1(x + ph*(j==seq(n)), ...) - func1(x, ...)
              a_forward[i,j,k,s]  <-    func2(x + ph*(k==seq(n)), ...) - func2(x, ...)
              func_1 <- function(x, ...) func(x, ...) -   func(x - mh*(i==seq(n)), ...)
              func_2 <- function(x, ...) func_1(x, ...) - func_1(x - mh*(j==seq(n)), ...)
              a_backward[i,j,k,s] <-     func_2(x, ...) - func_2(x - mh*(k==seq(n)), ...)
              a[i,j,k,s] <- (a_forward[i,j,k,s] + a_backward[i,j,k,s])/(2 * h[i]^3)
            }
          }
        }
      }
      if (any(is.na(a[,,,s])))
        stop("function returns NA at ", h, " distance from x.")
      h <- h/v
    }
    if (show.details) {
      cat("\n", "first order approximations",
          "\n")
      print(a, 12)
    }
    for (m in 1:(s - 1)) {
      a <- (a[,,,2:(s + 1 - m), drop = FALSE] * (4^m) -
              a[,,, 1:(s - m), drop = FALSE])/(4^m - 1)
      if (show.details & m != (s - 1)) {
        cat("\n", "Richarson improvement group No. ",
            m, "\n")
        print(a[,,, 1:(s - m), drop = FALSE], 12)
      }
    }
    if(output=="matrix"){
      return(matrix(a[,,,1], n^2, n))
    }
    if(output=="array"){
      return(a[,,,1])
    }
  }
  else stop("indicated method ", method, "not supported.")
}


#' 4th order vector derivative of a scalar function
#'
#' @description
#' Calculate a numerical approximation to the
#' 4th order vector derivative of a scalar function at a parameter value.
#' @inheritParams vector_deriv_3rd
#' @return numeric square matrix or 4-D array depending on \code{output} argument.
#' @examples
#' f <- function(x) norm(x, type = "2")
#' x <- 1:3
#' vector_deriv_4th(f, x)
#' vector_deriv_4th(f, x, output = "array")
#' @export
vector_deriv_4th <- function (func, x, method = "Richardson", side = NULL,
                              method.args = list(), output = "matrix", ...)
{
  f <- func(x, ...)
  n <- length(x)
  if (is.null(side))
    side <- rep(NA, n)
  else {
    if (n != length(side))
      stop("Non-NULL argument 'side' should have the same length as x")
    if (any(1 != abs(side[!is.na(side)])))
      stop("Non-NULL argument 'side' should have values NA, +1, or -1.")
  }
  if (1 != length(f))
    stop("func must be a scalar valued function.")
  if (method == "Richardson") {
    args <- list(eps = 0.1, d = 0.35, zero.tol = sqrt(.Machine$double.eps/7e-07),
                 r = 4, v = 2, show.details = FALSE)
    args[names(method.args)] <- method.args
    d <- args$d
    r <- args$r
    v <- args$v
    f0 <- func(x, ...)
    show.details <- args$show.details
    a_forward <- array(NA, c(n,n,n,n,r))
    a_backward <- array(NA, c(n,n,n,n,r))
    a <- array(NA, c(n,n,n,n,r))
    ones <- rep(1, n)
    h <- rep(min(abs(d * x) + args$eps * (abs(x) < args$zero.tol)), n)
    pna <- (side == 1) & !is.na(side)
    mna <- (side == -1) & !is.na(side)
    for (s in 1:r) {
      ph <- mh <- h
      ph[pna] <- 2 * ph[pna]
      ph[mna] <- 0
      mh[mna] <- 2 * mh[mna]
      mh[pna] <- 0
      for (i in 1:n) {
        for(j in 1:n) {
          for(k in 1:n) {
            for(l in 1:n){
              if ((s != 1) && (abs(a[i,j,k,l,s-1]) < 1e-20)){
                a[i,j,k,l,s] <- 0
              } else{
                func1 <- function(x, ...) func(x + ph*(i==seq(n)), ...) -  func(x, ...)
                func2 <- function(x, ...) func1(x + ph*(j==seq(n)), ...) - func1(x, ...)
                func3 <- function(x, ...) func2(x + ph*(k==seq(n)), ...) - func2(x, ...)
                a_forward[i,j,k,l,s]  <-    func3(x + ph*(l==seq(n)), ...) - func3(x, ...)
                func_1 <- function(x, ...) func(x, ...) -   func(x - mh*(i==seq(n)), ...)
                func_2 <- function(x, ...) func_1(x, ...) - func_1(x - mh*(j==seq(n)), ...)
                func_3 <- function(x, ...) func_2(x, ...) - func_2(x - mh*(k==seq(n)), ...)
                a_backward[i,j,k,l,s] <-     func_3(x, ...) - func_3(x - mh*(l==seq(n)), ...)
                a[i,j,k,l,s] <- (a_forward[i,j,k,l,s] + a_backward[i,j,k,l,s])/(2 * h[i]^4)
              }
            }
          }
        }
      }
      if (any(is.na(a[,,,,s])))
        stop("function returns NA at ", h, " distance from x.")
      h <- h/v
    }
    if (show.details) {
      cat("\n", "first order approximations",
          "\n")
      print(a, 12)
    }
    for (m in 1:(s - 1)) {
      a <- (a[,,,,2:(s + 1 - m), drop = FALSE] * (4^m) -
              a[,,,, 1:(s - m), drop = FALSE])/(4^m - 1)
      if (show.details & m != (s - 1)) {
        cat("\n", "Richarson improvement group No. ",
            m, "\n")
        print(a[,,,, 1:(s - m), drop = FALSE], 12)
      }
    }
    if(output=="matrix"){
      return(matrix(a[,,,,1], n^2, n^2))
    }
    if(output=="array"){
      return(a[,,,,1])
    }
  }
  else stop("indicated method ", method, "not supported.")
}
