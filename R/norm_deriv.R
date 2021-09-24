#' First Order Derivative of Euclidean norm
#'
#' First Order Derivative of Euclidean norm of vector \code{x} with respect
#' to \code{x}.
#' @param x numeric vector.
#' @return numeric vector of length x.
#' @seealso [D_norm_2nd], [D_norm_3rd], [D_norm_4th]
#' @examples
#' library(numDeriv)
#' func <- function(x) norm(x, type="2")
#' x <- 1:2
#' print(numDeriv::grad(func, x))
#' ## [1] 0.4472136 0.8944272
#' print(D_norm_1st(x))
#' ## [1] 0.4472136 0.8944272
#' @export
D_norm_1st <- function(x){
  f <- function(x) norm(x, type="2")
  return(x/f(x))
}

#' Second Order Derivative of Euclidean norm
#'
#' Second Order Derivative of Euclidean norm of vector \code{x} with respect
#' to \code{x}.
#' @inheritParams D_norm_1st
#' @return numeric vector of length x.
#' @seealso [D_norm_1st], [D_norm_3rd], [D_norm_4th]
#' @examples
#' library(numDeriv)
#' func <- function(x) norm(x, type="2")
#' x <- 1:2
#' print(numDeriv::hessian(func, x))
#' ##           [,1]        [,2]
#' ##[1,]  0.3577709 -0.17888544
#' ##[2,] -0.1788854  0.08944272
#' print(D_norm_2nd(x))
#' ##           [,1]        [,2]
#' ##[1,]  0.3577709 -0.17888544
#' ##[2,] -0.1788854  0.08944272
#' @export
D_norm_2nd <- function(x){
  f <- function(x) norm(x, type="2")
  d <- length(x)
  return((f(x)^2*diag(d) - x%*%t(x))/f(x)^3)
}

#' Third Order Derivative of Euclidean norm
#'
#' Third Order Derivative of Euclidean norm of vector \code{x} with respect
#' to \code{x}.
#' @inheritParams D_norm_1st
#' @return numeric vector of length x.
#' @seealso [D_norm_1st], [D_norm_2nd], [D_norm_4th]
#' @examples
#' func <- function(x) norm(x, type="2")
#' x <- 1:2
#' print(vector_deriv_3rd(func, x))
#' ##             [,1]        [,2]
#' ## [1,] -0.21466254 -0.07155418
#' ## [2,] -0.07155418  0.12521981
#' ## [3,] -0.07155418  0.12521981
#' ## [4,]  0.12521981 -0.10733128
#' print(D_norm_3rd(x))
#' ##             [,1]        [,2]
#' ## [1,] -0.21466253 -0.07155418
#' ## [2,] -0.07155418  0.12521981
#' ## [3,] -0.07155418  0.12521981
#' ## [4,]  0.12521981 -0.10733126
#' @export
D_norm_3rd <- function(x){
  f <- function(x) norm(x, type="2")
  d <- length(x)
  -(f(x)^2*{x}%x%diag(d) +
    f(x)^2*c(diag(d))%*%t(x) +
    f(x)^2*diag(d)%x%{x} +
    -3*{x}%x%{{x}%*%t(x)}
  )/f(x)^5
}

#' Fourth Order Derivative of Euclidean norm
#'
#' Fourth Order Derivative of Euclidean norm of vector \code{x} with respect
#' to \code{x}.
#' @inheritParams D_norm_1st
#' @return numeric vector of length x.
#' @seealso [D_norm_1st], [D_norm_2nd], [D_norm_3rd]
#' @examples
#' func <- function(x) norm(x, type="2")
#' x <- 1:2
#' print(vector_deriv_4th(func, x))
#' ##               [,1]        [,2]        [,3]        [,4]
#' ## [1,]  5.156024e-07  0.21466240  0.21466240 -0.03577708
#' ## [2,]  2.146624e-01 -0.03577708 -0.03577708 -0.10733140
#' ## [3,]  2.146624e-01 -0.03577708 -0.03577708 -0.10733140
#' ## [4,] -3.577708e-02 -0.10733140 -0.10733140  0.16099882
#' print(D_norm_4th(x))
#' ##               [,1]        [,2]        [,3]        [,4]
#' ## [1,] -2.081668e-17  0.21466253  0.21466253 -0.03577709
#' ## [2,]  2.146625e-01 -0.03577709 -0.03577709 -0.10733126
#' ## [3,]  2.146625e-01 -0.03577709 -0.03577709 -0.10733126
#' ## [4,] -3.577709e-02 -0.10733126 -0.10733126  0.16099689
#' @export
D_norm_4th <- function(x){
  f <- function(x) norm(x, type="2")
  f0 <- f(x)
  d <- length(x)
  term_3 <- -f0^(-3)*(
    commutation_matrix(d) +
    diag(d^2) +
    c(diag(d))%*%t(c(diag(d)))
  )
  term_5 <- 3/f0^5*(
    (x)%*%t(x)%x%diag(d) +
    diag(d)%x%(x%*%t(x)) +
    c(diag(d))%x%t(x)%x%t(x) +
    (x)%x%diag(d)%x%t(x) +
    t(x)%x%diag(d)%x%x +
    (x)%x%(x)%x%t(c(diag(d)))
  )
  term_7 <- -15/f0^7*
    (x)%x%(x)%x%t(x)%x%t(x)
  return(term_3 + term_5 + term_7)
}
