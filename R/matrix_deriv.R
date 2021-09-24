#' 1st order matrix derivative of a matrix function
#'
#' @description Let
#'   \eqn{X} be a matrix of size \eqn{r} by \eqn{c},
#'   \eqn{F} be a matrix function (wrt \eqn{X}) of size \eqn{p} by \eqn{q},
#'   then \eqn{dF/dX} is the matrix derivative of F of size \eqn{rp} by \eqn{cq}
#' @param func matrix function wrt matrix X.
#' @param X numeric matrix.
#' @param transpose_deriv_operator logical; apply transpose to the
#' derivative operator \eqn{d/dX}.
#' @inheritParams vector_deriv_3rd
#' @return matrix with \eqn{rp} rows and \eqn{cq} columns.
#' @seealso [vector_deriv_3rd], [vector_deriv_4th], [numDeriv::grad], [numDeriv::hessian]
#' @examples
#' func_11 <- function(X) X[1,1] * matrix(1:4,2,2)
#' func_12 <- function(X) X[1,2] * matrix(1:4,2,2)
#' func_21 <- function(X) X[2,1] * matrix(1:4,2,2)
#' func_22 <- function(X) X[2,2] * matrix(1:4,2,2)
#'
#' X <- matrix(1, 2, 2)
#' matrix_deriv(func_11, X)
#' ##      [,1] [,2] [,3] [,4]
#' ## [1,]    1    3    0    0
#' ## [2,]    2    4    0    0
#' ## [3,]    0    0    0    0
#' ## [4,]    0    0    0    0
#' matrix_deriv(func_12, X)
#' ##      [,1] [,2] [,3] [,4]
#' ## [1,]    0    0    1    3
#' ## [2,]    0    0    2    4
#' ## [3,]    0    0    0    0
#' ## [4,]    0    0    0    0
#' matrix_deriv(func_21, X)
#' ##      [,1] [,2] [,3] [,4]
#' ## [1,]    0    0    0    0
#' ## [2,]    0    0    0    0
#' ## [3,]    1    3    0    0
#' ## [4,]    2    4    0    0
#' matrix_deriv(func_22, X)
#' ##      [,1] [,2] [,3] [,4]
#' ## [1,]    0    0    0    0
#' ## [2,]    0    0    0    0
#' ## [3,]    0    0    1    3
#' ## [4,]    0    0    2    4
#' @export
matrix_deriv <- function(func, X, transpose_deriv_operator = FALSE,
                         method = "Richardson", method.args = list(), ...){
  f0 <- func(X)
  if(!is.matrix(f0)){
    stop("func(X) is not a matrix")
  }
  p <- nrow(f0)
  q <- ncol(f0)
  X <- as.matrix(X)
  m <- nrow(X)
  n <- ncol(X)
  zeros <- matrix(0,m,n)
  args <- list(eps = 1e-4, d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-07),
               r = 4, v = 2, show.details = FALSE)
  args[names(method.args)] <- method.args
  d <- args$d
  r <- args$r
  v <- args$v
  h <- max(abs(d*X), args$zero.tol*d)
  if(transpose_deriv_operator==FALSE){
    index_rowX <- matrix(1:m,m,n)%x%matrix(1,p,q)
    index_colX <- matrix(1:n,m,n,byrow=TRUE)%x%matrix(1,p,q)
    index_rowY <- matrix(1,m,n)%x%matrix(1:p,p,q)
    index_colY <- matrix(1,m,n)%x%matrix(1:q,p,q,byrow=TRUE)
  } else{
    index_rowX <- t(matrix(1:m,m,n))%x%matrix(1,p,q)
    index_colX <- t(matrix(1:n,m,n,byrow=TRUE))%x%matrix(1,p,q)
    index_rowY <- t(matrix(1,m,n))%x%matrix(1:p,p,q)
    index_colY <- t(matrix(1,m,n))%x%matrix(1:q,p,q,byrow=TRUE)
  }
  output <- array(NA, c(nrow(index_rowX), ncol(index_colX), r))
  for(k in 1:r){
    for(i in 1:nrow(index_rowX)){
      for(j in 1:ncol(index_colX)){
        e <- zeros
        e[index_rowX[i,j],index_colX[i,j]] <- 1
        output[i,j,k] <-
          (func(X+h*e)[index_rowY[i,j],index_colY[i,j]] +
            -func(X-h*e)[index_rowY[i,j],index_colY[i,j]])/(2*h)
      }
    }
    h <- h/v
  }
  for (m in 1:(k - 1)) {
    output <- (output[,, 2:(k + 1 - m), drop = FALSE] * (4^m) -
                output[,, 1:(k - m), drop = FALSE])/(4^m - 1)
  }
  return(output[,,1])
}
