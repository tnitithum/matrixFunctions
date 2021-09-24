#' Check if Matrix is invertible
#'
#' Checks if the matrix is numerically invertible.
#' @inheritParams MASS::ginv
#' @param X numeric square matrix.
#' @return TRUE if numerically invertible, FALSE otherwise.
#' @examples
#' X <- matrix(c(1,1,0,0), 2, 2)
#' is_invertible(X)
#'
#' Y <- diag(2)
#' is_invertible(Y)
#' @export
is_invertible <- function(X, tol=sqrt(.Machine$double.eps)){
	Xsvd <- svd(X)
	Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
	return(all(Positive))
}
