#' Half-vectorization
#'
#' Vectorising the lower triangular part of a symmetric matrix
#' @param A symmetric matrix.
#' @return vector with elements from the lower triangular matrix.
#' @seealso [symm_mat]
#' @examples
#' A <- symm_mat(m=2)
#' ##      [,1] [,2]
#' ## [1,]    1    2
#' ## [2,]    2    3
#'
#' vech(A)
#' ## [1] 1 2 3
#' @export
vech <- function(A){
  A[!upper.tri(A)]
}
