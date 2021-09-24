#' Symmetric matrix transform
#'
#' @description
#' Inverse of half-vectorisation function \code{vech}. In other words,
#' it creates a symmetric matrix using entries from the upper triangular matrix.
#' @param a vector of the upper matrix.
#' @param m dimension parameter of symmetric matrix.
#' @seealso [vech]
#' @examples
#' symm_mat(1:3)
#' ##      [,1] [,2]
#' ## [1,]    1    2
#' ## [2,]    2    3
#'
#' symm_mat(m=2)
#' ##      [,1] [,2]
#' ## [1,]    1    2
#' ## [2,]    2    3
#' @export
symm_mat <- function(a, m){
  if(missing(a) & !missing(m)){
    a <- 1:(m*(m+1)/2)
  } else{
    m <- {-1+sqrt(1+8*length(a))}/2
  }
  S <- matrix(NA,m,m)
  S[!upper.tri(S)] <- a
  S <- t(S)
  S[!upper.tri(S)] <- a
  return(S)
}
