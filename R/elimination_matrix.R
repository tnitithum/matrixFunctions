#' Elimination matrix L_m
#'
#' @description Let
#'   \eqn{A} be a \eqn{m} by \eqn{m} symmetric matrix, then
#'   \eqn{L_m} defines an elimination matrix such that
#'   \eqn{vech(A) = L_m vec(A)}
#' @inheritParams duplication_matrix
#' @return matrix with \eqn{m(m+1)/2} rows and \eqn{m^2} columns.
#' @keywords matrix function
#' @seealso \link{duplication_matrix}
#' @examples
#' elimination_matrix(2)
#' ##      [,1] [,2] [,3] [,4]
#' ## [1,]    1    0    0    0
#' ## [2,]    0    1    0    0
#' ## [3,]    0    0    0    1
#'
#' # verifying vech(A) = L_m vec(A) for m=2 case
#' m <- 2
#' A <- matrix(1:m^2, m, m)
#' print(vech(A))
#' ## [1] 1 2 4
#' print(c(elimination_matrix(m) %*% c(A)))
#' ## [1] 1 2 4
#' @export
elimination_matrix <- function(m){
  M <- matrix(1:m^2, m, m)
  output <- matrix(0, m*{m+1}/2, m^2)
  output[,1:m^2 %in% vech(M)] <- diag(m*{m+1}/2)
  return(output)
}
