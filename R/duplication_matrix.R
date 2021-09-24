#' Duplication matrix D_m
#'
#' @description Let
#'   \eqn{A} be a \eqn{m} by \eqn{m} symmetric matrix, then
#'   \eqn{D_m} defines a duplication matrix such that
#'   \eqn{vec(A) = D_m vech(A)}
#' @param m dimension parameter of \eqn{A}.
#' @return matrix with \eqn{m^2} rows and \eqn{m(m+1)/2} columns.
#' @keywords matrix function
#' @seealso \link{elimination_matrix}
#' @examples
#' duplication_matrix(2)
#' ##      [,1] [,2] [,3]
#' ## [1,]    1    0    0
#' ## [2,]    0    1    0
#' ## [3,]    0    1    0
#' ## [4,]    0    0    1
#'
#' # verifying vec(A) = D_m vech(A) for m=2 case
#' m <- 2
#' S <- symm_mat(1:3)
#' print(c(S))
#' ## [1] 1 2 2 3
#' print(c(duplication_matrix(m) %*% vech(S)))
#' ## [1] 1 2 2 3
#' @export
duplication_matrix <- function(m){
  diag(m*{m+1}/2)[c(symm_mat(m=m)),]
}
