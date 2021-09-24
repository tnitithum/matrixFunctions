#' Commutation matrix K
#'
#' @description Let
#'   \eqn{A} be a \eqn{r} by \eqn{c} matrix, then
#'   \eqn{K} defines a commutation matrix such that
#'   \eqn{vec(A) = K vec(A')} where prime denotes transpose.
#' @param r positive integer for row dimension.
#' @param c positive integer for column dimension.
#' @param checks logical, apply checks to arguments.
#' @return matrix with \eqn{rc} rows and \eqn{rc} columns.
#' @keywords matrix function
#' @examples
#' r = 2
#' c = 2
#' K <- commutation_matrix(r, c)
#' ##      [,1] [,2] [,3] [,4]
#' ## [1,]    1    0    0    0
#' ## [2,]    0    0    1    0
#' ## [3,]    0    1    0    0
#' ## [4,]    0    0    0    1
#'
#' # verifying vec(A) = K vec(A') for r=c=2 case
#' A <- matrix(1:(r*c), r, c)
#' vecA <- c(A)
#' vecAt <- c(t(A))
#' print(vecA)
#' ## [1] 1 2 3 4
#' print(c(K %*% vecAt))
#' ## [1] 1 2 3 4
#' @export
commutation_matrix <- function(r, c=r, checks=TRUE){
  if(checks){
    if((r%%1) > 0)
      stop(paste0("argument r:", r, " is not integer"))
    if((c%%1) > 0)
      stop(paste0("argument c:", c, " is not integer"))
  }
  perm_rows.v <- NULL
  for(i in 1:r){
    perm_rows.v <- c(perm_rows.v, seq.int(i, by=r, length.out=c))
  }
  return(diag(r*c)[perm_rows.v,])
}
