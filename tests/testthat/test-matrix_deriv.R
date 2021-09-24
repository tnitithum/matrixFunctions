test_that("matrix_deriv: examples", {
  func_11 <- function(X) X[1,1] * matrix(1:4,2,2)
  func_12 <- function(X) X[1,2] * matrix(1:4,2,2)
  func_21 <- function(X) X[2,1] * matrix(1:4,2,2)
  func_22 <- function(X) X[2,2] * matrix(1:4,2,2)
  X <- matrix(1, 2, 2)

  expect_equal(matrix_deriv(func_11, X),
               matrix(c(1,3,0,0,
                        2,4,0,0,
                        0,0,0,0,
                        0,0,0,0),
                      4, 4, byrow=TRUE))
  expect_equal(matrix_deriv(func_21, X),
               matrix(c(0,0,0,0,
                        0,0,0,0,
                        1,3,0,0,
                        2,4,0,0),
                      4, 4, byrow=TRUE))
  expect_equal(matrix_deriv(func_12, X),
               matrix(c(0,0,1,3,
                        0,0,2,4,
                        0,0,0,0,
                        0,0,0,0),
                      4, 4, byrow=TRUE))
  expect_equal(matrix_deriv(func_22, X),
               matrix(c(0,0,0,0,
                        0,0,0,0,
                        0,0,1,3,
                        0,0,2,4),
                      4, 4, byrow=TRUE))
})


test_that("matrix_deriv: norm 4th order derivative", {
  x <- 1:2
  expect_equal(matrix_deriv(D_norm_3rd, x, transpose_deriv_operator = TRUE),
               D_norm_4th(x))
})
