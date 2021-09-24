# Defining Functions ------------------------------------------------------
library(numDeriv)

vector_deriv_1st <- function(...){
  numDeriv::grad(...)
}

vector_deriv_2nd <- function(...) {
  numDeriv::hessian(...)
}


# Vector Polynomial (4th power case) --------------------------------------
test_that("vector polynomial (4th power case)", {
  func <- function(x) as.matrix(t(x)%*%x%*%t(x)%*%x)
  x <- 1

  # expect_equal(vector_deriv_1st(func, x), 4)
  # expect_equal(as.numeric(vector_deriv_2nd(func, x)), 12)
  expect_equal(as.numeric(vector_deriv_3rd(func, x)), 24)
  expect_equal(as.numeric(vector_deriv_4th(func, x)), 24)
})


# sin ---------------------------------------------------------------------
test_that("sin(sum(x))", {
  func <- function(x) sin(sum(x))
  x <- 0

  # expect_equal(vector_deriv_1st(func, x), 1)
  # expect_equal(as.numeric(vector_deriv_2nd(func, x)), 0)
  expect_equal(as.numeric(vector_deriv_3rd(func, x)), -1)
  expect_equal(as.numeric(vector_deriv_4th(func, x)), 0)
})


# Norm --------------------------------------------------------------------
test_that("norm", {
  func <- function(x) as.matrix(norm(x, type="2"))
  x <- c(1,5)

  expect_equal(vector_deriv_1st(func, x), D_norm_1st(x))
  expect_equal(vector_deriv_2nd(func, x), D_norm_2nd(x))
  expect_equal(vector_deriv_3rd(func, x), D_norm_3rd(x))
  expect_equal(vector_deriv_4th(func, x), D_norm_4th(x))
})


# Spherical Moments -------------------------------------------------------
test_that("spherical moments", {
  r <- 1.1
  func <- function(x) spherical_moment_0th(x, r)
  x <- c(1,3)

  expect_equal(vector_deriv_1st(func, x), spherical_moment_1st(x, r))
  expect_equal(vector_deriv_2nd(func, x), spherical_moment_2nd(x, r))
  expect_equal(vector_deriv_3rd(func, x), spherical_moment_3rd(x, r))
  expect_equal(vector_deriv_4th(func, x), spherical_moment_4th(x, r))
})
