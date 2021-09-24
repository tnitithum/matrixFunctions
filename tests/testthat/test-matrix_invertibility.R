test_that("is_invertible: examples", {
  X <- matrix(c(1,1,0,0), 2, 2)
  expect_false(is_invertible(X))

  Y <- diag(2)
  expect_true(is_invertible(Y))
})
