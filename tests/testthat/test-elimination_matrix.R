test_that("verifying vech(A) = L_m vec(A)", {
  m <- 2
  A <- matrix(1:m^2,m,m)
  expect_equal(vech(A),
               c(elimination_matrix(m) %*% c(A)))

  m <- 10
  A <- symm_mat(1:(10*11/2))
  expect_equal(vech(A),
               c(elimination_matrix(m) %*% c(A)))
})
