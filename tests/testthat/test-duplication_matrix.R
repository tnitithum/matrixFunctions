test_that("verifying vec(A) = D_m vech(A)", {
  # m=2 case
  m <- 2
  S <- symm_mat(1:3)
  expect_equal(c(S),
               c(duplication_matrix(m) %*% vech(S)))

  # m=10 case
  m <- 10
  S <- symm_mat(m=10) #symm_mat(1:(10*11/2))
  expect_equal(c(S),
               c(duplication_matrix(m) %*% vech(S)))
})
