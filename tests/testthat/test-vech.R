test_that("verifying definition", {
  # m=2 case
  A <- symm_mat(m=2)
  expect_equal(vech(A),
               c(1,2,3))

  # m=10 case
  m = 10
  A <- symm_mat(m=m)
  expect_equal(vech(A),
               1:(m*(m+1)/2))
})
