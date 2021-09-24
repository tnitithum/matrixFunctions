test_that("check example", {
  expect_equal(symm_mat(1:3),
               matrix(c(1,2,2,3), 2, 2))
})
