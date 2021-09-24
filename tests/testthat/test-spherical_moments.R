test_that("1-dim case", {
  x <- 2
  r <- 1.1
  a <- abs(x)
  expect_equal(spherical_moment_0th(x, r),
               2*cosh(r*a))
  expect_equal(spherical_moment_1st(x, r),
               2*r*sinh(r*a))
  expect_equal(as.numeric(spherical_moment_2nd(x, r)),
               2*r^2*cosh(r*a))
  expect_equal(as.numeric(spherical_moment_3rd(x, r)),
               2*r^3*sinh(r*a))
  expect_equal(as.numeric(spherical_moment_4th(x, r)),
               2*r^4*cosh(r*a))
})

test_that("spherical_moment_0th: compare with numeric integral", {
  x <- 1:2
  r <- 1.1
  d <- length(x)
  a <- norm(x, type="2")
  expect_equal(spherical_moment_0th(x, r),
               2 * r^(d-1) * pi^((d-1)/2) / gamma((d-1)/2) *
                 integrate(function(phi) exp(r*a*cos(phi)) * sin(phi)^(d-2), 0, pi)$value)
})
