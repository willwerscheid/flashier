context("low-rank matrix representations")

set.seed(666)

n <- 20
p <- 100

LF1 <- outer(rep(1, n), rep(1, p))
LL <- rnorm(n)
FF <- rnorm(p)
LF2 <- 5 * outer(LL, FF)
LF <- LF1 + LF2
M <- LF + 0.1 * rnorm(n * p)

test_that("Using a low-rank matrix representation works", {
  lr1 <- svd(M, nu = 3, nv = 3)
  f1  <- flash(lr1, verbose = 0)
  expect_equal(fitted(f1), LF1 + LF2, tol = 0.25, scale = 1)

  lr2 <- irlba::irlba(M, 3)
  f2  <- flash(lr2, verbose = 0)
  expect_equal(fitted(f2), LF1 + LF2, tol = 0.25, scale = 1)
})
