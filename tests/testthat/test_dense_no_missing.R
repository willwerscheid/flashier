context("flashr comparison (no missing data)")

set.seed(666)

n <- 20
p <- 100

LF1 <- outer(rep(1, n), rep(1, p))
LL <- rnorm(n)
FF <- rnorm(p)
LF2 <- 5 * outer(LL, FF)
LF <- LF1 + LF2
M <- LF + 0.1 * rnorm(n * p)

test_that("matrix factor initialization is correct", {
  f <- flash(M, greedy.Kmax = 2, verbose = 0)
  expect_equal(f$n.factors, 2)
  expect_equal(fitted(f), LF1 + LF2, tol = 0.25, scale = 1)
})

test_that("nullcheck works as expected", {
  EF1 <- cbind(rep(1, n), c(LL[1:5], rep(0, n - 5)), 1:n)
  EF2 <- t(solve(crossprod(EF1), crossprod(EF1, M)))

  is.fixed <- matrix(TRUE, nrow = n, ncol = 3)
  is.fixed[6:n, 2] <- FALSE

  f <- flash.init(M) %>%
    flash.set.verbose(0) %>%
    flash.init.factors(list(EF1, EF2)) %>%
    flash.fix.loadings(kset = 1:3, mode = 1, is.fixed = is.fixed) %>%
    flash.backfit(verbose = 0)

  f <- flash.nullcheck(f, remove = TRUE)

  expect_equal(f$n.factors, 2)
})
