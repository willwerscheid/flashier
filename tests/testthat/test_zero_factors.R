context("zero factors")

set.seed(666)

n <- 20
p <- 100

LF1 <- outer(rep(1, n), rep(1, p))
LL <- rnorm(n)
FF <- rnorm(p)
LF2 <- 5 * outer(LL, FF)
LF <- LF1 + LF2
M <- LF + 0.1 * rnorm(n * p)

# Factor 3 will be removed here:
f <- flashier(M, greedy.Kmax = 2,
              fixed.factors = fixed.factors(1, 1:n),
              backfit = "final", backfit.reltol = 100,
              nullchk.fixed = TRUE, verbose.lvl = 0)

test_that("sampler works correctly when there are zero factors", {
  samp <- f$sampler(1)[[1]]
  expect_equal(samp[[1]][, 1], rep(0, n))
  expect_equal(samp[[2]][, 1], rep(0, p))
})

test_that("zero factors are dealt with correctly when normalizing loadings", {
  expect_equal(f$loadings$scale.constant[1], 0)
  expect_equal(f$loadings$normalized.loadings[[1]][, 1], rep(0, n))
  expect_equal(f$loadings$normalized.loadings[[2]][, 1], rep(0, p))
})
