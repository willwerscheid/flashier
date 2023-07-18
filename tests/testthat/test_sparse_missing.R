context("flashr comparison (with missing data, sparse)")

library(Matrix)
set.seed(666)

n <- 40
p <- 100

LF1 <- outer(rep(1, n), rep(1, p))
LF2 <- 5 * outer(rnorm(n), rnorm(p))
LF <- LF1 + LF2
M <- LF + 0.1 * rnorm(n * p)

missing <- sample(1:length(M), floor(0.6 * length(M)))
M[missing] <- NA

f <- flash(M, greedy_Kmax = 2, verbose = 0)

test_that("matrix factor initialization is correct (sparse)", {
  expect_equal(f$n_factors, 2)
  expect_equal(lowrank.expand(get.EF(f$flash_fit)), LF1 + LF2, tol = 0.25, scale = 1)
})
