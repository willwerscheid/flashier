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

f <- flashier(M, greedy.Kmax = 2, use.R = TRUE, verbose.lvl = 0)

test_that("matrix factor initialization is correct (using R)", {
  expect_equal(f$n.factors, 2)
  expect_equal(lowrank.expand(get.EF(f$flash.fit)), LF1 + LF2, tol = 0.25, scale = 1)
})

f <- flashier(M, greedy.Kmax = 2, use.R = FALSE, verbose.lvl = 0)

test_that("matrix factor initialization is correct (using Y)", {
  expect_equal(f$n.factors, 2)
  expect_equal(lowrank.expand(get.EF(f$flash.fit)), LF1 + LF2, tol = 0.25, scale = 1)
})

f <- flashier(M, fix.dim = list(1, 1), fix.idx = list(1:n, 1:5),
              fix.vals = list(rep(1, n), LL[1:5]), backfit = "only",
              verbose.lvl = 0)

# Uncomment this test once ebnm can handle SEs equal to zero.
#
# f2 <- flashier(M, fix.dim = list(1, 1), fix.idx = list(1:n, 1:5),
#                fix.vals = list(rep(1, n), LL[1:5]),
#                backfit = "only",
#                use.fixed.to.est.g = TRUE,
#                verbose.lvl = 0)
#
# test_that("results are different if fixed elements are included in priors", {
#   expect_false(f$elbo == f2$elbo)
# })

f <- flashier(M, fix.dim = list(1, 1, 1), fix.idx = list(1:n, 1:5, 1:n),
              fix.vals = list(rep(1, n), LL[1:5], 1:n), backfit = "only",
              nullchk.fixed = TRUE, verbose.lvl = 0)

test_that("nullcheck works as expected", {
  expect_equal(is.zero(f$flash.fit), c(FALSE, FALSE, TRUE))
})
