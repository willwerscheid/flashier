context("Fixed factors interface")

set.seed(666)

n <- 20
p <- 100

LF1 <- outer(rep(1, n), rep(1, p))
LL <- rnorm(n)
FF <- rnorm(p)
LF2 <- 5 * outer(LL, FF)
LF <- LF1 + LF2
M <- LF + 0.1 * rnorm(n * p)

fl <- flashier(M, fixed.factors = c(fixed.ones(n = 1),
                                    fixed.factors(n = 1, vals = 1:n)),
               greedy.Kmax = 1, verbose.lvl = 0)

test_that("Fixed factors are correctly added to a new flash object", {
  expect_equal(fl$flash.fit$EF[[1]][, 1], rep(1, n))
  expect_equal(fl$flash.fit$EF[[1]][, 2], 1:n)
})

fl <- flashier(init = fl, fixed.factors = fixed.sparse(n = 1, 1:3),
               fit = "backfit", verbose.lvl = 0)

test_that("Fixed factors are correctly added to an existing flash object", {
  expect_equal(fl$flash.fit$EF[[1]][4:n, 4], rep(0, n - 3))
})
