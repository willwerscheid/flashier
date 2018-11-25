context("fitting tensors")

set.seed(666)

m <- 10
n <- 10
p <- 30

LF1 <- rep(1, m) %o% rep(1, n) %o% rep(1, p)
LF2 <- 5 * (rnorm(m) %o% rnorm(n) %o% rnorm(p))
LF <- LF1 + LF2
M <- LF + 0.1 * rnorm(m * n * p)

f <- flashier(M, greedy.Kmax = 2, use.R = TRUE)

test_that("matrix factor initialization is correct (using R)", {
  expect_equal(get.n.factors(f), 2)
  expect_equal(lowrank.expand(get.EF(f)), LF1 + LF2, tol = 0.25, scale = 1)
})

f.b <- flashier(M, flash.init = f, backfit = "only", backfit.maxiter = 1,
                final.nullchk = FALSE)

test_that ("the backfit objective improves after one iteration (using R)", {
  expect_true(f.b$obj > f$obj)
})

f.b2 <- flashier(M, flash.init = f, backfit = "only")

test_that ("the final backfit objective improves again (using R)", {
  expect_true(f.b2$obj > f.b$obj)
})

missing <- sample(1:length(M), floor(0.1 * length(M)))
M.missing <- M
M.missing[missing] <- NA

f <- flashier(M.missing, greedy.Kmax = 2, use.R = FALSE)

test_that("matrix factor initialization is correct (using Y, with missing)", {
  expect_equal(get.n.factors(f), 2)
  expect_equal(lowrank.expand(get.EF(f)), LF1 + LF2, tol = 0.25, scale = 1)
})

f.b <- flashier(M, flash.init = f, backfit = "only", backfit.maxiter = 1,
                final.nullchk = FALSE)

test_that ("the backfit objective improves after one iteration (using Y, with missing)", {
  expect_true(f.b$obj > f$obj)
})

f.b2 <- flashier(M, flash.init = f, backfit = "only")

test_that ("the final backfit objective improves again (using Y, with missing)", {
  expect_true(f.b2$obj > f.b$obj)
})
