context("flashr comparison (with missing data, sparse)")

library(flashr)
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

f <- flashier(M, greedy.Kmax = 2, use.R = TRUE, verbose.lvl = 0)

test_that("matrix factor initialization is correct (using R)", {
  expect_equal(f$n.factors, 2)
  expect_equal(lowrank.expand(get.EF(f$flash.fit)), LF1 + LF2, tol = 0.25, scale = 1)
})

old.obj <- f$elbo

test_that("the greedy objective approximately agrees with flashr (using R)", {
  flashr.greedy.res <- flashr::flash(M, var_type = "constant", Kmax = 2)
  expect_equal(f$elbo, flashr.greedy.res$objective, tol = 0.5, scale = 1)
})

f.b <- flashier(flash.init = f, backfit = "only", backfit.maxiter = 1,
                final.nullchk = FALSE, verbose.lvl = 0)

test_that ("the backfit objective agrees with flashr after one iteration (using R)", {
  expect_true(f.b$elbo > old.obj)
  flashr.res <- flashr:::flash_backfit_workhorse(M,
                                                 f_init = to.flashr(f$flash.fit),
                                                 var_type = "constant",
                                                 maxiter = 1, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(M, flashr.res),
               f.b$elbo)
})

f.b <- flashier(flash.init = f, backfit = "only", verbose.lvl = 0,
                greedy.tol = 0.01, backfit.reltol = 1)

test_that ("the final backfit objective approximately agrees with flashr (using R)", {
  flashr.res <- flashr::flash(M, var_type = "constant", Kmax = 2,
                              greedy = TRUE, backfit = TRUE, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(M, flashr.res),
               f.b$elbo, tol = 0.02, scale = 1)
})

f.sprs <- flashier(Matrix(M), greedy.Kmax = 2, use.R = FALSE, verbose.lvl = 0)

test_that("matrix factor initialization is correct (sparse, using Y)", {
  expect_equal(f.sprs$n.factors, 2)
  expect_equal(lowrank.expand(get.EF(f.sprs$flash.fit)), LF1 + LF2, tol = 0.25, scale = 1)
})

old.obj <- f.sprs$elbo

test_that("the greedy objective approximately agrees with flashr (sparse, using Y)", {
  flashr.greedy.res <- flashr::flash(M, var_type = "constant",
                                     Kmax = 2, nullcheck = FALSE)
  expect_equal(f.sprs$elbo, flashr.greedy.res$objective, tol = 0.5, scale = 1)
})

f.sprs.b <- flashier(flash.init = f.sprs, backfit = "only",
                     backfit.maxiter = 1, final.nullchk = FALSE,
                     verbose.lvl = 0)

test_that ("the backfit objective agrees with flashr after one iteration (sparse, using Y)", {
  expect_true(f.sprs.b$elbo > old.obj)
  flashr.res <- flashr:::flash_backfit_workhorse(M,
                                                 f_init = to.flashr(f.sprs$flash.fit),
                                                 var_type = "constant",
                                                 maxiter = 1, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(M, flashr.res),
               f.sprs.b$elbo)
})

f.sprs.b <- flashier(flash.init = f.sprs, backfit = "only", verbose.lvl = 0)

test_that ("the final backfit objective approximately agrees with flashr (sparse, using Y)", {
  flashr.res <- flashr::flash(M, var_type = "constant", Kmax = 2,
                              greedy = TRUE, backfit = TRUE, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(M, flashr.res),
               f.b$elbo, tol = 0.02, scale = 1)
})
