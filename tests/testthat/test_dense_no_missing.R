context("flashr comparison (no missing data)")

library(flashr)
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

old.obj <- f$elbo

test_that("the greedy objective approximately agrees with flashr (using R)", {
  flashr.greedy.res <- flashr::flash(M, var_type = "constant", Kmax = 2)
  expect_equal(f$elbo, flashr.greedy.res$objective, tol = 0.5, scale = 1)
})

f.b <- flashier(flash.init = f, backfit = "only", backfit.maxiter = 1,
                verbose.lvl = 0)

test_that ("the backfit objective agrees with flashr after one iteration (using R)", {
  expect_true(f.b$elbo > old.obj)
  flashr.res <- flashr:::flash_backfit_workhorse(M, f_init = to.flashr(f$flash.fit),
                                                 var_type = "constant",
                                                 maxiter = 1, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash_set_data(M), flashr.res),
               f.b$elbo, tol = 0.001)
})

f.b <- flashier(flash.init = f, backfit = "only", verbose.lvl = 0)

test_that ("the final backfit objective approximately agrees with flashr (using R)", {
  flashr.res <- flashr::flash(M, var_type = "constant", Kmax = 2, backfit = TRUE)
  expect_equal(flashr:::flash_get_objective(flash_set_data(M), flashr.res),
               f.b$elbo, tol = 0.1, scale = 1)
})

f <- flashier(M, greedy.Kmax = 2, use.R = FALSE, verbose.lvl = 0)

test_that("matrix factor initialization is correct (using Y)", {
  expect_equal(f$n.factors, 2)
  expect_equal(lowrank.expand(get.EF(f$flash.fit)), LF1 + LF2, tol = 0.25, scale = 1)
})

old.obj <- f$elbo

test_that("the greedy objective approximately agrees with flashr (using Y)", {
  flashr.greedy.res <- flashr::flash(M, var_type = "constant", Kmax = 2,
                                     nullcheck = FALSE)
  expect_equal(f$elbo, flashr.greedy.res$objective, tol = 0.5, scale = 1)
})

f.b <- flashier(flash.init = f, backfit = "only", backfit.maxiter = 1,
                verbose.lvl = 0)

test_that ("the backfit objective agrees with flashr after one iteration (using Y)", {
  expect_true(f.b$elbo > old.obj)
  flashr.res <- flashr:::flash_backfit_workhorse(M, f_init = to.flashr(f$flash.fit),
                                                 var_type = "constant",
                                                 maxiter = 1, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash_set_data(M), flashr.res),
               f.b$elbo, tol = 0.001)
})

f.b <- flashier(flash.init = f, backfit = "only", verbose.lvl = 0)

test_that ("the final backfit objective approximately agrees with flashr (using Y)", {
  flashr.res <- flashr::flash(M, f_init = to.flashr(f$flash.fit), var_type = "constant",
                              greedy = FALSE, backfit = TRUE, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash_set_data(M), flashr.res),
               f.b$elbo, tol = 0.01, scale = 1)
})

f <- flashier(M, fix.dim = list(1, 1), fix.idx = list(1:n, 1:5),
              fix.vals = list(rep(1, n), LL[1:5]), backfit = "only",
              verbose.lvl = 0)

test_that("the objective after adding fixed factors approximately agrees with flashr", {
  flashr.res <- flashr:::flash_backfit_workhorse(M, kset = 1:2,
                                                 f_init = to.flashr(f$flash.fit),
                                                 var_type = "constant",
                                                 maxiter = 100,
                                                 nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash_set_data(M), flashr.res),
               f$elbo, tol = 0.1, scale = 1)

  fixed_loadings <- matrix(1, nrow = n, ncol = 2)
  fixed_loadings[1:5, 2] <- LL[1:5]
  fixed_loadings[-(1:5), 2] <- 0
  fixl <- matrix(TRUE, nrow = n, ncol = 2)
  fixl[-(1:5), 2] <- FALSE
  flashr.res <- flashr::flash_add_fixed_loadings(M,
                                                 LL = fixed_loadings,
                                                 fixl = fixl,
                                                 var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash_set_data(M), flashr.res),
               f$elbo, tol = 0.1, scale = 1)
})

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
