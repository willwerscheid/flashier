context("flashr comparison (no missing data)")

library(flashr)
set.seed(666)

n <- 20
p <- 100

LF1 <- outer(rep(1, n), rep(1, p))
LF2 <- 5 * outer(rnorm(n), rnorm(p))
LF <- LF1 + LF2
M <- LF + 0.1 * rnorm(n * p)

f <- init.flash(M)
f <- add.greedy(f)
f <- add.greedy(f)

test_that("matrix factor initialization is correct (using R)", {
  expect_equal(get.n.factors(f), 2)
  expect_equal(lowrank.expand(get.EF(f)), LF1 + LF2, tol = 0.25, scale = 1)
})

old.obj <- f$obj

test_that("the greedy objective approximately agrees with flashr (using R)", {
  flashr.greedy.res <- flashr::flash(M, var_type = "constant", Kmax = 2,
                                     nullcheck = FALSE)
  expect_equal(f$obj, flashr.greedy.res$objective, tol = 0.5, scale = 1)
})

f.b <- backfit.once(f, 1:2)

test_that ("the backfit objective agrees with flashr after one iteration (using R)", {
  expect_true(f.b$obj > old.obj)
  flashr.res <- flashr:::flash_backfit_workhorse(M, f_init = to.flashr(f),
                                                 var_type = "constant",
                                                 maxiter = 1, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash_set_data(M), flashr.res),
               f.b$obj)
})

f.b <- backfit(f, 1:2)

test_that ("the final backfit objective approximately agrees with flashr (using R)", {
  flashr.res <- flashr::flash(M, f_init = to.flashr(f), var_type = "constant",
                              greedy = FALSE, backfit = TRUE, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash_set_data(M), flashr.res),
               f.b$obj, tol = 0.01, scale = 1)
})

f <- init.flash(M, use.R = FALSE)
f <- add.greedy(f)
f <- add.greedy(f)

test_that("matrix factor initialization is correct (using Y)", {
  expect_equal(get.n.factors(f), 2)
  expect_equal(lowrank.expand(get.EF(f)), LF1 + LF2, tol = 0.25, scale = 1)
})

old.obj <- f$obj

test_that("the greedy objective approximately agrees with flashr (using Y)", {
  flashr.greedy.res <- flashr::flash(M, var_type = "constant", Kmax = 2,
                                     nullcheck = FALSE)
  expect_equal(f$obj, flashr.greedy.res$objective, tol = 0.5, scale = 1)
})

f.b <- backfit.once(f, 1:2)

test_that ("the backfit objective agrees with flashr after one iteration (using Y)", {
  expect_true(f.b$obj > old.obj)
  flashr.res <- flashr:::flash_backfit_workhorse(M, f_init = to.flashr(f),
                                                 var_type = "constant",
                                                 maxiter = 1, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash_set_data(M), flashr.res),
               f.b$obj)
})

f.b <- backfit(f, 1:2)

test_that ("the final backfit objective approximately agrees with flashr (using Y)", {
  flashr.res <- flashr::flash(M, f_init = to.flashr(f), var_type = "constant",
                              greedy = FALSE, backfit = TRUE, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash_set_data(M),
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash_set_data(M), flashr.res),
               f.b$obj, tol = 0.01, scale = 1)
})
