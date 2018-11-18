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
Z <- matrix(1, nrow = nrow(M), ncol = ncol(M))

missing_idx <- sample(1:length(M), floor(0.6 * length(M)))
M[missing_idx] <- 0
Z[missing_idx] <- 0

flash.data <- M
flash.data[missing_idx] <- NA
flash.data <- flash_set_data(flash.data)

f <- init.flash(M, nonmissing = Z)
f <- add.next.factor(f)
f <- add.next.factor(f)

test_that("matrix factor initialization is correct (using R)", {
  expect_equal(get.n.factors(f), 2)
  expect_equal(lowrank.expand(get.EF(f)), LF1 + LF2, tol = 0.25, scale = 1)
})

old.obj <- f$obj

test_that("the greedy objective approximately agrees with flashr (using R)", {
  flashr.greedy.res <- flashr::flash(flash.data, var_type = "constant",
                                     Kmax = 2, nullcheck = FALSE)
  expect_equal(f$obj, flashr.greedy.res$objective, tol = 0.5, scale = 1)
})

f.b <- backfit.once(f, 1:2)

test_that ("the backfit objective agrees with flashr after one iteration (using R)", {
  expect_true(f.b$obj > old.obj)
  flashr.res <- flashr:::flash_backfit_workhorse(flash.data,
                                                 f_init = to.flashr(f),
                                                 var_type = "constant",
                                                 maxiter = 1, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash.data,
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash.data, flashr.res),
               f.b$obj)
})

f.b <- backfit(f, 1:2)

test_that ("the final backfit objective approximately agrees with flashr (using R)", {
  flashr.res <- flashr::flash(flash.data, var_type = "constant", Kmax = 2,
                              greedy = TRUE, backfit = TRUE, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash.data,
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash.data, flashr.res),
               f.b$obj, tol = 0.01, scale = 1)
})

f.sprs <- init.flash(Matrix(M), nonmissing = Matrix(Z), use.R = FALSE)
f.sprs <- add.next.factor(f.sprs)
f.sprs <- add.next.factor(f.sprs)

test_that("matrix factor initialization is correct (sparse, using Y)", {
  expect_equal(get.n.factors(f.sprs), 2)
  expect_equal(lowrank.expand(get.EF(f.sprs)), LF1 + LF2, tol = 0.25, scale = 1)
})

old.obj <- f.sprs$obj

test_that("the greedy objective approximately agrees with flashr (sparse, using Y)", {
  flashr.greedy.res <- flashr::flash(flash.data, var_type = "constant",
                                     Kmax = 2, nullcheck = FALSE)
  expect_equal(f.sprs$obj, flashr.greedy.res$objective, tol = 0.5, scale = 1)
})

f.sprs.b <- backfit.once(f.sprs, 1:2)

test_that ("the backfit objective agrees with flashr after one iteration (sparse, using Y)", {
  expect_true(f.sprs.b$obj > old.obj)
  flashr.res <- flashr:::flash_backfit_workhorse(flash.data,
                                                 f_init = to.flashr(f.sprs),
                                                 var_type = "constant",
                                                 maxiter = 1, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash.data,
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash.data, flashr.res),
               f.sprs.b$obj)
})

f.sprs.b <- backfit(f.sprs, 1:2)

test_that ("the final backfit objective approximately agrees with flashr (sparse, using Y)", {
  flashr.res <- flashr::flash(flash.data, var_type = "constant", Kmax = 2,
                              greedy = TRUE, backfit = TRUE, nullcheck = FALSE)
  flashr.res <- flashr:::flash_update_precision(flash.data,
                                                flashr.res$fit,
                                                var_type = "constant")
  expect_equal(flashr:::flash_get_objective(flash.data, flashr.res),
               f.b$obj, tol = 0.01, scale = 1)
})
