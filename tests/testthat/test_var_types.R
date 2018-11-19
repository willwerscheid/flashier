context("flashr comparison (variance types)")

library(flashr)
set.seed(666)

n <- 40
p <- 60

LF <- outer(rep(1, n), rep(1, p))
M <- LF + 0.1 * rnorm(n * p)
flash.M <- flash_set_data(M)

test_that("estimating variance by row produces identical estimates to flashr", {
  f <- init.flash(M, est.tau.dim = 1)
  f <- add.next.factor(f)

  flashr.res <- flashr:::flash_update_precision(flash.M, to.flashr(f), "by_row")
  expect_equal(f$tau, flashr.res$tau[, 1])
})

test_that("estimating variance by column produces identical estimates to flashr", {
  f <- init.flash(M, est.tau.dim = 2)
  f <- add.next.factor(f)

  flashr.res <- flashr:::flash_update_precision(flash.M, to.flashr(f), "by_column")
  expect_equal(f$tau, flashr.res$tau[1, ])
})

test_that("zero variance type (with S constant) produces same fit as flashr", {
  f <- init.flash(M, est.tau.dim = NULL, given.tau = 100, given.tau.dim = 0)
  f <- add.next.factor(f)
  expect_equal(f$tau, f$given.tau)

  flashr.res <- flashr::flash(flashr::flash_set_data(M, S = 0.1), Kmax = 1,
                              var_type = "zero", nullcheck = FALSE)
  expect_equal(f$obj, flashr.res$objective)
  expect_true(max(abs(flashr.res$fitted_values - lowrank.expand(f$EF))) < 1e-6)
})

test_that("zero variance type (with S low-rank) produces same fit as flashr", {
  S <- 0.1 + 0.01 * rnorm(n)
  f <- init.flash(M, est.tau.dim = NULL, given.tau = 1 / S^2, given.tau.dim = 1)
  f <- add.next.factor(f)
  expect_equal(f$tau, f$given.tau)

  flash.S <- matrix(S, nrow = n, ncol = p)
  flashr.res <- flashr::flash(flashr::flash_set_data(M, S = flash.S), Kmax = 1,
                              var_type = "zero", nullcheck = FALSE)
  expect_equal(f$obj, flashr.res$objective)
  expect_true(max(abs(flashr.res$fitted_values - lowrank.expand(f$EF))) < 1e-6)
})

test_that("zero variance type (with S a matrix) produces same fit as flashr", {
  S <- matrix(0.1 + 0.01 * rnorm(n * p), nrow = n, ncol = p)
  f <- init.flash(M, est.tau.dim = NULL, given.tau = 1 / S^2)
  f <- add.next.factor(f)
  expect_equal(f$tau, f$given.tau)

  flashr.res <- flashr::flash(flashr::flash_set_data(M, S = S), Kmax = 1,
                              var_type = "zero", nullcheck = FALSE)
  expect_equal(f$obj, flashr.res$objective)
  expect_true(max(abs(flashr.res$fitted_values - lowrank.expand(f$EF))) < 1e-6)
})
