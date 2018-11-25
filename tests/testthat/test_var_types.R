context("variance types")

library(flashr)
set.seed(666)

n <- 40
p <- 60

LF <- outer(rep(1, n), rep(1, p))
M <- LF + 0.1 * rnorm(n * p)
flash.M <- flash_set_data(M)

test_that("estimating variance by row produces identical estimates to flashr", {
  f <- flashier(M, greedy.Kmax = 1, var.type = 1)
  flashr.res <- flashr:::flash_update_precision(flash.M, to.flashr(f), "by_row")
  expect_equal(f$tau, flashr.res$tau[, 1])
})

test_that("estimating variance by column produces identical estimates to flashr", {
  f <- flashier(M, greedy.Kmax = 1, var.type = 2)
  flashr.res <- flashr:::flash_update_precision(flash.M, to.flashr(f), "by_column")
  expect_equal(f$tau, flashr.res$tau[1, ])
})

test_that("zero variance type (with S constant) produces same fit as flashr", {
  f <- flashier(M, S = 0.1, var.type = NULL)
  expect_equal(f$tau, f$given.tau)

  flashr.res <- flashr::flash(flashr::flash_set_data(M, S = 0.1), Kmax = 1,
                              var_type = "zero", nullcheck = FALSE)
  expect_equal(f$obj, flashr.res$objective)
  expect_true(max(abs(flashr.res$fitted_values - lowrank.expand(f$EF))) < 1e-6)
})

test_that("zero variance type (with S low-rank) produces same fit as flashr", {
  S <- 0.1 + 0.01 * rnorm(n)
  data <- set.flash.data(M, S, S.dim = 1)
  f <- flashier(data, greedy.Kmax = 1, var.type = NULL)
  expect_equal(f$tau, f$given.tau)

  flash.S <- matrix(S, nrow = n, ncol = p)
  flashr.res <- flashr::flash(flashr::flash_set_data(M, S = flash.S), Kmax = 1,
                              var_type = "zero", nullcheck = FALSE)
  expect_equal(f$obj, flashr.res$objective)
  expect_true(max(abs(flashr.res$fitted_values - lowrank.expand(f$EF))) < 1e-6)
})

test_that("zero variance type (with S a matrix) produces same fit as flashr", {
  S <- matrix(0.1 + 0.01 * rnorm(n * p), nrow = n, ncol = p)
  f <- flashier(M, S = S, var.type = NULL)
  expect_equal(f$tau, f$given.tau)

  flashr.res <- flashr::flash(flashr::flash_set_data(M, S = S), Kmax = 1,
                              var_type = "zero", nullcheck = FALSE)
  expect_equal(f$obj, flashr.res$objective)
  expect_true(max(abs(flashr.res$fitted_values - lowrank.expand(f$EF))) < 1e-6)
})

test_that("constant S + constant estimation works", {
  f <- flashier(M, S = 0.2, var.type = 0, greedy.Kmax = 1)
  expect_equal(f$tau, f$given.tau)

  f <- flashier(M, S = 0.05, var.type = 0, greedy.Kmax = 1)
  expect_equal(f$tau, f$est.tau)
})

test_that("by column S + by column estimation works", {
  tau = c(rep(50, 10), rep(250, p - 10))
  data <- set.flash.data(M, S = 1 / sqrt(tau), S.dim = 2)
  f <- flashier(data, var.type = 2, greedy.Kmax = 1)
  expect_equal(f$tau[1:10], rep(50, 10))
  expect_equal(f$tau[-(1:10)], f$est.tau[-(1:10)])
})
