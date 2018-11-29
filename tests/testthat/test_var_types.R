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
  flashr.res <- flashr:::flash_update_precision(flash.M, to.flashr(f$fit), "by_row")
  expect_equal(f$fit$tau, flashr.res$tau[, 1])
})

test_that("estimating variance by column produces identical estimates to flashr", {
  f <- flashier(M, greedy.Kmax = 1, var.type = 2)
  flashr.res <- flashr:::flash_update_precision(flash.M, to.flashr(f$fit), "by_column")
  expect_equal(f$fit$tau, flashr.res$tau[1, ])
})

test_that("zero variance type (with S constant) produces same fit as flashr", {
  f <- flashier(M, S = 0.1, var.type = NULL)
  expect_equal(f$fit$tau, f$fit$given.tau)

  flashr.res <- flashr::flash(flashr::flash_set_data(M, S = 0.1), Kmax = 1,
                              var_type = "zero", nullcheck = FALSE)
  expect_equal(f$obj, flashr.res$objective)
  expect_true(max(abs(flashr.res$fitted_values - lowrank.expand(f$fit$EF))) < 1e-6)
})

test_that("zero variance type (with S low-rank) produces same fit as flashr", {
  S <- 0.1 + 0.01 * rnorm(n)
  data <- set.flash.data(M, S, S.dim = 1)
  f <- flashier(data, greedy.Kmax = 1, var.type = NULL)
  expect_equal(f$fit$tau, f$fit$given.tau)

  flash.S <- matrix(S, nrow = n, ncol = p)
  flashr.res <- flashr::flash(flashr::flash_set_data(M, S = flash.S), Kmax = 1,
                              var_type = "zero", nullcheck = FALSE)
  expect_equal(f$obj, flashr.res$objective)
  expect_true(max(abs(flashr.res$fitted_values - lowrank.expand(f$fit$EF))) < 1e-6)
})

test_that("zero variance type (with S a matrix) produces same fit as flashr", {
  S <- matrix(0.1 + 0.01 * rnorm(n * p), nrow = n, ncol = p)
  f <- flashier(M, S = S, var.type = NULL)
  expect_equal(f$fit$tau, f$fit$given.tau)

  flashr.res <- flashr::flash(flashr::flash_set_data(M, S = S), Kmax = 1,
                              var_type = "zero", nullcheck = FALSE)
  expect_equal(f$obj, flashr.res$objective)
  expect_true(max(abs(flashr.res$fitted_values - lowrank.expand(f$fit$EF))) < 1e-6)
})

test_that("constant S + constant estimation works", {
  f <- flashier(M, S = 0.2, var.type = 0, greedy.Kmax = 1)
  expect_equal(f$fit$tau, f$fit$given.tau)

  f <- flashier(M, S = 0.05, var.type = 0, greedy.Kmax = 1)
  expect_equal(f$fit$tau, f$fit$est.tau)
})

test_that("by column S + by column estimation works", {
  tau = c(rep(50, 10), rep(250, p - 10))
  data <- set.flash.data(M, S = 1 / sqrt(tau), S.dim = 2)
  f <- flashier(data, var.type = 2, greedy.Kmax = 1)
  expect_equal(f$fit$tau[1:10], rep(50, 10))
  expect_equal(f$fit$tau[-(1:10)], f$fit$est.tau[-(1:10)])
})

test_that("kroncker variance estimation works", {
  Y <- matrix(10, nrow = 100, ncol = 100) + 0.1 * rnorm(100 * 100)
  f <- flashier(Y, var.type = c(1, 2), greedy.Kmax = 1)
  tau.mat <- r1.expand(f$fit$tau)
  expect_equal(mean(tau.mat), 100, tol = 0.1)

  R2 <- (Y - lowrank.expand(f$fit$EF))^2
  R2 <- R2 + lowrank.expand(f$fit$EF2) - lowrank.expand(lowrank.square(f$fit$EF))
  neg.llik <- function(x) {
    tau <- outer(x[1:100], x[101:200])
    return(-sum(log(tau)) + sum(R2 * tau))
  }
  optim.soln <- optim(rep(1, 200), neg.llik, method = "L-BFGS-B", lower = 0)
  optim.tau <- outer(optim.soln$par[1:100], optim.soln$par[101:200])
  expect_equal(tau.mat, optim.tau, tol = 0.1, scale = 1)
})

test_that("basic noisy variance estimation works", {
  f.const <- flashier(M, var.type = 0, greedy.Kmax = 1)
  f.noisy <- flashier(M, S = matrix(0.01, nrow = nrow(M), ncol = ncol(M)),
                      var.type = 0, greedy.Kmax = 1)
  expect_equal(f.const$fit$tau, f.noisy$fit$tau[1, 1], tol = 0.5, scale = 1)
  expect_equal(f.const$objective, f.noisy$objective, tol = 0.01, scale = 1)
})

test_that("fixed + by_column estimation works", {
  f.bycol <- flashier(M, var.type = 2, greedy.Kmax = 1)
  f.noisy <- flashier(M,
                      S = (matrix(0.01, nrow = nrow(M), ncol = ncol(M))
                           + 0.001 * rnorm(length(M))),
                      var.type = 2, greedy.Kmax = 1)
  expect_equal(f.bycol$fit$tau, f.noisy$fit$tau[1, ], tol = 0.5, scale = 1)
  expect_equal(f.bycol$objective, f.noisy$objective, tol = 0.1, scale = 1)
})

test_that("fixed + kronecker estimation works", {
  f.kron <- flashier(M, var.type = c(1, 2), greedy.Kmax = 0)
  f.noisy <- flashier(M, S = matrix(0.01, nrow = nrow(M), ncol = ncol(M)),
                      var.type = c(1, 2), greedy.Kmax = 0)

  expect_equal(r1.expand(f.kron$fit$tau), f.noisy$fit$tau, tol = 0.01, scale = 1)
  expect_equal(f.kron$objective, f.noisy$objective, tol = 0.01, scale = 1)

  f.kron <- flashier(M, var.type = c(1, 2), greedy.Kmax = 1)
  f.noisy <- flashier(M, S = matrix(0.01, nrow = nrow(M), ncol = ncol(M)),
                      var.type = c(1, 2), greedy.Kmax = 1)

  expect_equal(r1.expand(f.kron$fit$tau), f.noisy$fit$tau, tol = 1, scale = 1)
  expect_equal(f.kron$objective, f.noisy$objective, tol = 0.05, scale = 1)
})
