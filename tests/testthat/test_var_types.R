context("variance types")

set.seed(666)

n <- 40
p <- 60

LF <- outer(rep(1, n), rep(1, p))
M <- LF + 0.1 * rnorm(n * p)

test_that("constant S + constant estimation works", {
  f <- flash(M, S = 0.2, var_type = 0, greedy_Kmax = 1, verbose = 0)
  expect_equal(f$flash_fit$tau, f$flash_fit$given.tau)

  f <- flash(M, S = 0.05, var_type = 0, greedy_Kmax = 1, verbose = 0)
  expect_equal(f$flash_fit$tau, f$flash_fit$est.tau)
})

test_that("by column S + by column estimation works", {
  tau = c(rep(50, 10), rep(250, p - 10))
  f <- flash_init(M, S = 1 / sqrt(tau), S_dim = 2, var_type = 2) |>
    flash_greedy(1, verbose = 0)
  expect_equal(f$flash_fit$tau[1:10], rep(50, 10))
  expect_equal(f$flash_fit$tau[-(1:10)], f$flash_fit$est.tau[-(1:10)])
})

test_that("kronecker variance estimation works", {
  Y <- matrix(10, nrow = 100, ncol = 100) + 0.1 * rnorm(100 * 100)
  f <- flash(Y, var_type = c(1, 2), greedy_Kmax = 1, verbose = 0)
  tau.mat <- r1.expand(f$flash_fit$tau)
  expect_equal(mean(tau.mat), 100, tol = 0.1)

  R2 <- (Y - lowrank.expand(f$flash_fit$EF))^2
  R2 <- R2 + lowrank.expand(f$flash_fit$EF2) - lowrank.expand(lowrank.square(f$flash_fit$EF))
  neg.llik <- function(x) {
    tau <- outer(x[1:100], x[101:200])
    return(-sum(log(tau)) + sum(R2 * tau))
  }
  optim.soln <- optim(rep(1, 200), neg.llik, method = "L-BFGS-B", lower = 0)
  optim.tau <- outer(optim.soln$par[1:100], optim.soln$par[101:200])
  expect_equal(tau.mat, optim.tau, tol = 0.1, scale = 1)
})

test_that("basic noisy variance estimation works", {
  f.const <- flash(M, var_type = 0, greedy_Kmax = 1, verbose = 0)
  f.noisy <- flash(M, S = matrix(0.01, nrow = nrow(M), ncol = ncol(M)),
                   var_type = 0, greedy_Kmax = 1, verbose = 0)
  expect_equal(f.const$flash_fit$tau, f.noisy$flash_fit$tau[1, 1], tol = 0.5, scale = 1)
  expect_equal(f.const$elbo, f.noisy$elbo, tol = 0.01, scale = 1)
})

test_that("fixed + by_column estimation works", {
  f.bycol <- flash(M, var_type = 2, greedy_Kmax = 1, verbose = 0)
  f.noisy <- flash(M, S = (matrix(0.01, nrow = nrow(M), ncol = ncol(M))
                           + 0.001 * rnorm(length(M))),
                   var_type = 2, greedy_Kmax = 1, verbose = 0)
  expect_equal(f.bycol$flash_fit$tau, f.noisy$flash_fit$tau[1, ], tol = 0.5, scale = 1)
  expect_equal(f.bycol$elbo, f.noisy$elbo, tol = 0.1, scale = 1)
})

test_that("fixed + kronecker estimation works", {
  f.kron <- flash_init(M, var_type = c(1, 2))
  f.noisy <- flash_init(M, S = matrix(0.01, nrow = nrow(M), ncol = ncol(M)),
                        var_type = c(1, 2))

  expect_equal(r1.expand(f.kron$flash_fit$tau), f.noisy$flash_fit$tau, tol = 0.01, scale = 1)
  expect_equal(f.kron$elbo, f.noisy$elbo, tol = 0.01, scale = 1)

  f.kron <- flash(M, var_type = c(1, 2), greedy_Kmax = 1, verbose = 0)
  f.noisy <- flash(M, S = matrix(0.01, nrow = nrow(M), ncol = ncol(M)),
                   var_type = c(1, 2), greedy_Kmax = 1, verbose = 0)

  expect_equal(r1.expand(f.kron$flash_fit$tau), f.noisy$flash_fit$tau, tol = 1, scale = 1)
  expect_equal(f.kron$elbo, f.noisy$elbo, tol = 0.05, scale = 1)
})
