context("wrapup posterior sd")

# Regression test for L_psd / F_psd in the 2-D path:
# the fields must return posterior SD, not posterior variance.

test_that("L_psd / F_psd return posterior SD in 2-D path", {
  set.seed(1)
  n <- 200
  p <- 50
  Y <- matrix(rnorm(n * p), n, p) + tcrossprod(rnorm(n), rnorm(p))
  fit <- flash(Y, greedy_Kmax = 1, verbose = 0)

  samp <- fit$sampler(500)
  sd_emp_L <- apply(simplify2array(lapply(samp, `[[`, 1L)), c(1, 2), sd)
  sd_emp_F <- apply(simplify2array(lapply(samp, `[[`, 2L)), c(1, 2), sd)

  expect_equal(median(as.numeric(fit$L_psd / sd_emp_L)), 1, tolerance = 0.10)
  expect_equal(median(as.numeric(fit$F_psd / sd_emp_F)), 1, tolerance = 0.10)
})
