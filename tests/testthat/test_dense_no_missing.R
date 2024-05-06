context("flashr comparison (no missing data)")

set.seed(666)

n <- 20
p <- 100

LF1 <- outer(rep(1, n), rep(1, p))
LL <- rnorm(n)
FF <- rnorm(p)
LF2 <- 5 * outer(LL, FF)
LF <- LF1 + LF2
M <- LF + 0.1 * rnorm(n * p)

test_that("matrix factor initialization is correct", {
  f <- flash(M, greedy_Kmax = 2, verbose = 0)
  expect_equal(f$n_factors, 2)
  expect_equal(fitted(f), LF1 + LF2, tol = 0.25, scale = 1)
})

test_that("nullcheck works as expected", {
  EF1 <- cbind(rep(1, n), c(LL[1:5], rep(0, n - 5)), 1:n)
  EF2 <- t(solve(crossprod(EF1), crossprod(EF1, M)))

  is.fixed <- matrix(TRUE, nrow = n, ncol = 3)
  is.fixed[6:n, 2] <- FALSE

  f <- flash_init(M) |>
    flash_set_verbose(0) |>
    flash_factors_init(list(EF1, EF2)) |>
    flash_factors_fix(kset = 1:3, which_dim = "loadings", fixed_idx = is.fixed) |>
    flash_backfit(verbose = 0)

  f <- flash_nullcheck(f, remove = TRUE)

  expect_equal(f$n_factors, 2)
})
