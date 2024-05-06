context("Fixed factors interface")

set.seed(666)

n <- 20
p <- 100

LF1 <- outer(rep(1, n), rep(1, p))
LL <- rnorm(n)
FF <- rnorm(p)
LF2 <- 5 * outer(LL, FF)
LF <- LF1 + LF2
M <- LF + 0.1 * rnorm(n * p)

EF1 <- cbind(1, 1:n)
EF2 <- t(solve(crossprod(EF1), crossprod(EF1, M)))

fl <- flash_init(M) |>
  flash_set_verbose(0) |>
  flash_factors_init(list(EF1, EF2)) |>
  flash_factors_fix(1:2, which_dim = "loadings") |>
  flash_backfit()

test_that("Fixed factors are correctly added to a new flash object", {
  expect_equal(fl$flash_fit$EF[[1]][, 1], rep(1, n))
  expect_equal(fl$flash_fit$EF[[1]][, 2], 1:n)
})

fl <- flash_init(M) |>
  flash_set_verbose(0) |>
  flash_greedy(1) |>
  flash_factors_init(list(EF1, EF2)) |>
  flash_factors_fix(2:3, which_dim = "loadings") |>
  flash_backfit()

test_that("Fixed factors are correctly added to an existing flash object", {
  expect_equal(fl$flash_fit$EF[[1]][, 2], rep(1, n))
  expect_equal(fl$flash_fit$EF[[1]][, 3], 1:n)
})
