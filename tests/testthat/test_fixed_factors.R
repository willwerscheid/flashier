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

fl <- flash.init(M) %>%
  flash.set.verbose(0) %>%
  flash.init.factors(list(EF1, EF2)) %>%
  flash.fix.loadings(1:2, mode = 1) %>%
  flash.backfit()

test_that("Fixed factors are correctly added to a new flash object", {
  expect_equal(fl$flash.fit$EF[[1]][, 1], rep(1, n))
  expect_equal(fl$flash.fit$EF[[1]][, 2], 1:n)
})

fl <- flash.init(M) %>%
  flash.set.verbose(0) %>%
  flash.add.greedy(1) %>%
  flash.init.factors(list(EF1, EF2)) %>%
  flash.fix.loadings(2:3, mode = 1) %>%
  flash.backfit()

test_that("Fixed factors are correctly added to an existing flash object", {
  expect_equal(fl$flash.fit$EF[[1]][, 2], rep(1, n))
  expect_equal(fl$flash.fit$EF[[1]][, 3], 1:n)
})
