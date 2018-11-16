context("init_factor")

set.seed(666)

M <- matrix(5, nrow = 2, ncol = 3) + 0.1*rnorm(6)
M.svd <- svd(M)
correct.result <- M.svd$d[1] * outer(M.svd$u[, 1], M.svd$v[, 1])

test_that("matrix factor initialization is correct for R (no missing data)", {
  f <- init.flash(M)
  f.init <- init.r1(f)
  expect_equal(outer(f.init[[1]], f.init[[2]]), correct.result, tol = 1e-3)
})

Z <- matrix(1, nrow = nrow(M), ncol = ncol(M))

test_that("matrix factor initialization is correct for R with missing data", {
  f <- init.flash(M, nonmissing = Z)
  f.init <- init.r1(f)
  expect_equal(outer(f.init[[1]], f.init[[2]]), correct.result, tol = 1e-3)
})

F.init <- list(matrix(rnorm(nrow(M)), ncol = 1),
               matrix(rnorm(ncol(M)), ncol = 1))
class(F.init) <- "lowrank"
Y <- M + lowrank.expand(F.init)

test_that("matrix factor initialization is correct for Y (no missing data)", {
  f <- init.flash(Y, F.init = F.init, use.R = FALSE)
  f.init <- init.r1(f)
  expect_equal(outer(f.init[[1]], f.init[[2]]), correct.result, tol = 1e-3)
})

test_that("matrix factor initialization is correct for Y with missing data", {
  f <- init.flash(Y, nonmissing = Z, F.init = F.init, use.R = FALSE)
  f.init <- init.r1(f)
  expect_equal(outer(f.init[[1]], f.init[[2]]), correct.result, tol = 1e-3)
})

A <- array(5, dim = c(4, 3, 2)) + 0.001 * rnorm(24)

test_that("tensor factor initialization is correct for R (no missing data)", {
  g <- init.flash(A)
  g.init <- init.r1(g)
  expect_equal(g.init[[1]] %o% g.init[[2]] %o% g.init[[3]], A, tol = 1e-2)
})

Z <- array(1, dim = dim(A))

test_that("tensor factor initialization is correct for R with missing data", {
  g <- init.flash(A, nonmissing = Z)
  g.init <- init.r1(g)
  expect_equal(g.init[[1]] %o% g.init[[2]] %o% g.init[[3]], A, tol = 1e-2)
})

F.init <- list(matrix(rnorm(dim(A)[1]), ncol = 1),
               matrix(rnorm(dim(A)[2]), ncol = 1),
               matrix(rnorm(dim(A)[3]), ncol = 1))
class(F.init) <- "lowrank"
Y <- A + lowrank.expand(F.init)

test_that("tensor factor initialization is correct for Y (no missing data)", {
  g <- init.flash(Y, F.init = F.init, use.R = FALSE)
  g.init <- init.r1(g)
  expect_equal(g.init[[1]] %o% g.init[[2]] %o% g.init[[3]], A, tol = 1e-2)
})

test_that("tensor factor initialization is correct for Y with missing data", {
  g <- init.flash(Y, nonmissing = Z, F.init = F.init, use.R = FALSE)
  g.init <- init.r1(g)
  expect_equal(g.init[[1]] %o% g.init[[2]] %o% g.init[[3]], A, tol = 1e-2)
})
