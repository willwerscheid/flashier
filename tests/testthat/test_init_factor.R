context("init_factor")

set.seed(666)

M <- matrix(5, nrow = 2, ncol = 3) + 0.1*rnorm(6)
M.svd <- svd(M)
correct.result <- M.svd$d[1] * outer(M.svd$u[, 1], M.svd$v[, 1])

test_that("matrix factor initialization is correct for R (no missing data)", {
  f <- init.flash(M)
  f.init <- init.next.EF(f)
  expect_equal(outer(f.init[[1]], f.init[[2]]), correct.result, tol = 1e-3)
})

Z <- matrix(1, nrow = nrow(M), ncol = ncol(M))

test_that("matrix factor initialization is correct for R with missing data", {
  f <- init.flash(M, nonmissing = Z)
  f.init <- init.next.EF(f)
  expect_equal(outer(f.init[[1]], f.init[[2]]), correct.result, tol = 1e-3)
})

EF.init <- list(matrix(rnorm(nrow(M)), ncol = 1),
               matrix(rnorm(ncol(M)), ncol = 1))
class(EF.init) <- "lowrank"
Y <- M + lowrank.expand(EF.init)

test_that("matrix factor initialization is correct for Y (no missing data)", {
  f <- init.flash(Y, EF.init = EF.init, use.R = FALSE)
  f.init <- init.next.EF(f)
  expect_equal(outer(f.init[[1]], f.init[[2]]), correct.result, tol = 1e-3)
})

test_that("matrix factor initialization is correct for Y with missing data", {
  f <- init.flash(Y, nonmissing = Z, EF.init = EF.init, use.R = FALSE)
  f.init <- init.next.EF(f)
  expect_equal(outer(f.init[[1]], f.init[[2]]), correct.result, tol = 1e-3)
})

A <- array(5, dim = c(4, 3, 2)) + 0.001 * rnorm(24)

test_that("tensor factor initialization is correct for R (no missing data)", {
  g <- init.flash(A)
  g.init <- init.next.EF(g)
  expect_equal(g.init[[1]] %o% g.init[[2]] %o% g.init[[3]], A, tol = 1e-2)
})

Z <- array(1, dim = dim(A))

test_that("tensor factor initialization is correct for R with missing data", {
  g <- init.flash(A, nonmissing = Z)
  g.init <- init.next.EF(g)
  expect_equal(g.init[[1]] %o% g.init[[2]] %o% g.init[[3]], A, tol = 1e-2)
})

EF.init <- list(matrix(rnorm(dim(A)[1]), ncol = 1),
               matrix(rnorm(dim(A)[2]), ncol = 1),
               matrix(rnorm(dim(A)[3]), ncol = 1))
class(EF.init) <- "lowrank"
Y <- A + lowrank.expand(EF.init)

test_that("tensor factor initialization is correct for Y (no missing data)", {
  g <- init.flash(Y, EF.init = EF.init, use.R = FALSE)
  g.init <- init.next.EF(g)
  expect_equal(g.init[[1]] %o% g.init[[2]] %o% g.init[[3]], A, tol = 1e-2)
})

test_that("tensor factor initialization is correct for Y with missing data", {
  g <- init.flash(Y, nonmissing = Z, EF.init = EF.init, use.R = FALSE)
  g.init <- init.next.EF(g)
  expect_equal(g.init[[1]] %o% g.init[[2]] %o% g.init[[3]], A, tol = 1e-2)
})

M <- outer(1:4, 1:6) + 0.1 * rnorm(24)

test_that("fixed factor initialization works as expected", {
  f <- init.flash(M,
                  fix.dim = list(1),
                  fix.idx = list(1:4),
                  fix.vals = list(1:4))
  f.init <- init.next.EF(f)
  expect_equal(f.init[[1]], 1:4)

  f <- init.flash(M,
                  fix.dim = list(2),
                  fix.idx = list(1:3),
                  fix.vals = list(1:3))
  f.init <- init.next.EF(f)
  expect_equal(f.init[[2]][1:3], 1:3)

  f <- init.flash(M,
                  fix.dim = list(1),
                  fix.idx = list(1),
                  fix.vals = list(0),
                  use.R = FALSE)
  f.init <- init.next.EF(f)
  expect_equal(f.init[[1]][1], 0)
  expect_false(any(f.init[[1]][2:4] == 0))
})

M <- outer(-1:3, -2:6)
M <- M + 0.1 * rnorm(length(M))

test_that("nonnegative factor initialization works as expected", {
  f <- init.flash(M, dim.signs = list(c(1, 0)))
  f.init <- init.next.EF(f)
  expect_true(all(f.init[[1]] >= 0))
  expect_true(any(f.init[[2]] < 0))

  f <- init.flash(M, dim.signs = list(c(1, -1)))
  f.init <- init.next.EF(f)
  expect_true(all(f.init[[1]] >= 0))
  expect_true(all(f.init[[2]] <= 0))

  f <- init.flash(M, dim.signs = list(c(0, 1)),
                  fix.dim = list(2),
                  fix.idx = list(4:5),
                  fix.vals = list(rep(0, 2)))
  f.init <- init.next.EF(f)
  expect_true(all(f.init[[2]] >= 0))
  expect_true(any(f.init[[1]] < 0))
  expect_true(all(f.init[[2]][4:5] == 0))
})
