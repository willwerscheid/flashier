context("ops")

M <- matrix(1:6, nrow = 2, ncol = 3)
A <- array(M, dim = c(2, 3, 1))
v1 <- 7:8
v2 <- 9:11
Y <- matrix(1:12, nrow = 3, ncol = 4)

test_that("r1 product is the same for general matrices and arrays", {
  expect_equal(nmode.prod.r1(M, v2, 1),
               nmode.prod.r1(A, list(v2, 1), 1))
  expect_equal(nmode.prod.r1(M, v1, 2),
               nmode.prod.r1(A, list(v1, 1), 2))
})

lowrank.M <- list(matrix(1:10, ncol = 2),
                  matrix(11:16, ncol = 2))
class(lowrank.M) <- "lowrank"
lowrank.A <- c(lowrank.M,
               list(matrix(1, nrow = 1, ncol = 2)))
class(lowrank.A) <- "lowrank"
v1 <- 1:5
v2 <- 6:8

test_that("r1 product is the same for low-rank matrices and arrays", {
  correct.result.1 <- as.vector(lowrank.M[[1]] %*% t(lowrank.M[[2]]) %*% v2)
  expect_equal(nmode.prod.r1(lowrank.M, list(v2), 1),
               correct.result.1)
  expect_equal(nmode.prod.r1(lowrank.A, list(v2, 1), 1),
               correct.result.1)

  correct.result.2 <- as.vector(v1 %*% lowrank.M[[1]] %*% t(lowrank.M[[2]]))
  expect_equal(nmode.prod.r1(lowrank.M, list(v1), 2),
               correct.result.2)
  expect_equal(nmode.prod.r1(lowrank.A, list(v1, 1), 2),
               correct.result.2)
})

Z <- matrix(rnorm(15), nrow = 5, ncol = 3)

test_that("r1 product with premultiplication is correct for matrices", {
  correct.result.1 <- (Z * (lowrank.M[[1]] %*% t(lowrank.M[[2]]))) %*% v2
  expect_equal(premult.lowrank.nmode.prod.r1(Z, lowrank.M, list(v2), 1),
               as.vector(correct.result.1))

  correct.result.2 <- v1 %*% (Z * (lowrank.M[[1]] %*% t(lowrank.M[[2]])))
  expect_equal(premult.lowrank.nmode.prod.r1(Z, lowrank.M, list(v1), 2),
               as.vector(correct.result.2))
})

Z.array <- array(1:24, dim = c(4, 3, 2))
lowrank.A2 <- list(matrix(1:8, ncol = 2),
                   matrix(1:6, ncol = 2),
                   matrix(1:4, ncol = 2))
A2 <- lowrank.A2[[1]][, 1] %o% lowrank.A2[[2]][, 1] %o% lowrank.A2[[3]][, 1] +
  lowrank.A2[[1]][, 2] %o% lowrank.A2[[2]][, 2] %o% lowrank.A2[[3]][, 2]
r1 <- list(rnorm(4), rnorm(3), rnorm(2))

test_that("r1 product with premultiplication is correct for arrays", {
  for (n in 1:3) {
    correct.result <- nmode.prod.r1(Z.array * A2, r1[-n], n)
    expect_equal(premult.lowrank.nmode.prod.r1(Z.array, lowrank.A2, r1[-n], n),
                 as.vector(correct.result))
  }
})
