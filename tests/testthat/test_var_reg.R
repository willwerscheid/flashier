context("variance types")

ebpm_exponential <- function(x, s = 1) {
  neg_llik_fn <- function(log.lambda) {
    return(-sum(dexp(x, rate = exp(log.lambda) / s, log = TRUE)))
  }
  opt_res <- optimize(
    neg_llik_fn,
    lower = -100,
    upper = 100,
    tol = 1e-6
  )

  est.lambda <- exp(opt_res$minimum)
  llik <- -opt_res$objective

  posterior <- data.frame(
    mean = (x + 1) / (s + est.lambda),
    mean_log = digamma(x + 1) - log(s + est.lambda)
  )

  return(list(posterior = posterior, log_likelihood = llik))
}

set.seed(666)

n <- 40
p <- 60

LF <- outer(rep(1, n), rep(1, p))
M <- LF + 0.1 * rnorm(n * p)

test_that("variance regularization only works with simple variance types", {
  expect_error(flash.init(M, S = 0.2, var.type = 0, var.reg.fn = ebpm_exponential))
  expect_error(flash.init(M, var.type = c(1, 2), var.reg.fn = ebpm_exponential))
})

f.noreg <- flash.init(M, var.type = 1) %>%
  flash.set.verbose(3) %>%
  flash.add.greedy(Kmax = 1)
f.reg <- flash.init(M, var.type = 1, var.reg.fn = ebpm_exponential) %>%
  flash.set.verbose(3) %>%
  flash.add.greedy(Kmax = 1)

test_that("variance regularization changes estimates", {
  expect_true(all(f.noreg$residuals.sd != f.reg$residuals.sd))
})

test_that("variance regularization estimates make sense", {
  expect_true(max(abs(f.noreg$residuals.sd - f.reg$residuals.sd)) < 0.01)
  expect_true(sum(f.reg$flash.fit$R2) - sum(f.noreg$flash.fit$R2) < 0.01)
  expect_true(f.reg$elbo - f.reg$flash.fit$KL.tau < f.noreg$elbo)
})
