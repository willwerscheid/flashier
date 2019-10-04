context("backfit parameter")

data(gtex)

fl.none <- flash.init(gtex) %>%
  flash.set.verbose(0) %>%
  flash.add.greedy(K = 5)

fl.p <- flash.backfit(fl.none, method = "parallel")

test_that("ELBO for parallel backfit makes sense", {
  expect_true(fl.p$elbo > fl.none$elbo)
})

test_that("parallel backfit goes to convergence", {
  fl.pb <- flash.backfit(fl.p)
  expect_equal(fl.p$elbo, fl.pb$elbo, tol = 0.001)
})
