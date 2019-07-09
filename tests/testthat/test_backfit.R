context("backfit parameter")

data(gtex)

# fl.none <- flashier(gtex, greedy.Kmax = 3, fit = "greedy", verbose.lvl = 0)
# fl.alt  <- flashier(gtex, greedy.Kmax = 3, backfit = "alternating", verbose.lvl = 0)
# fl.ff   <- flashier(gtex, greedy.Kmax = 3, backfit = "first.factor", verbose.lvl = 0)

# test_that("ELBOs for different backfits make sense", {
#   # unlike 'alternating', 'first.factor' doesn't include factor 2 in final backfit
#   expect_true(fl.alt$elbo > fl.ff$elbo)
#   expect_true(fl.ff$elbo > fl.none$elbo)
# })

fl.none <- flashier(gtex, greedy.Kmax = 5, fit = "greedy", verbose.lvl = 0)
fl.p <- flashier(flash.init = fl.none, fit = "backfit", backfit.order = "parallel",
                 verbose.lvl = 0)

test_that("ELBO for parallel backfit makes sense", {
  expect_true(fl.p$elbo > fl.none$elbo)
})

test_that("parallel backfit goes to convergence", {
  fl.pb <- flashier(flash.init = fl.p, fit = "backfit", verbose.lvl = 0)
  expect_equal(fl.p$elbo, fl.pb$elbo, tol = 0.001)
})
