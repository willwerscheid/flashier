context("backfit parameter")

data(gtex)

fl.none <- flashier(gtex, greedy.Kmax = 3, backfit = "none", verbose.lvl = 0)
fl.alt  <- flashier(gtex, greedy.Kmax = 3, backfit = "alternating", verbose.lvl = 0)
fl.ff   <- flashier(gtex, greedy.Kmax = 3, backfit = "first.factor", verbose.lvl = 0)

test_that("ELBOs for different backfits make sense", {
  # unlike 'alternating', 'first.factor' doesn't include factor 2 in final backfit
  expect_true(fl.alt$elbo > fl.ff$elbo)
  expect_true(fl.ff$elbo > fl.none$elbo)
})
