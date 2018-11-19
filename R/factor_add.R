# TODO remove default tol and maxiter after tests are removed
add.next.factor <- function(flash,
                            greedy.tol = 1e-2,
                            greedy.maxiter = 200,
                            init.tol = 1e-2,
                            init.maxiter = 100) {
  factor <- init.factor(flash, init.tol, init.maxiter)
  if (is.fixed(factor)) {
    flash <- add.new.factor.to.flash(factor, flash)
  } else {
    flash <- add.greedy(factor, flash, greedy.tol, greedy.maxiter)
  }
  return(flash)
}

add.greedy <- function(factor, flash, tol, maxiter) {
  factor <- optimize.it(factor,
                        update.fn = update.factor,
                        update.args = list(flash = flash),
                        obj.fn = calc.obj.diff,
                        tol = tol,
                        maxiter = maxiter)
  if (ok.to.add.factor(factor, flash)) {
    flash <- add.new.factor.to.flash(factor, flash)
  } else {
    flash <- set.greedy.fail.flag(flash)
  }
  # TODO: what if objective decreases?

  return(flash)
}

ok.to.add.factor <- function(factor, flash) {
  return(get.obj(factor) > get.obj(flash) || !is.obj.valid(flash, factor))
}

add.new.factor.to.flash <- function(factor, flash) {
  flash <- add.factor.to.EF(flash, get.EF(factor))
  flash <- add.factor.to.EF2(flash, get.EF2(factor))
  flash <- add.factor.to.KL(flash, get.KL(factor))
  flash <- add.factor.to.g(flash, get.g(factor))
  if (!is.null(get.delta.R2(factor)))
    flash <- set.R2(flash, get.R2(flash) + get.delta.R2(factor))
  flash <- set.tau(flash, get.tau(factor))
  flash <- set.obj(flash, get.obj(factor))
  flash <- add.is.zero(flash, FALSE)
  flash <- add.is.valid(flash, is.valid(factor))
  flash <- update.residuals(flash, factor)

  return(flash)
}
