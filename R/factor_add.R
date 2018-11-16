add.greedy <- function(flash, tol = 1e-2) {
  factor <- init.factor(flash)
  factor <- optimize.it(factor,
                        update.fn = update.factor,
                        update.args = list(flash = flash),
                        obj.fn = calc.obj.diff,
                        tol = tol)
  if (get.obj(factor) > get.obj(flash)) {
    flash <- add.new.factor.to.flash(factor, flash)
  } else {
    flash <- set.greedy.fail.flag(flash)
  }
  # TODO: what if objective decreases?

  return(flash)
}

add.new.factor.to.flash <- function(factor, flash) {
  flash <- add.factor.to.EF(flash, get.EF(factor))
  flash <- add.factor.to.EF2(flash, get.EF2(factor))
  flash <- add.factor.to.KL(flash, get.KL(factor))
  flash <- add.factor.to.g(flash, get.g(factor))
  flash <- set.R2(flash, get.R2(flash) + get.delta.R2(factor))
  flash <- set.est.tau(flash, get.est.tau(factor))
  flash <- set.obj(flash, get.obj(factor))
  flash <- set.is.zero(flash, c(is.zero(flash), FALSE))

  if (uses.R(flash))
    flash <- set.R(flash, get.R(flash) - r1.expand(get.EF(factor)))

  return(flash)
}
