add.next.factor <- function(flash, tol = 1e-2) {
  factor <- init.factor(flash)
  if (is.fixed(factor)) {
    flash <- add.new.factor.to.flash(factor, flash)
  } else {
    flash <- add.greedy(factor, flash, tol)
  }
  return(flash)
}

add.greedy <- function(factor, flash, tol = 1e-2) {
  factor <- optimize.it(factor,
                        update.fn = update.factor,
                        update.args = list(flash = flash),
                        obj.fn = calc.obj.diff,
                        tol = tol)
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
  flash <- set.R2(flash, get.R2(flash) + get.delta.R2(factor))
  flash <- set.est.tau(flash, get.est.tau(factor))
  flash <- set.obj(flash, get.obj(factor))
  flash <- add.is.zero(flash, FALSE)
  flash <- add.is.valid(flash, is.valid(factor))

  if (uses.R(flash)) {
    R <- get.R(flash) - get.nonmissing(flash) * r1.expand(get.EF(factor))
    flash <- set.R(flash, R)
  }

  return(flash)
}
