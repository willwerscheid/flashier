add.new.factor.to.flash <- function(factor, flash) {
  flash <- add.factor.to.EF(flash, get.EF(factor))
  flash <- add.factor.to.EF2(flash, get.EF2(factor))
  flash <- add.factor.to.KL(flash, get.KL(factor))
  flash <- add.factor.to.g(flash, get.g(factor))
  flash <- add.factor.to.ebnm.fn(flash, get.ebnm.fn(flash, factor))
  flash <- set.tau(flash, get.tau(factor))
  flash <- set.obj(flash, get.obj(factor))
  flash <- add.is.zero(flash, FALSE)
  flash <- add.is.valid(flash, is.valid(factor))
  flash <- add.exclusions(flash, get.exclusions(factor))

  if (uses.R(flash))
    flash <- update.R(flash, factor)

  if (is.tau.simple(flash)) {
    flash <- set.R2(flash, get.R2(flash) + get.delta.R2(factor))
    flash <- set.est.tau(flash, get.est.tau(factor))
  }

  return(flash)
}
