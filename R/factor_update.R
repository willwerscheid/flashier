update.factor <- function(factor, flash) {
  if (is.zero(factor))
    return(factor)

  for (n in 1:get.dim(flash))
    if (!is.zero(factor))
      factor <- update.factor.one.n(factor, n, flash)

  delta.R2 <- calc.delta.R2(factor, flash)

  factor <- set.delta.R2(factor, delta.R2)
  factor <- set.est.tau(factor, calc.est.tau(flash, delta.R2))
  factor <- set.obj(factor, calc.obj(flash, factor))

  return(factor)
}

update.factor.one.n <- function(factor, n, flash) {
  ebnm.res <- solve.ebnm(factor, n, flash)

  factor <- set.EF(factor, ebnm.res$postmean, n)
  factor <- set.EF2(factor, ebnm.res$postmean2, n)
  factor <- set.KL(factor, ebnm.res$KL, n)
  factor <- set.g(factor, ebnm.res$fitted_g, n)

  if (all(ebnm.res$postmean == 0))
    factor <- set.to.zero(factor)

  return(factor)
}

solve.ebnm <- function(factor, n, flash) {
  ebnm.args   <- calc.ebnm.args(factor, n, flash)
  ebnm.res    <- do.call(get.ebnm.fn(flash),
                         c(ebnm.args, list(get.ebnm.param(flash))))
  ebnm.res$KL <- (ebnm.res$penloglik
                  - normal.means.loglik(ebnm.args$x, ebnm.args$s,
                                        ebnm.res$postmean, ebnm.res$postmean2))
  return(ebnm.res)
}

calc.ebnm.args <- function(factor, n, flash) {
  s2 <- calc.s2(factor, n, flash)
  x  <- calc.x(factor, n, flash, s2)

  return(list(x = x, s = sqrt(s2)))
}

calc.s2 <- function(factor, n, flash) {
  Z          <- get.nonmissing(flash)
  factor.EF2 <- get.EF2(factor)
  tau        <- get.tau.lowrank(flash, est.tau = get.est.tau(factor))

  s2 <- 1 / premult.nmode.prod.r1(Z, tau, factor.EF2[-n], n)

  return(pmax(s2, 0))
}

calc.x <- function(factor, n, flash, s2) {
  R         <- get.R(flash)
  Y         <- get.Y(flash)
  Z         <- get.nonmissing(flash)
  factor.EF <- get.EF(factor)
  flash.EF  <- get.EF(flash)
  tau       <- get.tau.lowrank(flash, est.tau = get.est.tau(factor))
  k         <- get.k(factor) # set during backfitting

  if (uses.R(flash)) {
    x <- premult.nmode.prod.r1(R, tau, factor.EF[-n], n)
    if (!is.null(k))
      x <- x + nmode.prod.r1(elemwise.prod.lowrank.r1(tau, get.EFk(flash, k)),
                             factor.EF[-n], n)
    x <- s2 * x
  } else {
    flash.EF.tau <- lowranks.prod(tau, flash.EF, broadcast = TRUE)
    if (!is.null(k))
      flash.EF.tau <- lowrank.drop.k(flash.EF.tau, k)
    x <- s2 * (premult.nmode.prod.r1(Y, tau, factor.EF[-n], n)
               - premult.nmode.prod.r1(Z, flash.EF.tau, factor.EF[-n], n))
  }

  return(x)
}
