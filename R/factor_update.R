update.factor <- function(factor, flash) {
  if (is.zero(factor))
    return(factor)

  for (n in 1:get.dim(flash))
    if (!is.zero(factor) && !are.all.fixed(factor, n))
      factor <- update.factor.one.n(factor, n, flash)

  delta.R2 <- calc.delta.R2(factor, flash)

  factor <- set.delta.R2(factor, delta.R2)
  factor <- set.est.tau(factor, calc.est.tau(flash, delta.R2))
  factor <- set.obj(factor, calc.obj(flash, factor))
  factor <- set.to.valid(factor)

  return(factor)
}

update.factor.one.n <- function(factor, n, flash) {
  ebnm.res <- solve.ebnm(factor, n, flash)

  if (identical(get.fix.dim(factor), n) && !incl.fixed.in.prior.est(flash)) {
    new.EF  <- get.EF(factor, n)
    new.EF2 <- get.EF2(factor, n)
    new.EF[get.idx.subset(factor)]  <- ebnm.res$postmean
    new.EF2[get.idx.subset(factor)] <- ebnm.res$postmean2
  } else {
    new.EF  <- ebnm.res$postmean
    new.EF2 <- ebnm.res$postmean2
  }

  factor <- set.EF(factor, new.EF, n)
  factor <- set.EF2(factor, new.EF2, n)
  factor <- set.KL(factor, ebnm.res$KL, n)
  factor <- set.g(factor, ebnm.res$fitted_g, n)

  if (all(ebnm.res$postmean == 0))
    factor <- set.to.zero(factor)

  return(factor)
}

solve.ebnm <- function(factor, n, flash) {
  fix.dim <- get.fix.dim(factor)
  if (identical(fix.dim, n) && !are.all.fixed(factor, n))
    factor <- add.subset.data(factor, flash, fix.dim, get.idx.subset(factor))

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

  if (identical(get.fix.dim(factor), n) && incl.fixed.in.prior.est(flash)) {
    all.x <- get.EF(factor)[[n]]
    all.x[get.idx.subset(factor)] <- x
    x <- all.x
    all.s2 <- rep(0, length(all.x))
    all.s2[get.idx.subset(factor)] <- s2
    s <- all.s2
  }

  return(list(x = x, s = sqrt(s2)))
}

calc.s2 <- function(factor, n, flash) {
  if (identical(get.fix.dim(factor), n) && !are.all.fixed(factor, n)) {
    Z <- get.Z.subset(factor)
  } else {
    Z <- get.nonmissing(flash)
  }

  factor.EF2 <- get.EF2(factor)
  tau        <- get.tau.lowrank(flash, est.tau = get.est.tau(factor))
  if (identical(get.fix.dim(factor), n)) {
    factor.EF2 <- r1.subset(factor.EF2, n, get.idx.subset(factor))
    tau        <- lowrank.subset(tau, n, get.idx.subset(factor))
  }

  s2 <- 1 / premult.nmode.prod.r1(Z, tau, factor.EF2[-n], n)

  return(pmax(s2, 0))
}

calc.x <- function(factor, n, flash, s2) {
  k <- get.k(factor) # set during backfitting

  if (identical(get.fix.dim(factor), n) && !are.all.fixed(factor, n)) {
    R        <- get.R.subset(factor)
    Y        <- get.Y.subset(factor)
    Z        <- get.Z.subset(factor)
    flash.EF <- get.EF.subset(factor)
    if (uses.R(flash) && !is.null(k)) {
      flash.EFk <- get.EFk(flash, k)
      flash.EFk <- r1.subset(flash.EFk, n, get.idx.subset(factor))
    }
  } else {
    R        <- get.R(flash)
    Y        <- get.Y(flash)
    Z        <- get.nonmissing(flash)
    flash.EF <- get.EF(flash)
    if (uses.R(flash) && !is.null(k))
      flash.EFk <- get.EFk(flash, k)
  }

  factor.EF <- get.EF(factor)
  tau       <- get.tau.lowrank(flash, est.tau = get.est.tau(factor))
  if (identical(get.fix.dim(factor), n)) {
    factor.EF <- r1.subset(factor.EF, n, get.idx.subset(factor))
    tau       <- lowrank.subset(tau, n, get.idx.subset(factor))
  }

  if (uses.R(flash)) {
    x <- premult.nmode.prod.r1(R, tau, factor.EF[-n], n)
    if (!is.null(k)) {
      EFk.tau <- elemwise.prod.lowrank.r1(tau, flash.EFk)
      x <- x + premult.nmode.prod.r1(Z, EFk.tau, factor.EF[-n], n)
    }
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
