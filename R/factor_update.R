update.factor <- function(factor, flash) {
  if (is.zero(factor))
    return(factor)

  for (n in 1:get.dim(flash))
    if (!is.zero(factor) && !all.fixed(factor, n))
      factor <- update.factor.one.n(factor, n, flash)

  factor <- update.R2.tau.and.obj(factor, flash)

  return(factor)
}

update.factor.one.n <- function(factor, n, flash) {
  ebnm.res <- solve.ebnm(factor, n, flash)

  if (only.update.subset(factor, n, flash)) {
    idx.subset          <- get.idx.subset(factor)
    new.EF              <- get.EF(factor, n)
    new.EF[idx.subset]  <- ebnm.res$postmean
    new.EF2             <- get.EF2(factor, n)
    new.EF2[idx.subset] <- ebnm.res$postmean2
  } else {
    new.EF  <- ebnm.res$postmean
    new.EF2 <- ebnm.res$postmean2
  }

  factor <- set.EF(factor, new.EF, n)
  factor <- set.EF2(factor, new.EF2, n)
  factor <- set.KL(factor, ebnm.res$KL, n)
  factor <- set.g(factor, ebnm.res$fitted_g, n)

  if (all(new.EF == 0))
    factor <- set.to.zero(factor)

  return(factor)
}

update.R2.tau.and.obj <- function(factor, flash) {
  factor <- update.tau(factor, flash)
  factor <- set.obj(factor, calc.obj(flash, factor))
  factor <- set.to.valid(factor)
  return(factor)
}

solve.ebnm <- function(factor, n, flash, return.sampler = FALSE) {
  fix.dim <- get.fix.dim(factor)
  if (use.subsetted.flash.data(factor, n))
    factor <- add.subset.data(factor, flash, fix.dim, get.idx.subset(factor))

  ebnm.fn     <- get.ebnm.fn(flash, factor, n)
  ebnm.args   <- calc.ebnm.args(factor, n, flash, return.sampler)
  ebnm.param  <- get.ebnm.param(flash, factor, n)

  g <- get.g(factor, n)
  if (return.sampler && !is.null(g)) {
    ebnm.param <- modifyList(ebnm.param, list(g = g,
                                              fixg = TRUE,
                                              output = "post_sampler"))
  } else if (!is.new(factor) && warmstart.backfits(flash) && !is.null(g)) {
    ebnm.param <- c(ebnm.param, list(g = g))
  }
  ebnm.res <- do.call(ebnm.fn, c(ebnm.args, list(ebnm.param)))

  if (!return.sampler) {
    ebnm.res$KL <- (ebnm.res$penloglik
                    - normal.means.loglik(ebnm.args$x,
                                          ebnm.args$s,
                                          ebnm.res$postmean,
                                          ebnm.res$postmean2))
  }

  return(ebnm.res)
}

calc.ebnm.args <- function(factor, n, flash, return.sampler) {
  tau <- get.tau.for.ebnm.calc(flash, tau = get.tau(factor))
  if (n %in% get.fix.dim(factor))
    tau <- full.or.lowrank.subset(tau, n, get.idx.subset(factor))

  s2 <- calc.s2(factor, n, flash, tau)
  x  <- calc.x(factor, n, flash, s2, tau)

  if (add.fixed.to.ebnm.args(factor, n, flash, return.sampler)) {
    idx.subset         <- get.idx.subset(factor)
    all.x              <- get.EF(factor, n)
    all.x[idx.subset]  <- x
    all.s2             <- rep(0, length(all.x))
    all.s2[idx.subset] <- s2
  } else {
    all.x  <- x
    all.s2 <- s2
  }

  return(list(x = all.x, s = sqrt(all.s2)))
}

calc.s2 <- function(factor, n, flash, tau) {
  if (use.subsetted.flash.data(factor, n)) {
    Z <- get.Z.subset(factor)
  } else {
    Z <- get.nonmissing(flash)
  }

  factor.EF2 <- get.EF2(factor)
  if (n %in% get.fix.dim(factor))
    factor.EF2 <- r1.subset(factor.EF2, n,  get.idx.subset(factor))

  if (is.tau.lowrank(flash)) {
    s2 <- 1 / premult.nmode.prod.r1(Z, tau, factor.EF2[-n], n)
  } else {
    # If tau is full-rank, then it has already been multiplied by Z:
    s2 <- 1 / nmode.prod.r1(tau, factor.EF2[-n], n)
  }

  return(pmax(s2, 0))
}

calc.x <- function(factor, n, flash, s2, tau) {
  if (use.subsetted.flash.data(factor, n)) {
    R        <- get.R.subset(factor)
    Y        <- get.Y.subset(factor)
    Z        <- get.Z.subset(factor)
    flash.EF <- get.EF.subset(factor)
  } else {
    R        <- get.R(flash)
    Y        <- get.Y(flash)
    Z        <- get.nonmissing(flash)
    flash.EF <- get.EF(flash)
  }

  k <- get.k(factor)
  if (uses.R(flash) && !is.new(factor)) {
    flash.EFk <- get.EF.k(flash, k)
    if (use.subsetted.flash.data(factor, n))
      flash.EFk <- r1.subset(flash.EFk, n, get.idx.subset(factor))
  }

  factor.EF <- get.EF(factor)
  if (n %in% get.fix.dim(factor))
    factor.EF <- r1.subset(factor.EF, n, get.idx.subset(factor))

  if (uses.R(flash)) {
    x <- premult.nmode.prod.r1(R, tau, factor.EF[-n], n)
    if (!is.new(factor)) {
      if (is.tau.lowrank(flash)) {
        EFk.tau <- elemwise.prod.lowrank.r1(tau, flash.EFk)
        x <- x + premult.nmode.prod.r1(Z, EFk.tau, factor.EF[-n], n)
      } else {
        EFk.tau <- elemwise.prod.fullrank.r1(tau, flash.EFk)
        x <- x + nmode.prod.r1(EFk.tau, factor.EF[-n], n)
      }
    }
  } else {
    x <- premult.nmode.prod.r1(Y, tau, factor.EF[-n], n)
    # Note that whenever Y is used, tau is low-rank.
    flash.EF.tau <- lowranks.prod(tau, flash.EF, broadcast = TRUE)
    if (!is.new(factor))
      flash.EF.tau <- lowrank.drop.k(flash.EF.tau, k)
    x <- x - premult.nmode.prod.r1(Z, flash.EF.tau, factor.EF[-n], n)
  }
  x <- s2 * x

  return(x)
}

only.update.subset <- function(factor, n, flash) {
  return((n %in% get.fix.dim(factor)) && !use.fixed.to.est.g(flash))
}

use.subsetted.flash.data <- function(factor, n) {
  return((n %in% get.fix.dim(factor)) && !all.fixed(factor, n))
}

add.fixed.to.ebnm.args <- function(factor, n, flash, return.sampler) {
  return((n %in% get.fix.dim(factor))
          && (use.fixed.to.est.g(flash) || return.sampler))
}
