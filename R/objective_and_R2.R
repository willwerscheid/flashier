update.R2.tau.and.obj <- function(factor, flash) {
  delta.R2 <- calc.delta.R2(factor, flash)
  factor   <- set.delta.R2(factor, delta.R2)
  factor   <- set.est.tau(factor, calc.est.tau(flash, delta.R2))
  factor   <- set.obj(factor, calc.obj(flash, factor))
  factor   <- set.to.valid(factor)

  return(factor)
}

calc.delta.R2 <- function(factor, flash) {
  R <- get.R(flash)
  Y <- get.Y(flash)
  Z <- get.nonmissing(flash)
  n <- get.tau.n(flash)
  k <- get.k(factor)

  is.new.factor <- is.new(factor)

  new.EF  <- as.lowrank(get.EF(factor))
  new.EF2 <- as.lowrank(get.EF2(factor))
  if (is.new.factor) {
    EF.delta.mat  <- new.EF
  } else {
    old.EF        <- as.lowrank(get.EFk(flash, k))
    EF.delta.mat  <- lowrank.delta.mat(new.EF, old.EF)
    old.EF2       <- as.lowrank(get.EF2k(flash, k))
    EF2.delta.mat <- lowrank.delta.mat(new.EF2, old.EF2)
  }

  if (uses.R(flash) && is.new.factor) {
    ugly.mat <- new.EF2
  } else if (uses.R(flash)) { # && !is.new.factor
    EFprod.delta.mat <- lowrank.delta.mat(lowranks.prod(old.EF, new.EF),
                                          lowrank.square(old.EF))
    ugly.mat <- lowranks.combine(EF2.delta.mat,
                                 lowrank.sc.mult(EFprod.delta.mat, -2))
  } else if (is.new.factor) { # && !uses.R(flash)
    EFprod.mat <- lowranks.prod(new.EF, get.EF(flash))
    ugly.mat <- lowranks.combine(new.EF2,
                                 lowrank.sc.mult(EFprod.mat, 2))
  } else { # !is.new.factor && !uses.R(flash)
    EF.less.k        <- lowrank.drop.k(get.EF(flash), k)
    EFprod.delta.mat <- lowrank.delta.mat(lowranks.prod(new.EF, EF.less.k),
                                          lowranks.prod(old.EF, EF.less.k))
    ugly.mat <- lowranks.combine(EF2.delta.mat,
                                 lowrank.sc.mult(EFprod.delta.mat, 2))
  }

  if (uses.R(flash)) {
    delta.R2 <- -2 * premult.nmode.prod.r1(R, EF.delta.mat, r1.ones(flash), n)
  } else {
    delta.R2 <- -2 * premult.nmode.prod.r1(Y, EF.delta.mat, r1.ones(flash), n)
  }
  delta.R2 <- delta.R2 + premult.nmode.prod.r1(Z, ugly.mat, r1.ones(flash), n)

  if (store.R2.as.scalar(flash))
    delta.R2 <- sum(delta.R2)

  return(delta.R2)
}

calc.est.tau <- function(flash, delta.R2 = 0) {
  return(get.n.nonmissing(flash) / (get.R2(flash) + delta.R2))
}

calc.obj <- function(flash, factor = NULL) {
  n.nonmissing <- get.n.nonmissing(flash)
  k            <- get.k(factor)

  KL <- sum(unlist(get.KL(flash)))
  if (!is.null(factor)) {
    KL <- KL + sum(unlist(get.KL(factor)))
    if (!is.new(factor))
      KL <- KL - sum(get.KLk(flash, k))
    est.tau <- get.est.tau(factor)
  } else {
    est.tau <- get.est.tau(flash)
  }

  # TODO: fixed S
  return(KL - 0.5 * sum(n.nonmissing * (log(2 * pi / est.tau) + 1)))
}

normal.means.loglik <- function(x, s, Et, Et2) {
  idx <- is.finite(s) && s > 0
  x   <- x[idx]
  s   <- s[idx]
  Et  <- Et[idx]
  Et2 <- Et2[idx]
  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}
