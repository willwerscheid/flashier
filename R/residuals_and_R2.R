update.residuals <- function(flash, factor) {
  if (uses.R(flash)) {
    R <- get.R(factor)
    if (is.null(R))
      R <- calc.residuals(flash, factor)
    flash <- set.R(flash, R)
  }

  return(flash)
}

calc.residuals <- function(flash, factor = NULL) {
  old.R <- get.R(flash)

  if (!is.null(factor) && !is.new(factor)) {
    new.EF       <- as.lowrank(get.EF(factor))
    old.EF       <- as.lowrank(get.EF.k(flash, get.k(factor)))
    EF.delta.mat <- lowrank.delta.mat(new.EF, old.EF)
    new.R <- old.R - get.nonmissing(flash) * lowrank.expand(EF.delta.mat)
  } else {
    new.R <- old.R - get.nonmissing(flash) * r1.expand(get.EF(factor))
  }

  return(new.R)
}

calc.delta.R2.for.lowrank.tau <- function(factor, flash) {
  R <- get.R(flash)
  Y <- get.Y(flash)
  Z <- get.nonmissing(flash)
  n <- get.R2.n(flash)
  k <- get.k(factor)

  is.new.factor <- is.new(factor)

  new.EF  <- as.lowrank(get.EF(factor))
  new.EF2 <- as.lowrank(get.EF2(factor))
  if (is.new.factor) {
    EF.delta.mat  <- new.EF
  } else {
    old.EF        <- as.lowrank(get.EF.k(flash, k))
    EF.delta.mat  <- lowrank.delta.mat(new.EF, old.EF)
    old.EF2       <- as.lowrank(get.EF2.k(flash, k))
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
    EFprod.mat <- lowranks.prod(new.EF, get.EF(flash), broadcast = TRUE)
    ugly.mat <- lowranks.combine(new.EF2,
                                 lowrank.sc.mult(EFprod.mat, 2))
  } else { # !is.new.factor && !uses.R(flash)
    EF.less.k        <- lowrank.drop.k(get.EF(flash), k)
    newprod.mat      <- lowranks.prod(new.EF, EF.less.k, broadcast = TRUE)
    oldprod.mat      <- lowranks.prod(old.EF, EF.less.k, broadcast = TRUE)
    EFprod.delta.mat <- lowrank.delta.mat(newprod.mat, oldprod.mat)
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

calc.sum.tau.R2 <- function(flash, factor) {
  if (is.null(factor)) {
    R   <- get.R(flash)
    tau <- get.tau(flash)
  } else {
    R   <- get.R(factor)
    tau <- get.tau(factor)
    k   <- get.k(factor)
  }
  n   <- get.R2.n(flash)

  EF  <- get.EF(flash)
  EF2 <- get.EF2(flash)
  if (!is.new(factor)) {
    EF  <- lowrank.drop.k(EF, k)
    EF2 <- lowrank.drop.k(EF2, k)
  }
  EF  <- lowranks.combine(EF, as.lowrank(get.EF(factor)))
  EF2 <- lowranks.combine(EF2, as.lowrank(get.EF2(factor)))

  EFsquared <- lowrank.square(EF)

  return(sum(tau * R^2)
         + sum(premult.nmode.prod.r1(tau, EF2, r1.ones(flash), n)
               - premult.nmode.prod.r1(tau, EFsquared, r1.ones(flash), n)))
}
