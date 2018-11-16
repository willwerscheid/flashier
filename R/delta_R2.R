calc.delta.R2 <- function(factor, flash) {
  R <- get.R(flash)
  Y <- get.Y(flash)
  Z <- get.nonmissing(flash)
  n <- get.tau.n(flash)
  k <- get.k(factor)

  new.EF  <- as.lowrank(get.EF(factor))
  new.EF2 <- as.lowrank(get.EF2(factor))
  if (!is.null(k)) {
    old.EF  <- as.lowrank(get.EFk(flash, k))
    old.EF2 <- as.lowrank(get.EF2k(flash, k))
  }

  if (is.null(k)) {
    EF.delta.mat <- new.EF
  } else {
    EF.delta.mat <- lowrank.delta.mat(new.EF, old.EF)
  }

  if (uses.R(flash)) {
    delta.R2 <- -2 * premult.nmode.prod.r1(R, EF.delta.mat, r1.ones(flash), n)
  } else {
    delta.R2 <- -2 * premult.nmode.prod.r1(Y, EF.delta.mat, r1.ones(flash), n)
  }

  if (uses.R(flash) && is.null(k)) {
    delta.mat2 <- new.EF2
  } else if (uses.R(flash)) {
    EF2.delta.mat <- lowrank.delta.mat(new.EF2, old.EF2)
    EFs.delta.mat <- lowrank.delta.mat(lowranks.prod(old.EF, new.EF),
                                       lowrank.square(old.EF))
    delta.mat2 <- lowranks.combine(EF2.delta.mat,
                                   lowrank.sc.mult(EFs.delta.mat, -2))
  } else {
    EF.less.k <- get.EF(flash)
    if (!is.null(k))
      EF.less.k <- lowrank.drop.k(EF.less.k, k)
    EFk.mat <- lowranks.prod(new.EF, EF.less.k)
    if (is.null(k)) {
      delta.mat2 <- lowranks.combine(new.EF2,
                                     lowrank.sc.mult(EFk.mat, -2))
    } else {
      EF2.delta.mat <- lowrank.delta.mat(new.EF2, old.EF2)
      EFk.delta.mat <- lowrank.delta.mat(EFk.mat,
                                         lowranks.prod(old.EF, EF.less.k))
      delta.mat2 <- lowranks.combine(EF2.delta.mat,
                                     lowrank.sc.mult(EFk.delta.mat, -2))
    }
  }

  delta.R2 <- (delta.R2
               + premult.nmode.prod.r1(Z, delta.mat2, r1.ones(flash), n))

  if (get.tau.dim(flash) < 1)
    delta.R2 <- sum(delta.R2)

  return(delta.R2)
}
