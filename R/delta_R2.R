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
