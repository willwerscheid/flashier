#' @exportS3Method NULL
update.R <- function(flash, factor) {
  R <- get.R(factor)
  if (is.null(R))
    R <- calc.residuals(flash, factor)
  flash <- set.R(flash, R)

  return(flash)
}

calc.residuals <- function(flash, factor = NULL) {
  R <- get.R(flash)

  if (!is.null(factor) && !is.new(factor)) {
    new.EF <- as.lowrank(get.EF(factor))
    old.EF <- as.lowrank(get.EF.k(flash, get.k(factor)))
    EF.delta.mat <- lowrank.delta.mat(new.EF, old.EF)
    R <- R - get.nonmissing(flash) * lowrank.expand(EF.delta.mat)
  } else if (!is.null(factor)) {
    R <- R - get.nonmissing(flash) * r1.expand(get.EF(factor))
  }

  return(R)
}

calc.Y2 <- function(flash) {
  Y <- get.Y(flash)
  n  <- get.R2.n(flash)

  Y2 <- sq.nmode.prod.r1(Y, r1.ones(flash), n)
  if (store.R2.as.scalar(flash))
    Y2 <- sum(Y2)

  return(Y2)
}

# Mainly used to initialize tau, but also used for parallel backfits.
calc.R2 <- function(flash) {
  R   <- get.R(flash)
  Y   <- get.Y(flash)
  Y2  <- get.Y2(flash)
  Z   <- get.nonmissing(flash)
  EF  <- get.EF(flash)
  EF2 <- get.EF2(flash)
  n   <- get.R2.n(flash)

  if (uses.R(flash)) {
    R2 <- nmode.prod.r1(R^2, r1.ones(flash), n)

    if (store.R2.as.scalar(flash)) {
      R2 <- sum(R2)
    }
  } else {
    if (is.null(EF)) {
      Y.EF <- 0
      EFsq <- 0
    } else {
      Y.EF <- premult.nmode.prod.r1(Y, EF, r1.ones(flash), n)
      if (!any_missing(flash) && store.R2.as.scalar(flash)) {
        EFsq <- sum(Reduce(`*`, lapply(EF, crossprod)))
      } else if (get.dim(flash) == 2 && identical(Z, 1)) {
        EFsq <- apply(EF[[n]], 1, tcrossprod) * as.vector(crossprod(EF[[-n]]))
        if (is.matrix(EFsq)) {
          EFsq <- colSums(EFsq)
        }
      } else {
        # TODO: fix for tensors
        EFsq <- premult.nmode.prod.r1(Z, lowrank.expand(EF)^2, r1.ones(flash), n)
      }
    }

    if (store.R2.as.scalar(flash)) {
      R2 <- Y2 - 2 * sum(Y.EF) + sum(EFsq)
    } else {
      R2 <- Y2 - 2 * Y.EF + EFsq
    }
  }

  tmp <- premult.nmode.prod.r1(Z, EF2, r1.ones(flash), n)
  tmp <- tmp - premult.nmode.prod.r1(Z, lowrank.square(EF), r1.ones(flash), n)
  if (store.R2.as.scalar(flash)) {
    R2 <- R2 + sum(tmp)
  } else {
    R2 <- R2 + tmp
  }

  return(R2)
}

# Used to update tau when tau is simple.
calc.delta.R2 <- function(factor, flash) {
  R <- get.R(flash)
  Y <- get.Y(flash)
  Z <- get.nonmissing(flash)
  n <- get.R2.n(flash)
  k <- get.k(factor)

  is.new.factor <- is.new(factor)

  new.EF  <- as.lowrank(get.EF(factor))
  new.EF2 <- as.lowrank(get.EF2(factor))
  if (is.new.factor) {
    EF.delta.mat <- new.EF
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
    if (is.null(EF.less.k)) {
      ugly.mat <- NULL
    } else {
      newprod.mat      <- lowranks.prod(new.EF, EF.less.k, broadcast = TRUE)
      oldprod.mat      <- lowranks.prod(old.EF, EF.less.k, broadcast = TRUE)
      EFprod.delta.mat <- lowrank.delta.mat(newprod.mat, oldprod.mat)
      ugly.mat <- lowranks.combine(EF2.delta.mat,
                                   lowrank.sc.mult(EFprod.delta.mat, 2))
    }
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

# Used to calculate objective when variance type is NULL.
calc.tau.R2 <- function(flash, factor, n) {
  if (is.null(factor)) {
    R   <- get.R(flash)
    tau <- get.tau(flash)
  } else {
    R   <- get.R(factor)
    tau <- get.tau(factor)
    k   <- get.k(factor)
  }

  EF  <- get.EF(flash)
  EF2 <- get.EF2(flash)
  if (!is.new(factor)) {
    EF  <- lowrank.drop.k(EF, k)
    EF2 <- lowrank.drop.k(EF2, k)
  }
  EF  <- lowranks.combine(EF, as.lowrank(get.EF(factor)))
  EF2 <- lowranks.combine(EF2, as.lowrank(get.EF2(factor)))

  EFsquared <- lowrank.square(EF)

  if (is.tau.lowrank(flash)) {
    tau <- as.r1.tau(tau, flash, n)
    tau.R2 <- (nmode.prod.r1(R^2, tau, n)
               + nmode.prod.r1(EF2, tau, n)
               - nmode.prod.r1(EFsquared, tau, n))
  } else {
    tau.R2 <- (nmode.prod.r1(tau * R^2, r1.ones(flash), n)
               + premult.nmode.prod.r1(tau, EF2, r1.ones(flash), n)
               - premult.nmode.prod.r1(tau, EFsquared, r1.ones(flash), n))
  }
  return(tau.R2)
}
