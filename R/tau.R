# Functions for estimating tau ------------------------------------------------

init.tau <- function(flash) {
  if (is.tau.simple(flash)) {
    flash <- init.simple.tau(flash)
  } else if (is.var.type.zero(flash)) {
    flash <- init.zero.tau(flash)
  } else if (is.var.type.kronecker(flash)) {
    flash <- init.kronecker.tau(flash)
  } else if (is.var.type.noisy(flash)) {
    flash <- init.noisy.tau(flash)
  } else if (is.var.type.noisy.kron(flash)) {
    flash <- init.noisy.kron.tau(flash)
  } else {
    # This error should never occur:
    stop("The requested variance structure has not yet been implemented.")
  }

  return(flash)
}

update.tau <- function(factor, flash) {
  if (is.tau.simple(flash)) {
    factor <- update.simple.tau(factor, flash)
  } else if (is.var.type.zero(flash)) {
    factor <- update.zero.tau(factor, flash)
  } else if (is.var.type.kronecker(flash)) {
    factor <- update.kronecker.tau(factor, flash)
  } else if (is.var.type.noisy(flash)) {
    factor <- update.noisy.tau(factor, flash)
  } else if (is.var.type.noisy.kron(flash)) {
    factor <- update.noisy.kron.tau(factor, flash)
  } else {
    stop("The requested variance structure has not yet been implemented.")
  }

  return(factor)
}

# Mainly used to initialize tau, but also used for parallel backfits.
init.simple.tau <- function(flash) {
  if (is.null(get.Y2(flash)) && !uses.R(flash)) {
    flash <- set.Y2(flash, calc.Y2(flash))
  }
  flash   <- set.R2(flash, calc.R2(flash))
  est.tau <- estimate.simple.tau(flash)
  flash   <- set.est.tau(flash, est.tau)
  flash   <- set.tau(flash, tau.from.given.and.est(flash, est.tau))

  return(flash)
}

update.simple.tau <- function(factor, flash) {
  delta.R2 <- calc.delta.R2(factor, flash)
  est.tau  <- estimate.simple.tau(flash, delta.R2)
  tau      <- tau.from.given.and.est(flash, est.tau)
  factor   <- set.delta.R2(factor, delta.R2)
  factor   <- set.est.tau(factor, est.tau)
  factor   <- set.tau(factor, tau)

  return(factor)
}

init.zero.tau <- function(flash) {
  if (is.null(get.tau(flash)))
    flash <- set.tau(flash, get.given.tau(flash))

  return(flash)
}

update.zero.tau <- function(factor, flash) {
  if (uses.R(flash))
    factor <- set.R(factor, calc.residuals(flash, factor))

  factor <- set.tau(factor, get.given.tau(flash))

  return(factor)
}

init.tau.at.one <- function(flash) {
  return(lapply(get.dims(flash), function(dim) rep(1, dim)))
}

init.kronecker.tau <- function(flash) {
  if (is.null(get.tau(flash)))
    flash <- set.tau(flash, init.tau.at.one(flash))

  flash <- set.tau(flash, estimate.kronecker.tau(flash))

  return(flash)
}

update.kronecker.tau <- function(factor, flash) {
  factor <- set.tau(factor, estimate.kronecker.tau(flash, factor))

  return(factor)
}

init.noisy.tau <- function(flash) {
  noisy.tau <- estimate.noisy.tau(flash)
  flash     <- set.sum.tau.R2(flash, noisy.tau$sum.tau.R2)
  flash     <- set.tau(flash, noisy.tau$tau)

  return(flash)
}

update.noisy.tau <- function(factor, flash) {
  if (uses.R(flash))
    factor <- set.R(factor, calc.residuals(flash, factor))

  noisy.tau <- estimate.noisy.tau(flash, factor)
  factor    <- set.tau(factor, noisy.tau$tau)
  factor    <- set.sum.tau.R2(factor, noisy.tau$sum.tau.R2)

  return(factor)
}

init.noisy.kron.tau <- function(flash) {
  if (is.null(get.est.S2(flash)))
    flash   <- set.est.S2(flash, init.tau.at.one(flash))

  noisy.tau <- estimate.noisy.kron.tau(flash)
  flash     <- set.sum.tau.R2(flash, noisy.tau$sum.tau.R2)
  flash     <- set.est.S2(flash, noisy.tau$est.S2)
  flash     <- set.tau(flash, noisy.tau$tau)

  return(flash)
}

update.noisy.kron.tau <- function(factor, flash) {
  if (uses.R(flash))
    factor <- set.R(factor, calc.residuals(flash, factor))

  noisy.tau <- estimate.noisy.kron.tau(flash, factor)
  factor <- set.tau(factor, noisy.tau$tau)
  factor <- set.est.S2(factor, noisy.tau$est.S2)
  factor <- set.sum.tau.R2(factor, noisy.tau$sum.tau.R2)

  return(factor)
}

estimate.simple.tau <- function(flash, delta.R2 = 0) {
  R2 <- pmax(get.R2(flash), .Machine$double.eps)
  return(get.n.nonmissing(flash) / (R2 + delta.R2))
}

tau.from.given.and.est <- function(flash, est.tau) {
  given.tau <- get.given.tau(flash)

  if (is.var.type.zero(flash)) {
    # All variance is pre-specified ("zero" variance type):
    tau <- given.tau
  } else if (is.null(given.tau)) {
    # All variance is estimated:
    tau <- est.tau
  } else {
    # Otherwise both types of variance are used. In the simple case here, the
    #   total variance is estimated and the given variance serves as a floor:
    tau <- pmin(given.tau, est.tau)
  }

  return(tau)
}

estimate.kronecker.tau <- function(flash, factor = NULL) {
  # tol and maxiter are hardcoded.
  tol <- sqrt(.Machine$double.eps)
  maxiter <- 100

  tau <- get.tau(factor)
  if (is.null(tau))
    tau <- get.tau(flash)
  tau <- as.r1.tau(tau, flash)
  tau.dim <- get.est.tau.dim(flash)
  kron.nonmissing <- get.kron.nonmissing(flash)

  Z   <- get.nonmissing(flash)
  EF  <- get.new.EF(flash, factor)
  EF2 <- get.new.EF2(flash, factor)

  R2 <- get.new.Rsquared(flash, factor, EF)
  R2 <- pmax(R2, sqrt(.Machine$double.eps))

  max.chg <- Inf
  iter <- 0
  while (max.chg > tol && iter < maxiter) {
    iter <- iter + 1
    old.tau <- tau
    for (i in 1:length(tau.dim)) {
      n <- tau.dim[i]
      tau.R2 <- (nmode.prod.r1(R2, tau[-n], n)
                 + premult.nmode.prod.r1(Z, EF2, tau[-n], n)
                 - premult.nmode.prod.r1(Z, lowrank.square(EF), tau[-n], n))
      tau[[i]] <- kron.nonmissing[[i]] / tau.R2
    }
    max.chg <- calc.max.abs.chg(tau, old.tau) / max(unlist(tau))
  }

  return(tau)
}

estimate.noisy.tau <- function(flash, factor = NULL) {
  tau.dim <- get.est.tau.dim(flash)

  Z   <- get.nonmissing(flash)
  EF  <- get.new.EF(flash, factor)
  EF2 <- get.new.EF2(flash, factor)

  R2 <- get.new.Rsquared(flash, factor, EF, set.missing.to.zero = FALSE)
  R2 <- R2 + lowrank.expand(EF2) - lowrank.expand(lowrank.square(EF))
  S2 <- get.given.S2(flash)

  # TODO: handle vector-valued S2.
  R2.slices <- slice.data(R2, tau.dim, Z)
  S2.slices <- slice.data(S2, tau.dim, Z)

  est.S2 <- mapply(optimize.noisy, R2.slices, S2.slices)
  each <- prod(c(1, get.dims(flash)[1:get.dim(flash) < tau.dim]))

  tau <- Z / (S2 + rep(est.S2, each = each))
  sum.tau.R2 <- sum(tau * R2)
  return(list(tau = tau, sum.tau.R2 = sum.tau.R2))
}

estimate.noisy.kron.tau <- function(flash, factor = NULL) {
  # tol and maxiter are again hardcoded.
  tol <- sqrt(.Machine$double.eps)
  maxiter <- 100

  est.S2 <- get.est.S2(factor)
  if (is.null(est.S2))
    est.S2 <- get.est.S2(flash)
  est.S2 <- as.r1.tau(est.S2, flash)
  S2.dim <- get.est.tau.dim(flash)

  Z   <- get.nonmissing(flash)
  EF  <- get.new.EF(flash, factor)
  EF2 <- get.new.EF2(flash, factor)

  # TODO: handle vector-valued S2.
  R2 <- get.new.Rsquared(flash, factor, EF, set.missing.to.zero = FALSE)
  R2 <- R2 + lowrank.expand(EF2) - lowrank.expand(lowrank.square(EF))
  S2 <- get.given.S2(flash)

  max.chg <- Inf
  iter <- 0
  while (max.chg > tol && iter < maxiter) {
    iter <- iter + 1
    old.est.S2 <- est.S2
    for (i in 1:length(S2.dim)) {
      n <- S2.dim[i]
      R2.slices <- slice.data(R2, n, Z)
      S2.slices <- slice.data(S2, n, Z)
      wts <- r1.expand(est.S2[-n])
      est.S2[[i]] <- mapply(R2.slices, S2.slices,
                            FUN = function(r2, s2) optimize.noisy(r2, s2, wts))
    }
    max.chg <- calc.max.abs.chg(est.S2, old.est.S2) / max(unlist(est.S2))
  }

  tau <- Z / (S2 + r1.expand(est.S2))
  sum.tau.R2 <- sum(tau * R2)
  return(list(est.S2 = est.S2, tau = tau, sum.tau.R2 = sum.tau.R2))
}

slice.data <- function(data, dim, Z) {
  if (all(dim == 0)) {
    if (identical(Z, 1)) {
      slices <- list(data)
    } else {
      slices <- list(data[as.logical(Z)])
    }
  } else {
    if (identical(Z, 1)) {
      slices <- lapply(apply(data, dim, list), unlist)
    } else {
      slices <- data
      is.na(slices) <- !Z
      slices <- lapply(apply(slices, dim, list.with.no.NAs), unlist)
    }
  }
  return(slices)
}

list.with.no.NAs <- function(x) list(x[!is.na(x)])

#' @importFrom stats optimize
optimize.noisy <- function(R2, S2, wts = 1) {
  interval.max <- max((R2 - S2) / wts)
  if (interval.max <= 0)
    return(0)
  opt.res <- optimize(function(x) {
    sum(log(wts * x + S2)) + sum(R2 / (wts * x + S2))
  }, interval = c(0, interval.max), tol = sqrt(.Machine$double.eps))
  return(opt.res$minimum)
}

# Helper functions to deal with various possible storage modes ----------------

get.tau.for.ebnm.calc <- function(flash, tau = NULL) {
  if (is.null(tau))
    tau <- get.tau(flash)
  if (inherits(tau, "r1"))
    return(r1.to.lowrank(tau, flash))
  if (!is.tau.simple(flash))
    return(tau)
  return(as.lowrank.tau(tau, flash))
}

as.lowrank.tau <- function(tau, flash) {
  tau.n <- get.est.tau.dim(flash)
  tau.lowrank <- lapply(as.list(get.dims(flash)),
                        function(dim) {matrix(1, nrow = dim, ncol = 1)})
  if (is.null(tau.n) || (length(tau.n) == 1)) {
    n <- get.R2.n(flash)
    tau.lowrank[[n]] <- tau * tau.lowrank[[n]]
  } else {
    tau.lowrank[tau.n] <- mapply(`*`, tau, tau.lowrank[tau.n])
  }

  class(tau.lowrank) <- "lowrank"
  return(tau.lowrank)
}

as.r1.tau <- function(tau, flash, n = NULL) {
  r1 <- as.list(rep(1, get.dim(flash)))
  r1[get.est.tau.dim(flash)] <- tau
  if (!is.null(n))
    r1 <- r1[-n]

  class(r1) <- "r1"
  return(r1)
}
