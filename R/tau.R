update.tau <- function(factor, flash) {
  if (is.tau.simple(flash)) {
    factor <- update.simple.tau(factor, flash)
  } else if (is.var.type.zero(flash)) {
    factor <- update.zero.tau(factor, flash)
  } else if (is.var.type.kronecker(flash)) {
    factor <- update.kronecker.tau(factor, flash)
  }
  return(factor)
}

update.simple.tau <- function(factor, flash) {
  delta.R2 <- calc.delta.R2.for.simple.tau(factor, flash)
  est.tau  <- estimate.simple.tau(flash, delta.R2)
  tau      <- tau.from.given.and.est(flash, est.tau)
  factor   <- set.delta.R2(factor, delta.R2)
  factor   <- set.est.tau(factor, est.tau)
  factor   <- set.tau(factor, tau)

  return(factor)
}

update.zero.tau <- function(factor, flash) {
  factor <- set.R(factor, calc.residuals(flash, factor))
  factor <- set.tau(factor, get.given.tau(flash))

  return(factor)
}

update.kronecker.tau <- function(factor, flash) {
  factor <- set.tau(factor, estimate.kronecker.tau(flash, factor))

  return(factor)
}

estimate.simple.tau <- function(flash, delta.R2 = 0) {
  return(get.n.nonmissing(flash) / (get.R2(flash) + delta.R2))
}

tau.from.given.and.est <- function(flash, est.tau) {
  given.tau <- get.given.tau(flash)

  # All variance is pre-specified ("zero" variance type):
  if (is.var.type.zero(flash)) {
    tau <- given.tau
  # All variance is estimated:
  } else if (is.null(given.tau)) {
    tau <- est.tau
  # Otherwise both types of variance are used. In the simple case here, the
  #   total variance is estimated and the pre-specified variance functions
  #   as a minimum:
  } else {
    tau <- pmin(given.tau, est.tau)
  }

  return(tau)
}

estimate.kronecker.tau <- function(flash, factor = NULL) {
  # Tol and maxiter are hardcoded for now.
  tol     <- 1e-3
  maxiter <- 100

  tau <- get.tau(factor)
  if (is.null(tau))
    tau <- get.tau(flash)
  tau <- as.r1.tau(tau, flash)
  tau.dim <- get.est.tau.dim(flash)
  kron.nonmissing <- get.kron.nonmissing(flash)

  EF  <- get.EF(flash)
  EF2 <- get.EF2(flash)
  if (!is.new(factor)) {
    EF  <- lowrank.drop.k(EF, get.k(factor))
    EF2 <- lowrank.drop.k(EF2, get.k(factor))
  }
  EF  <- lowranks.combine(EF, as.lowrank(get.EF(factor)))
  EF2 <- lowranks.combine(EF2, as.lowrank(get.EF2(factor)))
  EFsquared <- lowrank.square(EF)

  if (uses.R(flash) && !is.null(factor)) {
    R2 <- get.R(factor)^2
  } else if (uses.R(flash)) {
    R2 <- get.R(flash)^2
  } else {
    R2 <- (get.Y(flash) - lowrank.expand(EF))^2
  }

  max.chg <- Inf
  iter <- 0
  while (max.chg > tol && iter < maxiter) {
    iter <- iter + 1
    old.tau <- tau
    for (i in 1:length(tau.dim)) {
      n <- tau.dim[i]
      tau.R2 <- (nmode.prod.r1(R2, tau[-n], n)
                 + nmode.prod.r1(EF2, tau[-n], n)
                 - nmode.prod.r1(EFsquared, tau[-n], n))
      tau[[i]] <- kron.nonmissing[[i]] / tau.R2
    }
    max.chg <- calc.max.abs.chg(tau, old.tau)
  }

  return(tau)
}

get.tau.for.ebnm.calc <- function(flash, tau = NULL) {
  if (is.null(tau))
    tau <- get.tau(flash)
  if (is(tau, "r1"))
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

as.r1.tau <- function(tau, flash) {
  r1 <- as.list(rep(1, get.dim(flash)))
  r1[get.est.tau.dim(flash)] <- tau

  class(r1) <- "r1"
  return(r1)
}
