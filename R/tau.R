update.tau <- function(factor, flash) {
  if (is.tau.lowrank(flash)) {
    factor <- update.lowrank.tau(factor, flash)
  } else {
    factor <- update.fullrank.tau(factor, flash)
  }
  return(factor)
}

update.lowrank.tau <- function(factor, flash) {
  delta.R2 <- calc.delta.R2.for.lowrank.tau(factor, flash)
  est.tau  <- estimate.lowrank.tau(flash, delta.R2)
  tau      <- reconcile.given.and.est.tau(flash, est.tau)
  factor   <- set.delta.R2(factor, delta.R2)
  factor   <- set.est.tau(factor, est.tau)
  factor   <- set.tau(factor, tau)

  return(factor)
}

update.fullrank.tau <- function(factor, flash) {
  factor <- set.R(factor, calc.residuals(flash, factor))
  tau    <- estimate.fullrank.tau(flash, factor)
  factor <- set.tau(factor, tau)

  return(factor)
}

estimate.lowrank.tau <- function(flash, delta.R2 = 0) {
  # TODO: alternating algorithm for "Kronecker" variance type

  return(get.n.nonmissing(flash) / (get.R2(flash) + delta.R2))
}

reconcile.given.and.est.tau <- function(flash, est.tau = NULL) {
  given.tau <- get.given.tau(flash)
  if (is.null(est.tau))
    est.tau <- get.est.tau(flash)

  # All variance is pre-specified ("zero" variance type):
  if (is.tau.zero(flash)) {
    tau <- given.tau
  # All variance is estimated:
  } else if (is.null(given.tau)) {
    tau <- est.tau
  # Otherwise both types of variance are used. In the low-rank case here, the
  #   total variance is estimated and the pre-specified variance functions
  #   as a minimum:
  } else {
    tau <- pmin(given.tau, est.tau)
  }

  return(tau)
}

estimate.fullrank.tau <- function(flash, factor = NULL) {
  if (is.tau.zero(flash))
    return(get.given.tau(flash))

  # TODO: estimation when est.tau.dim is not NULL
}

get.tau.for.ebnm.calc <- function(flash, tau = NULL) {
  if (is.null(tau))
    tau <- get.tau(flash)
  if (!is.tau.lowrank(flash))
    return(tau)

  n <- get.R2.n(flash)

  tau.lowrank <- lapply(as.list(get.dims(flash)),
                        function(dim) {matrix(1, nrow = dim, ncol = 1)})
  tau.lowrank[[n]] <- tau * tau.lowrank[[n]]
  class(tau.lowrank) <- "lowrank"

  return(tau.lowrank)
}
