# TODO: remove defaults once tests have been removed
#
init.flash <- function(flash.init,
                       data = data,
                       EF.init = NULL,
                       est.tau.dim = 0,
                       dim.signs = NULL,
                       ebnm.fn = ebnm.pn,
                       ebnm.param = list(list(prior_type = "point_normal")),
                       warmstart.backfits = TRUE,
                       fix.dim = NULL,
                       fix.idx = NULL,
                       fix.vals = NULL,
                       use.fixed.to.est.g = FALSE,
                       use.R = TRUE) {
  flash <- flash.init
  if (is.null(flash))
    flash <- list()

  flash$EF       <- lowranks.combine(flash$EF, EF.init)
  flash$EF2      <- lowranks.combine(flash$EF2, lowrank.square(EF.init))
  flash$is.valid <- c(flash$is.valid, rep(FALSE, length(EF.init)))
  flash$is.zero  <- c(flash$is.zero, rep(FALSE, length(EF.init)))

  Y             <- get.Y(data)
  nonmissing    <- get.nonmissing(data)
  given.tau     <- get.given.tau(data)
  given.tau.dim <- get.given.tau.dim(data)

  if (use.R) {
    flash$R <- nonmissing * (Y - lowrank.expand(flash$EF))
  } else {
    flash$Y <- Y
  }

  if (is.null(nonmissing))
    nonmissing <- 1
  flash$Z <- nonmissing

  flash$est.tau.dim <- est.tau.dim
  if (!is.null(est.tau.dim) && !is.null(given.tau.dim) && est.tau.dim > 0) {
    if (given.tau.dim == 0) {
      # Convert given.tau to a vector:
      flash$given.tau.dim <- est.tau.dim
      flash$given.tau <- rep(given.tau, dim(Y)[est.tau.dim])
    } else if (given.tau.dim != est.tau.dim) {
      # Can't do efficient var estimation, so use full given.tau matrix/array:
      flash$given.tau.dim <- NULL
      dims <- c(given.tau.dim, (1:length(dim(Y)))[-given.tau.dim])
      flash$given.tau <- array(given.tau, dim = dim(Y)[dims])
      flash$given.tau <- aperm(flash$given.tau, perm = order(dims))
    } else {
      flash$given.tau.dim <- given.tau.dim
      flash$given.tau     <- given.tau
    }
  } else {
    flash$given.tau.dim <- given.tau.dim
    flash$given.tau     <- given.tau
  }
  if (is.matrix(flash$given.tau))
    flash$given.tau <- flash$Z * flash$given.tau

  flash$dim.signs          <- dim.signs
  flash$ebnm.fn            <- ebnm.fn
  flash$ebnm.param         <- ebnm.param
  flash$warmstart.backfits <- warmstart.backfits
  flash$fix.dim            <- fix.dim
  flash$fix.idx            <- fix.idx
  flash$fix.vals           <- fix.vals
  flash$use.fixed.to.est.g <- use.fixed.to.est.g

  if (is.tau.lowrank(flash)) {
    flash$n.nonmissing <- init.n.nonmissing(flash)
    flash$R2           <- init.R2(flash)
    flash$est.tau      <- estimate.lowrank.tau(flash)
    flash$tau          <- reconcile.given.and.est.tau(flash)
  } else {
    if (is.tau.zero(flash))
      flash$log.2pi.s2 <- precompute.log.2pi.s2(given.tau)
    flash$tau          <- estimate.fullrank.tau(flash)
  }
  flash$obj            <- calc.obj(flash)

  return(flash)
}

force.use.R <- function(data) {
  tau <- get.given.tau(data)
  return(!(is.null(tau) || is.vector(tau)))
}

init.n.nonmissing <- function(flash) {
  Z     <- get.nonmissing(flash)
  dims  <- get.dims(flash)
  n     <- get.R2.n(flash)

  if (identical(Z, 1)) {
    n.nonmissing <- rep(prod(dims[-n]), dims[n])
  } else {
    n.nonmissing <- nmode.prod.r1(Z, r1.ones(flash), n)
  }

  if (store.R2.as.scalar(flash))
    n.nonmissing <- sum(n.nonmissing)

  return(n.nonmissing)
}

init.R2 <- function(flash) {
  R   <- get.R(flash)
  Y   <- get.Y(flash)
  Z   <- get.nonmissing(flash)
  EF  <- get.EF(flash)
  EF2 <- get.EF2(flash)
  n   <- get.R2.n(flash)

  if (!uses.R(flash))
    R <- Z * (Y - lowrank.expand(EF))

  R2 <- nmode.prod.r1(R^2, r1.ones(flash), n)
  if (get.n.factors(flash) > 0)
    R2 <- (R2 + premult.nmode.prod.r1(Z, EF2, r1.ones(flash), n)
           - premult.nmode.prod.r1(Z, lowrank.square(EF), r1.ones(flash), n))

  if (store.R2.as.scalar(flash))
    R2 <- sum(R2)

  return(R2)
}
