# TODO: remove defaults once tests have been removed
#
init.new.flash <- function(Y,
                           nonmissing = NULL,
                           EF.init = NULL,
                           given.tau = NULL,
                           given.tau.dim = NULL,
                           est.tau.dim = 0,
                           dim.signs = NULL,
                           fix.dim = NULL,
                           fix.idx = NULL,
                           fix.vals = NULL,
                           use.fixed.to.est.g = FALSE,
                           ebnm.fn = ebnm.pn,
                           ebnm.param = list(list(prior_type = "point_normal")),
                           warmstart.backfits = TRUE,
                           use.R = TRUE) {
  flash <- list()

  if (use.R) {
    flash$R <- Y - lowrank.expand(EF.init$EF)
  } else {
    flash$Y <- Y
  }

  if (is.null(nonmissing))
    nonmissing <- 1
  flash$Z <- nonmissing

  if (!is.null(EF.init)) {
    flash$EF  <- EF.init
    flash$EF2 <- lowrank.square(EF.init)
  } else {
    flash$EF  <- NULL
    flash$EF2 <- NULL
  }

  flash$given.tau          <- given.tau
  flash$given.tau.dim      <- given.tau.dim
  flash$est.tau.dim        <- est.tau.dim
  flash$dim.signs          <- dim.signs
  flash$fix.dim            <- fix.dim
  flash$fix.idx            <- fix.idx
  flash$fix.vals           <- fix.vals
  flash$use.fixed.to.est.g <- use.fixed.to.est.g
  flash$ebnm.fn            <- ebnm.fn
  flash$ebnm.param         <- ebnm.param
  flash$warmstart.backfits <- warmstart.backfits

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

force.use.R <- function(given.tau) {
  !(is.null(given.tau) || is.vector(given.tau))
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
    R <- Y - lowrank.expand(EF)

  R2 <- nmode.prod.r1(R^2, r1.ones(flash), n)
  if (get.n.factors(flash) > 0)
    R2 <- (R2 + premult.nmode.prod.r1(Z, EF2, r1.ones(flash), n)
           - premult.nmode.prod.r1(Z, lowrank.square(EF), r1.ones(flash), n))

  if (store.R2.as.scalar(flash))
    R2 <- sum(R2)

  return(R2)
}
