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
  if (is.null(flash.init)) {
    flash <- list()
  } else {
    flash <- flash.init
  }

  Y             <- get.Y(data)
  nonmissing    <- get.nonmissing(data)
  given.tau     <- get.given.tau(data)
  given.tau.dim <- get.given.tau.dim(data)

  flash$EF  <- lowranks.combine(flash$EF, EF.init)
  flash$EF2 <- lowranks.combine(flash$EF2, lowrank.square(EF.init))
  if (use.R) {
    flash$R <- nonmissing * (Y - lowrank.expand(flash$EF))
  } else {
    flash$Y <- Y
  }

  if (is.null(nonmissing))
    nonmissing <- 1
  flash$Z <- nonmissing

  flash$est.tau.dim <- est.tau.dim
  flash$given.tau.dim <- given.tau.dim
  if (!is.null(flash$given.tau) && !is.vector(flash$given.tau)) {
    flash$given.tau <- nonmissing * given.tau
  } else {
    flash$given.tau <- given.tau
  }

  if (is.tau.simple(flash)) {
    flash$n.nonmissing    <- init.n.nonmissing(flash, get.R2.n(flash))
    flash$R2              <- init.R2(flash)
    flash$est.tau         <- estimate.simple.tau(flash)
    flash$tau             <- tau.from.given.and.est(flash, flash$est.tau)
  } else if (is.var.type.zero(flash)) {
    flash$log.2pi.s2      <- init.log.2pi.s2(given.tau)
    flash$tau             <- given.tau
  } else if (is.var.type.kronecker(flash)) {
    flash$kron.nonmissing <- init.kron.nonmissing(flash)
    flash$tau             <- init.kronecker.tau(flash)
    flash$tau             <- estimate.kronecker.tau(flash)
  } else if (is.var.type.noisy(flash)) {
    noisy.tau             <- estimate.noisy.tau(flash)
    flash$sum.tau.R2      <- noisy.tau$sum.tau.R2
    flash$tau             <- noisy.tau$tau
  } else {
    stop("The requested variance structure has not yet been implemented.")
  }

  flash$obj                <- calc.obj(flash)
  flash$dim.signs          <- dim.signs
  flash$ebnm.fn            <- ebnm.fn
  flash$ebnm.param         <- ebnm.param
  flash$warmstart.backfits <- warmstart.backfits
  flash$fix.dim            <- fix.dim
  flash$fix.idx            <- fix.idx
  flash$fix.vals           <- fix.vals
  flash$use.fixed.to.est.g <- use.fixed.to.est.g

  flash$is.valid <- c(flash$is.valid, rep(FALSE, length(EF.init)))
  flash$is.zero  <- c(flash$is.zero, rep(FALSE, length(EF.init)))

  return(flash)
}

init.kronecker.tau <- function(flash) {
  return(lapply(get.dims(flash), function(dim) rep(1, dim)))
}

# Precomputations for estimating variance and calculating objective -----------

init.n.nonmissing <- function(flash, n) {
  Z    <- get.nonmissing(flash)
  dims <- get.dims(flash)

  if (identical(Z, 1)) {
    n.nonmissing <- rep(prod(dims[-n]), dims[n])
  } else {
    n.nonmissing <- nmode.prod.r1(Z, r1.ones(flash), n)
  }

  if (store.R2.as.scalar(flash))
    n.nonmissing <- sum(n.nonmissing)

  return(n.nonmissing)
}

init.kron.nonmissing <- function(flash) {
  return(lapply(get.est.tau.dim(flash),
                function(n) {init.n.nonmissing(flash, n)}))
}

init.R2 <- function(flash) {
  R   <- get.R(flash)
  Y   <- get.Y(flash)
  Z   <- get.nonmissing(flash)
  EF  <- get.EF(flash)
  EF2 <- get.EF2(flash)

  if (!uses.R(flash))
    R <- Z * (Y - lowrank.expand(EF))

  n  <- get.R2.n(flash)
  R2 <- (nmode.prod.r1(R^2, r1.ones(flash), n)
         + premult.nmode.prod.r1(Z, EF2, r1.ones(flash), n)
         - premult.nmode.prod.r1(Z, lowrank.square(EF), r1.ones(flash), n))

  if (store.R2.as.scalar(flash))
    R2 <- sum(R2)

  return(R2)
}

init.log.2pi.s2 <- function(tau) {
  return(sum(log(2 * pi / tau)))
}
