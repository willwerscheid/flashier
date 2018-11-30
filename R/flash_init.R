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
  given.S2      <- get.given.S2(data)
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

  flash$est.tau.dim   <- est.tau.dim
  flash$given.S2      <- given.S2
  flash$given.tau     <- given.tau
  flash$given.tau.dim <- given.tau.dim

  # Precomputations.
  if (is.tau.simple(flash)) {
    flash$n.nonmissing    <- init.n.nonmissing(flash, get.R2.n(flash))
  } else if (is.var.type.zero(flash)) {
    flash$log.2pi.s2      <- init.log.2pi.s2(given.tau)
  } else if (is.var.type.kronecker(flash)) {
    flash$kron.nonmissing <- init.kron.nonmissing(flash)
  }
  flash <- estimate.tau(flash)

  flash$obj                <- calc.obj(flash)
  flash$dim.signs          <- dim.signs
  flash$ebnm.fn            <- ebnm.fn
  flash$ebnm.param         <- ebnm.param
  flash$warmstart.backfits <- warmstart.backfits
  flash$fix.dim            <- fix.dim
  flash$fix.idx            <- fix.idx
  flash$fix.vals           <- fix.vals
  flash$use.fixed.to.est.g <- use.fixed.to.est.g

  if (!is.null(EF.init)) {
    EF.init.k <- ncol(EF.init[[1]])

    EF.init.KL <- rep(list(rep(0, EF.init.k)), get.dim(flash))
    if (!is.null(flash$KL)) {
      flash$KL <- mapply(c, flash$KL, EF.init.KL)
    } else {
      flash$KL <- EF.init.KL
    }

    flash$g <- c(flash$g, rep(list(rep(list(NULL), get.dim(flash))), EF.init.k))
    flash$is.valid <- c(flash$is.valid, rep(FALSE, EF.init.k))
    flash$is.zero  <- c(flash$is.zero, rep(FALSE, EF.init.k))
  }

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

init.log.2pi.s2 <- function(tau) {
  return(sum(log(2 * pi / tau)))
}
