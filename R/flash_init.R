init.flash <- function(flash.init,
                       data,
                       EF.init,
                       est.tau.dim,
                       dim.signs,
                       ebnm.fn,
                       ebnm.param,
                       warmstart.backfits,
                       fix.dim,
                       fix.idx,
                       fix.vals,
                       use.fixed.to.est.g,
                       use.R) {
  if (is.null(flash.init)) {
    flash <- list()
  } else {
    flash <- flash.init
  }

  if (!is.null(EF.init)) {
    flash$EF  <- lowranks.combine(flash$EF, EF.init)
    flash$EF2 <- lowranks.combine(flash$EF2, lowrank.square(EF.init))
  }

  flash$est.tau.dim <- est.tau.dim

  if (bypass.init(flash)) {
    flash <- clear.bypass.init.flag(flash)
  } else {
    Y <- get.Y(data)
    nonmissing <- get.nonmissing(data)

    if (use.R) {
      flash$R <- nonmissing * (Y - lowrank.expand(flash$EF))
    } else {
      flash$Y <- Y
    }

    if (is.null(nonmissing))
      nonmissing <- 1
    flash$Z <- nonmissing

    flash$given.S2      <- get.given.S2(data)
    flash$given.tau     <- get.given.tau(data)
    flash$given.tau.dim <- get.given.tau.dim(data)

    # Precomputations.
    if (is.tau.simple(flash)) {
      flash$n.nonmissing <- init.n.nonmissing(flash, get.R2.n(flash))
    } else if (is.var.type.zero(flash)) {
      flash$log.2pi.s2 <- init.log.2pi.s2(get.given.tau(data))
    } else if (is.var.type.kronecker(flash)) {
      flash$kron.nonmissing <- init.kron.nonmissing(flash)
    }
    flash <- init.tau(flash)
    flash$obj <- calc.obj(flash)
  }

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

  flash$dim.signs  <- dim.signs
  flash$ebnm.fn    <- ebnm.fn
  flash$ebnm.param <- ebnm.param

  flash <- extend.ebnm.lists(flash)

  flash$fix.dim  <- fix.dim
  flash$fix.idx  <- fix.idx
  flash$fix.vals <- fix.vals

  flash$warmstart.backfits <- warmstart.backfits
  flash$use.fixed.to.est.g <- use.fixed.to.est.g

  return(flash)
}

init.tau.at.one <- function(flash) {
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
