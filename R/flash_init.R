init.flash <- function(flash.init,
                       data,
                       EF.init,
                       est.tau.dim,
                       dim.signs,
                       ebnm.fn,
                       warmstart.backfits,
                       fix.dim,
                       fix.idx,
                       fix.vals,
                       use.fixed.to.est.g,
                       nonmissing.thresh,
                       use.R) {
  if (is.null(flash.init)) {
    flash <- list()
  } else {
    flash <- flash.init
  }

  if (!is.null(EF.init)) {
    # Allow an svd object here.
    if (all(c("u", "d", "v") %in% names(EF.init))) {
      EF.init <- list(t(t(EF.init$u) * sqrt(EF.init$d)),
                      t(t(EF.init$v) * sqrt(EF.init$d)))
    }

    class(EF.init) <- "lowrank"
    flash$EF  <- lowranks.combine(get.EF(flash), EF.init)
    flash$EF2 <- lowranks.combine(get.EF2(flash), lowrank.square(EF.init))
  }

  flash$est.tau.dim <- est.tau.dim

  if (!bypass.init(flash) || !is.null(EF.init)) {
    flash$Y <- get.Y(data)
    nonmissing <- get.nonmissing(data)

    if (use.R) {
      flash$R <- nonmissing * (flash$Y - lowrank.expand(get.EF(flash)))
    }

    if (is.null(nonmissing))
      nonmissing <- 1
    flash$Z <- nonmissing

    flash$given.S2      <- get.given.S2(data)
    flash$given.tau     <- get.given.tau(data)
    flash$given.tau.dim <- get.given.tau.dim(data)

    if (any.missing(flash) && is.var.type.noisy.kron(flash)) {
      stop("The noisy Kronecker variance structure has not yet been implemented ",
           "for missing data.")
    }

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

  flash <- clear.bypass.init.flag(flash)

  if (!is.null(EF.init)) {
    EF.init.k <- ncol(EF.init[[1]])

    # For each factor, initialize KL at zero and g at NULL.
    EF.init.KL <- rep(list(rep(0, EF.init.k)), get.dim(flash))
    if (!is.null(get.KL(flash))) {
      flash$KL <- mapply(c, get.KL(flash), EF.init.KL)
    } else {
      flash$KL <- EF.init.KL
    }
    EF.init.g <- rep(list(rep(list(NULL), get.dim(flash))), EF.init.k)
    flash$g <- c(get.g(flash), EF.init.g)

    # Initialize is.valid, is.zero, and exclusions.
    flash$is.valid   <- c(is.valid(flash), rep(FALSE, EF.init.k))
    flash$is.zero    <- c(is.zero(flash), rep(FALSE, EF.init.k))
    flash$exclusions <- c(get.exclusions(flash),
                          rep(list(rep(list(integer(0)), get.dim(flash))),
                              EF.init.k))
  }

  flash$dim.signs  <- dim.signs
  flash$ebnm.fn    <- ebnm.fn

  flash <- extend.ebnm.lists(flash)

  flash$fix.dim  <- fix.dim
  flash$fix.idx  <- fix.idx
  flash$fix.vals <- fix.vals

  if (is.null(nonmissing.thresh))
    nonmissing.thresh <- get.default.nonmissing.thresh(flash)
  flash$nonmissing.thresh <- nonmissing.thresh

  if (is.null(flash$exclusions))
    flash$exclusions <- list()

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
                function(n) init.n.nonmissing(flash, n)))
}

init.log.2pi.s2 <- function(tau) {
  return(sum(log(2 * pi / tau[tau > 0])))
}

# Default threshold of nonmissingness required to estimate loadings -----------

get.default.nonmissing.thresh <- function(flash) {
  return(rep(0, get.dim(flash)))
}
