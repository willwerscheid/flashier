init.flash <- function(Y,
                       nonmissing = NULL,
                       F.init = NULL,
                       tau.dim = 0,
                       tau.n = 1,
                       ebnm.fn = flashr:::ebnm_pn,
                       ebnm.param = list(),
                       use.R = TRUE) {
  flash <- list()

  if (use.R) {
    flash$R <- Y - lowrank.expand(F.init)
  } else {
    flash$Y <- Y
  }

  if (is.null(nonmissing))
    nonmissing <- 1
  flash$Z <- nonmissing

  if (!is.null(F.init)) {
    flash$EF  <- F.init
    flash$EF2 <- lowrank.square(F.init)
    # TODO: calc KL, g
  }

  flash$tau.dim      <- tau.dim
  flash$tau.n        <- tau.n
  flash$n.nonmissing <- init.n.nonmissing(flash)
  flash$R2           <- init.R2(flash)
  flash$est.tau      <- calc.est.tau(flash)
  flash$obj          <- calc.obj(flash)

  flash$ebnm.fn    <- ebnm.fn
  flash$ebnm.param <- ebnm.param

  return(flash)
}

init.n.nonmissing <- function(flash) {
  Z     <- get.nonmissing(flash)
  dims  <- get.dims(flash)
  n     <- get.tau.n(flash)

  if (identical(Z, 1)) {
    n.nonmissing <- rep(prod(dims[-n]), dims[n])
  } else {
    n.nonmissing <- nmode.prod.r1(Z, r1.ones(flash), n)
  }

  # If tau is constant or zero:
  if (get.tau.dim(flash) < 1)
    n.nonmissing <- sum(n.nonmissing)

  return(n.nonmissing)
}

init.R2 <- function(flash) {
  R   <- get.R(flash)
  Y   <- get.Y(flash)
  Z   <- get.nonmissing(flash)
  EF  <- get.EF(flash)
  EF2 <- get.EF2(flash)
  n   <- get.tau.n(flash)

  if (!uses.R(flash))
    R <- Y - lowrank.expand(EF)

  R2 <- nmode.prod.r1(R^2, r1.ones(flash), n)
  if (get.n.factors(flash) > 0)
    R2 <- (R2 + premult.nmode.prod.r1(Z, EF2, r1.ones(flash), n)
           - premult.nmode.prod.r1(Z, lowrank.square(EF), r1.ones(flash), n))

  # If tau is constant or zero:
  if (get.tau.dim(flash) < 1)
    R2 <- sum(R2)

  return(R2)
}
