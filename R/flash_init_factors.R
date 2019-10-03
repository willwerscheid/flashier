# TODO: document and export.
flash.init.factors <- function(flash, EF, EF2 = NULL) {
  # TODO: allow elements of EF to be NA (need to write init fns!)

  if (is.udv(EF)) {
    EF <- udv.to.lowrank(EF)
  }

  if (!inherits(EF, "lowrank")) {
    stop("EF must be an SVD-like object or an object of class 'lowrank'.")
  }
  dims.must.match(EF, get.Y(flash))

  if (is.null(EF2)) {
    EF2 <- lowrank.square(EF)
  } else {
    if (!inherits(EF2, "lowrank")) {
      stop("If supplied, EF2 must be an object of class 'lowrank'.")
    }
    dims.must.match(EF2, get.Y(flash))
  }

  flash <- set.EF(flash, lowranks.combine(get.EF(flash), EF))
  flash <- set.EF2(flash, lowranks.combine(get.EF2(flash), EF2))

  if (uses.R(flash)) {
    R <- get.Y(flash) - lowrank.expand(get.EF(flash))
    flash <- set.R(flash, get.nonmissing(flash) * R)
  }

  flash <- init.tau(flash)
  flash <- set.obj(flash, calc.obj(flash))

  K <- ncol(EF[[1]])

  # Initialize KL at zero and g at NULL.
  for (n in 1:get.dim(flash)) {
    flash <- set.KL(flash, c(get.KL(flash, n), rep(0, K)), n)
  }
  EF.g <- rep(list(rep(list(NULL), get.dim(flash))), K)
  flash <- set.g(flash, c(get.g(flash), EF.g))

  # Initialize is.valid and is.zero.
  flash <- set.is.valid(flash, c(is.valid(flash), rep(FALSE, K)))
  flash <- set.is.zero(flash, c(is.zero(flash), rep(FALSE, K)))

  # Note that ebnm.fn doesn't get set here. This will need to be done during
  #   backfitting.

  return(flash)
}

#flash$warmstart.backfits <- warmstart.backfits
#flash$use.fixed.to.est.g <- use.fixed.to.est.g
