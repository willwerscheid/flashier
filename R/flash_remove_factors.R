# TODO: handle conv.stat better... probably remove from final object
#' @export
#
flash.remove.factors <- function(flash, kset, output.lvl = 3) {
  if (inherits(flash, "flash")) {
    conv.stat <- get.conv.stat(flash)
  }

  flash <- flash.set.factors.to.zero(flash, kset, output.lvl = 0)

  flash <- set.EF(flash, lowrank.drop.k(get.EF(flash), kset))
  flash <- set.EF2(flash, lowrank.drop.k(get.EF2(flash), kset))
  flash <- set.KL(flash, lapply(get.KL(flash), `[`, -kset))
  flash <- set.g(flash, get.g(flash)[-kset])
  flash <- set.ebnm.fn(flash, get.ebnm.fn(flash)[-kset])
  flash <- set.is.zero(flash, is.zero(flash)[-kset])
  flash <- set.is.valid(flash, is.valid(flash)[-kset])

  flash <- wrapup.flash(flash, output.lvl, is.converged = TRUE)
  flash$convergence.status <- conv.stat

  return(flash)
}

#' @export
#
flash.set.factors.to.zero <- function(flash, kset, output.lvl = 3) {
  if (inherits(flash, "flash")) {
    conv.stat <- get.conv.stat(flash)
    flash <- get.fit(flash)
  } else {
    conv.stat <- NULL
  }

  for (k in kset) {
    flash <- nullcheck.factor(flash, k, verbose.lvl = 0, tol = Inf)
  }

  flash <- wrapup.flash(flash, output.lvl, is.converged = TRUE)
  flash$convergence.status <- conv.stat

  return(flash)
}
