#' Remove factors from a flash object
#'
#' Sets flash factors to zero and then removes them from the flash object.
#'   Note that this will change the indices of existing flash factors.
#'
#' @param flash A \code{flash} or \code{flash.fit} object.
#'
#' @param kset A vector of integers specifying which factors to set to zero.
#'
#' @seealso \code{\link{flash.set.factors.to.zero}}
#'
#' @export
#'
flash.remove.factors <- function(flash, kset) {
  flash <- flash.set.factors.to.zero(flash, kset)
  flash <- get.fit(flash)

  flash <- set.EF(flash, lowrank.drop.k(get.EF(flash), kset))
  flash <- set.EF2(flash, lowrank.drop.k(get.EF2(flash), kset))
  flash <- set.KL(flash, lapply(get.KL(flash), `[`, -kset))
  flash <- set.g(flash, get.g(flash)[-kset])
  flash <- set.ebnm.fn(flash, get.ebnm.fn(flash)[-kset])
  flash <- set.is.zero(flash, is.zero(flash)[-kset])
  flash <- set.is.valid(flash, is.valid(flash)[-kset])

  fix.dim <- get.fix.dim(flash)
  fix.idx <- get.fix.idx(flash)
  fix.kset <- intersect(1:length(fix.dim), kset)

  if (length(fix.kset) > 0) {
    flash <- set.fix.idx(flash, fix.idx[-fix.kset])
    flash <- set.fix.dim(flash, fix.dim[-fix.kset])
  }

  flash <- wrapup.flash(flash, output.lvl = 3L)

  return(flash)
}

#' Set flash factors to zero
#'
#' Sets flash factors to zero but does not remove them (so as to keep the
#'   indices of existing factors the same).
#'
#' @param flash A \code{flash} or \code{flash.fit} object.
#'
#' @param kset A vector of integers specifying which factors to remove.
#'
#' @seealso \code{\link{flash.remove.factors}}
#'
#' @export
#
flash.set.factors.to.zero <- function(flash, kset) {
  flash <- get.fit(flash)
  must.be.valid.kset(flash, kset)

  for (k in kset) {
    flash <- nullcheck.factor(flash, k, verbose.lvl = 0, tol = Inf)
  }

  flash <- wrapup.flash(flash, output.lvl = 3L)

  return(flash)
}
