#' Remove factors from a flash object
#'
#' Sets factor/loadings pairs to zero and then removes them from the
#'   \code{\link{flash}} object. Note that this will change the indices of
#'   existing pairs.
#'
#' @param flash A \code{flash} or \code{flash_fit} object.
#'
#' @param kset A vector of integers specifying which factor/loadings pairs to
#'   remove.
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with the
#'   factors specified by \code{kset} removed.
#'
#' @seealso \code{\link{flash_factors_set_to_zero}}
#'
#' @export
#'
flash_factors_remove <- function(flash, kset) {
  flash <- flash_factors_set_to_zero(flash, kset)
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
#' Sets factor/loadings pairs to zero but does not remove them from the
#'   \code{\link{flash}} object (so as to keep the indices of existing pairs
#'   the same).
#'
#' @param flash A \code{flash} or \code{flash_fit} object.
#'
#' @param kset A vector of integers specifying which factor/loadings pairs to
#'   set to zero.
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with the
#'   factors specified by \code{kset} set to zero.
#'
#' @seealso \code{\link{flash_factors_remove}}
#'
#' @export
#
flash_factors_set_to_zero <- function(flash, kset) {
  flash <- get.fit(flash)
  must.be.valid.kset(flash, kset)

  for (k in kset) {
    flash <- nullcheck.factor(flash, k, verbose.lvl = 0, tol = Inf)
  }

  flash <- wrapup.flash(flash, output.lvl = 3L)

  return(flash)
}
