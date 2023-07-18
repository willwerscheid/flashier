#' Reorder factors in a flash object
#'
#' Reorders the factor/loadings pairs in a \code{\link{flash}} object.
#'
#' @param flash A \code{flash} or \code{flash_fit} object.
#'
#' @param kset A vector of integers specifying the new order of the
#'   factor/loadings pairs. All existing factors must be included in
#'   \code{kset}; to drop factors, use \code{\link{flash_factors_remove}}.
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with the
#'   factors reordered according to argument \code{kset}.
#'
#' @export
#'
flash_factors_reorder <- function(flash, kset) {
  flash <- get.fit(flash)
  must.be.valid.kset(flash, kset)

  if (length(kset) != get.n.factors(flash)) {
    stop("kset must include all factors.")
  }

  flash <- set.EF(flash, lapply(get.EF(flash), function(X) X[, kset, drop = FALSE]))
  flash <- set.EF2(flash, lapply(get.EF2(flash), function(X) X[, kset, drop = FALSE]))
  flash <- set.KL(flash, lapply(get.KL(flash), `[`, kset))
  flash <- set.g(flash, get.g(flash)[kset])
  flash <- set.ebnm.fn(flash, get.ebnm.fn(flash)[kset])
  flash <- set.is.zero(flash, is.zero(flash)[kset])
  flash <- set.is.valid(flash, is.valid(flash)[kset])

  fix.dim <- get.fix.dim(flash)
  fix.dim <- c(fix.dim, rep(list(NULL), length(kset) - length(fix.dim)))
  fix.idx <- get.fix.idx(flash)
  fix.idx <- c(fix.idx, rep(list(NULL), length(kset) - length(fix.dim)))

  flash <- set.fix.idx(flash, fix.idx[kset])
  flash <- set.fix.dim(flash, fix.dim[kset])

  flash <- wrapup.flash(flash, output.lvl = 3L)

  return(flash)
}
