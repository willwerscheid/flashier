set.default.tol <- function(flash) {
  flash <- get.fit(flash)
  return(sqrt(.Machine$double.eps) * prod(get.dims(flash)))
}

get.conv.crit <- function(update.info) {
  return(update.info[[length(update.info)]])
}

#' Calculate the difference in ELBO
#'
#' The default objective function used to determine convergence when fitting
#'   a \code{\link{flash}} object. Calculates the difference in the
#'   variational lower bound from one iteration to the next.
#'
#' @param new The \code{flash.fit} object from the current iteration.
#'
#' @param old The \code{flash.fit} object from the previous iteration.
#'
#' @param k Ignored.
#'
#' @seealso \code{\link{conv.crit.loadings}}, \code{\link{conv.crit.factors}}
#'
#' @export
#'
conv.crit.elbo <- function(new, old, k) {
  return(calc.obj.diff(new, old, k))
}

#' Calculate the maximum absolute difference in scaled loadings
#'
#' An alternative objective function that can be used to determine
#'   convergence when fitting a \code{\link{flash}} object. Calculates the
#'   maximum absolute difference in the L2-scaled loadings,
#'   \eqn{\max_{i, k} | \frac{\ell_ik^{new}}{\| \ell_k^{new} \|_2} -
#'   \frac{\ell_ik^{old}}{\| \ell_k^{old} \|_2} |}.
#'
#' @inheritParams conv.crit.elbo
#'
#' @seealso \code{\link{conv.crit.elbo}}, \code{\link{conv.crit.factors}}
#'
#' @export
#'
conv.crit.loadings <- function(new, old, k) {
  return(calc.max.abs.chg.EF(new, old, k, n = 1))
}

#' Calculate the maximum absolute difference in scaled factors
#'
#' An alternative objective function that can be used to determine
#'   convergence when fitting a \code{\link{flash}} object. Calculates the
#'   maximum absolute difference in the L2-scaled factors,
#'   \eqn{\max_{j, k} | \frac{f_jk^{new}}{\| f_jk^{new} \|_2} -
#'   \frac{f_jk^{old}}{\| f_jk^{old} \|_2} |}.
#'
#' @inheritParams conv.crit.elbo
#'
#' @seealso \code{\link{conv.crit.elbo}}, \code{\link{conv.crit.loadings}}
#'
#' @export
#'
conv.crit.factors <- function(new, old, k) {
  return(calc.max.abs.chg.EF(new, old, k, n = 2))
}
