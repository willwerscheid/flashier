#' Initialize flash object
#'
#' Sets up a \code{\link{flash}} object with no factors. Since all other
#' \code{flash_xxx} functions take a \code{flash} or \code{flash_fit} object
#' as their first argument, calling \code{flash_init} should be the first step
#' in any \code{flash} pipeline. See \code{\link{flash}} for examples of usage.
#'
#' @inheritParams flash
#'
#' @param S_dim If the argument to \code{S} is a vector and the data matrix is
#'   square, then \code{S_dim} must specify whether \code{S} encodes row-wise or
#'   column-wise standard errors. More precisely,
#'   if \code{S_dim = 1}, then \code{S} will be interpreted as giving
#'   standard errors that vary across rows but are constant within any particular
#'   row; if \code{S_dim = 2}, then it will be interpreted as giving
#'   standard errors that vary across columns but are constant within any
#'   particular column. If \code{S} is a matrix or scalar, or if the data
#'   matrix is not square, then \code{S_dim} should be left unspecified
#'   (\code{NULL}).
#'
#' @return An initialized \code{\link{flash}} object (with no factors).
#'
#' @export
#'
flash_init <- function(data, S = NULL, var_type = 0L, S_dim = NULL) {
  flash <- set.flash.data(data, S = S, S.dim = S_dim, var.type = var_type)

  if (is.var.type.zero(flash) && !is.tau.simple(flash)) {
    flash$R <- flash$Y
  }

  # Precomputations.
  if (is.tau.simple(flash)) {
    flash$n.nonmissing <- init.n.nonmissing(flash, get.R2.n(flash))
  } else if (is.var.type.zero(flash)) {
    flash$log.2pi.s2 <- init.log.2pi.s2(get.given.tau(flash))
  } else if (is.var.type.kronecker(flash)) {
    flash$kron.nonmissing <- init.kron.nonmissing(flash)
  }

  # Calculate initial residual variances and ELBO.
  flash <- init.tau(flash)
  flash$obj <- calc.obj(flash)

  # Fields used to fix factors.
  flash$fix.dim <- list()
  flash$fix.idx <- list()

  # Some 'hidden' global options.
  flash$use.fixed.to.est.g <- FALSE

  # Fields that are no longer used.
  flash$nonmissing.thresh <- rep(0, get.dim(flash))
  flash$exclusions <- list()

  flash <- wrapup.flash(flash, output.lvl = 3L)
  flash <- flash_set_verbose(flash, verbose = 1L)

  tol <- sqrt(.Machine$double.eps) * prod(get.dims(flash_fit(flash)))
  flash <- flash_set_conv_crit(flash, flash_conv_crit_elbo_diff, tol)

  return(flash)
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
