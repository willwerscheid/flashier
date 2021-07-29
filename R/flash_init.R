#' Initialize flash object
#'
#' Sets up a \code{flash} object with no factors. Since all other
#' \code{flash.xxx} functions take a \code{flash} or \code{flash.fit} object
#' as their first argument, calling \code{flash.init} will be the first step
#' in any \code{flash} pipeline.
#'
#' @inheritParams flash
#'
#' @param var.reg.fn The EBPM function used to regularize estimated
#'   variances. Currently only available for "simple" variance structures.
#'   The function should take two parameters called "x" and "s", which are
#'   the inputs to the EBPM problem \eqn{x_i \sim \text{Poisson}(s_i\theta_i)},
#'   \eqn{\theta_i \sim g}, and should return (at least) outputs
#'   "posterior$mean", "posterior$mean_log", and "log_likelihood".
#'
#' @param S.dim The dimension along which \code{S} lies when \code{S} is a
#'   vector. It is only necessary to specify \code{S.dim} when it cannot be
#'   inferred from the data (when, for example, \code{data} is a square
#'   matrix).
#'
#' @return A \code{flash.fit} object.
#'
#' @export
#'
flash.init <- function(data,
                       S = NULL,
                       var.type = 0L,
                       var.reg.fn = NULL,
                       S.dim = NULL) {
  flash <- set.flash.data(data, S = S, S.dim = S.dim, var.type = var.type)

  if (is.tau.simple(flash) && is.null(S)) {
    flash$var.reg.fn <- var.reg.fn
  } else if (!is.null(var.reg.fn)) {
    stop("Parameter 'var.reg.fn' can only be used with simple variance ",
         "structures.")
  }

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
  flash <- flash.set.verbose(flash, verbose = 1L)

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
