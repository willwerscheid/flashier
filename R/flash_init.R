#' Initialize flash fit object
#'
#' Sets up a \code{flash.fit} object with no factors. Since all other
#' \code{flash.xxx} functions take a \code{flash.fit} object as their first
#' argument, initializing a \code{flash.fit} object will be the first step
#' in any \code{flash} pipeline.
#'
#' @inheritParams flash
#'
#' @param S.dim The dimension along which \code{S} lies when \code{S} is a
#'   vector. Only necessary when it cannot be inferred from the data (when,
#'   for example, \code{data} is a square matrix).
#'
#' @return A \code{flash.fit} object.
#'
#' @export
#'
flash.init <- function(data, var.type, S = NULL, S.dim = NULL) {
  flash <- set.flash.data(data, S = S, S.dim = S.dim, var.type = var.type)

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

  # These fields are used to fix factors.
  flash$fix.dim <- list()
  flash$fix.idx <- list()

  # These fields are no longer used.
  flash$nonmissing.thresh <- rep(0, get.dim(flash))
  flash$exclusions <- list()

  class(flash) <- c("flash.fit", "list")

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
