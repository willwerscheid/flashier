#' Extract a flash.fit object
#'
#' \code{flash.fit} objects are the "internal" objects used by \code{flash}
#'   functions to fit an EBMF model. Whereas \code{flash} objects
#'   (the end results of the fitting process) include user-friendly fields and
#'   methods, \code{flash.fit} objects were not designed for public
#'   consumption and can be unwieldy. Nonetheless, some advanced
#'   \code{flash} functionality requires the wielding of
#'   \code{flash.fit} objects. In particular, initialization, convergence,
#'   and "verbose" display functions all take one or more \code{flash.fit}
#'   objects as input (see parameter \code{init.fn} in function
#'   \code{\link{flash.add.greedy}}; parameter \code{conv.crit.fn} in
#'   \code{\link{flash.add.greedy}} and \code{\link{flash.backfit}};
#'   and parameter \code{disp.fns} in \code{\link{flash.set.verbose}}).
#'   For users who would like to write custom functions, the getter functions
#'   and methods enumerated below may prove useful.
#'
#' The following S3 methods are available for \code{flash.fit} objects at all
#'   times except while optimizing new factor/loadings pairs as part of a
#'   "greedy" fit:
#'   \describe{
#'     \item{\code{\link{fitted.flash.fit}}}{Returns the "fitted values"
#'       \eqn{E(LF') = E(L) E(F)'}.}
#'     \item{\code{\link{residuals.flash.fit}}}{Returns the expected residuals
#'       \eqn{Y - E(LF') = Y - E(L) E(F)'}.}
#'     \item{\code{\link{ldf.flash.fit}}}{Returns an \eqn{LDF} decomposition,
#'       with columns of \eqn{L} and \eqn{F} scaled as specified by the user.}
#'   }
#'
#' @param flash A \code{flash} object.
#'
#' @export
#'
flash.fit <- function(flash) {
  return(get.fit(flash))
}


# Getters ---------------------------------------------------------------------

#' @describeIn flash.fit The posterior means for the loadings matrix \eqn{L}
#'   (when parameter \code{n} is equal to \code{1}) or factor matrix \eqn{F}
#'   (when \code{n = 2}). While optimizing new factor/loadings pairs as part of
#'   a "greedy" fit, only the posterior means for the new loadings \eqn{\ell_k}
#'   or factor \eqn{f_k} will be returned.
#' @param f A \code{flash.fit} object.
#' @param n The dimension (set \code{n = 1} for loadings \eqn{L} and \code{n = 2}
#'   for factors \eqn{F}).
#' @export
ff.pm <- function(f, n) get.EF(f, n)

#' @describeIn flash.fit The posterior second moments for the loadings matrix
#'   \eqn{L} (when parameter \code{n} is equal to \code{1}) or factor matrix
#'   \eqn{F} (when \code{n = 2}). While optimizing new factor/loadings pairs,
#'   only the posterior second moments for the new loadings \eqn{\ell_k} or
#'   factor \eqn{f_k} will be returned.
#' @export
ff.p2m <- function(f, n) get.EF2(f, n)

#' @describeIn flash.fit Equal to \eqn{1 / \sigma^2}, where \eqn{\sigma^2}
#'   is the estimated portion of the residual variance (total, by row, or
#'   by column, depending on the variance type).
#' @export
ff.est.tau <- function(f) get.est.tau(f)

#' @describeIn flash.fit Equal to \eqn{1 / s^2}, where \eqn{s^2} is the
#'   fixed portion of the residual variance (total, by row, or by column).
#' @export
ff.fixed.tau <- function(f) get.given.tau(f)

#' @describeIn flash.fit The overall precision \eqn{1 / (\sigma^2 + s^2)}.
#' @export
ff.tau <- function(f) get.tau(f)

#' @describeIn flash.fit The variational lower bound (ELBO).
#' @export
ff.elbo <- function(f) get.obj(f)

#' @describeIn flash.fit A vector containing the KL-divergence portions of
#'   the ELBO, with one value for each factor (when \code{n = 2}) or set of
#'   loadings (when \code{n = 1}). While optimizing new factor/loadings pairs,
#'   only the KL-divergence for the new loadings or factor will be returned.
#' @export
ff.KL <- function(f, n) get.KL(f, n)

#' @describeIn flash.fit A list containing estimated priors on loadings
#'   \eqn{\hat{g}_\ell} (when \code{n = 1}) or on factors \eqn{\hat{g}_f} (when
#'   \code{n = 2}). While optimizing new factor/loadings pairs, only the
#'   estimated prior on the new loadings or factor will be returned.
#' @export
ff.g <- function(f, n) {
  if (inherits(f, "flash.fit")) {
    return(get.g.k(f, k = NULL, n = n))
  } else {
    return(get.g(f, n))
  }
}
