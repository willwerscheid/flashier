#' Extract a flash_fit object
#'
#' \code{flash_fit} objects are the "internal" objects used by \code{flash}
#'   functions to fit an EBMF model. Whereas \code{flash} objects
#'   (the end results of the fitting process) include user-friendly fields and
#'   methods, \code{flash_fit} objects were not designed for public
#'   consumption and can be unwieldy. Nonetheless, some advanced
#'   \code{flash} functionality requires the wielding of
#'   \code{flash_fit} objects. In particular, initialization, convergence,
#'   and "verbose" display functions all take one or more \code{flash_fit}
#'   objects as input (see parameter \code{init_fn} in function
#'   \code{\link{flash_greedy}}; parameter \code{fn} in
#'   \code{\link{flash_set_conv_crit}};
#'   and parameter \code{fns} in \code{\link{flash_set_verbose}}).
#'   For users who would like to write custom functions, the accessor functions
#'   and methods enumerated below may prove useful. See
#'   \code{\link{flash_set_verbose}} for an example.
#'
#' The following S3 methods are available for \code{flash_fit} objects at all
#'   times except while optimizing new factor/loadings pairs as part of a
#'   "greedy" fit:
#'   \describe{
#'     \item{\code{\link{fitted.flash_fit}}}{Returns the "fitted values"
#'       \eqn{E(LF') = E(L) E(F)'}.}
#'     \item{\code{\link{residuals.flash_fit}}}{Returns the expected residuals
#'       \eqn{Y - E(LF') = Y - E(L) E(F)'}.}
#'     \item{\code{\link{ldf.flash_fit}}}{Returns an \eqn{LDF} decomposition,
#'       with columns of \eqn{L} and \eqn{F} scaled as specified by the user.}
#'   }
#'
#' @param flash A \code{flash} object.
#'
#' @return See function descriptions below.
#'
#' @export
#'
flash_fit <- function(flash) {
  return(get.fit(flash))
}


# Getters ---------------------------------------------------------------------

#' @describeIn flash_fit The posterior means for the loadings matrix \eqn{L}
#'   (when parameter \code{n} is equal to \code{1}) or factor matrix \eqn{F}
#'   (when \code{n = 2}). While optimizing new factor/loadings pairs as part of
#'   a "greedy" fit, only the posterior means for the new loadings
#'   \eqn{\ell_{\cdot k}} or factor \eqn{f_{\cdot k}} will be returned.
#' @param f A \code{flash_fit} object.
#' @param n Set \code{n = 1} to access loadings \eqn{L} and \code{n = 2} to
#'   access factors \eqn{F}).
#' @export
flash_fit_get_pm <- function(f, n) get.EF(f, n)

#' @describeIn flash_fit The posterior second moments for the loadings matrix
#'   \eqn{L} (when parameter \code{n} is equal to \code{1}) or factor matrix
#'   \eqn{F} (when \code{n = 2}). While optimizing new factor/loadings pairs,
#'   only the posterior second moments for the new loadings \eqn{\ell_{\cdot k}}
#'   or factor \eqn{f_{\cdot k}} will be returned.
#' @export
flash_fit_get_p2m <- function(f, n) get.EF2(f, n)

#' @describeIn flash_fit Equal to \eqn{1 / \sigma^2}, where \eqn{\sigma^2}
#'   is the estimated portion of the residual variance (total, by row, or
#'   by column, depending on the variance type).
#' @export
flash_fit_get_est_tau <- function(f) get.est.tau(f)

#' @describeIn flash_fit Equal to \eqn{1 / s^2}, where \eqn{s^2} is the
#'   fixed portion of the residual variance (total, by row, or by column).
#' @export
flash_fit_get_fixed_tau <- function(f) get.given.tau(f)

#' @describeIn flash_fit The overall precision \eqn{1 / (\sigma^2 + s^2)}.
#' @export
flash_fit_get_tau <- function(f) get.tau(f)

#' @describeIn flash_fit The variational lower bound (ELBO).
#' @export
flash_fit_get_elbo <- function(f) get.obj(f)

#' @describeIn flash_fit A vector containing the KL-divergence portions of
#'   the ELBO, with one value for each factor (when \code{n = 2}) or set of
#'   loadings (when \code{n = 1}). While optimizing new factor/loadings pairs,
#'   only the KL-divergence for the new factor or loadings will be returned.
#' @export
flash_fit_get_KL <- function(f, n) get.KL(f, n)

#' @describeIn flash_fit A list containing estimated priors on loadings
#'   \eqn{\hat{g}_\ell} (when \code{n = 1}) or factors \eqn{\hat{g}_f} (when
#'   \code{n = 2}). While optimizing new factor/loadings pairs, only the
#'   estimated prior on the new factor or loadings will be returned.
#' @export
flash_fit_get_g <- function(f, n) {
  if (inherits(f, "flash_fit")) {
    return(get.g.k(f, k = NULL, n = n))
  } else {
    return(get.g(f, n))
  }
}
