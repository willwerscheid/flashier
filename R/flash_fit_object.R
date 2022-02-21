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
#' The following S3 methods are available for \code{flash.fit} objects:
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

#' @describeIn flash.fit The matrix of posterior means for the loadings matrix
#'   \eqn{L} (when parameter \code{n} is equal to \code{1}) or factors matrix
#'   \eqn{F} (when \code{n = 2}).
#' @param f A \code{flash.fit} object.
#' @param n The "mode" (set \code{n = 1} for loadings \eqn{L} and \code{n = 2}
#'   for factors \eqn{F}).
#' @param k Only used during sequential backfits (that is, calls to
#'   \code{flash.backfit} where \code{extrapolate = FALSE}). It then gives the
#'   index of the factor/loadings pair currently being optimized.
#' @export
ff.pm <- function(f, n) get.EF(f, n)

#' @describeIn flash.fit The matrix of posterior second moments for the
#'   loadings matrix \eqn{L} (when \code{n = 1}) or factors matrix \eqn{F}
#'   (when \code{n = 2}).
#' @export
ff.p2m <- function(f, n) get.EF2(f, n)

#' @describeIn flash.fit Describes the structure of the estimated portion of
#'   the residual variance. Identical to parameter \code{var.type} in function
#'   \code{\link{flash}}.
#' @export
ff.var.type <- function(f) get.est.tau.dim(f)

#' @describeIn flash.fit Describes the structure of the fixed portion of
#'   the residual variance. Analogous to \code{ff.var.type}.
#' @export
ff.fixed.var.type <- function(f) get.given.tau.dim(f)

#' @describeIn flash.fit Equal to \eqn{1 / \sigma^2}, where \eqn{\sigma^2}
#'   is the estimated portion of the residual variance (total, by row, or
#'   by column, depending on the variance type).
#' @export
ff.est.tau <- function(f) get.est.tau(f)

#' @describeIn flash.fit Equal to \eqn{1 / s^2}, where \eqn{s^2} is the
#'   fixed portion of the residual variance (total, by row, or by column).
#' @export
ff.fixed.tau <- function(f) get.fixed.tau(f)

#' @describeIn flash.fit The overall precision \eqn{1 / (\sigma^2 + s^2)}.
#' @export
ff.tau <- function(f) get.tau(f)

#' @describeIn flash.fit The variational lower bound (ELBO).
#' @export
ff.elbo <- function(f) get.obj(f)

#' @describeIn flash.fit The KL-divergence portion of the ELBO for the
#'   \eqn{k}th set of loadings (when \code{n = 1}) or factor (when
#'   \code{n = 2}). During extrapolated backfits, a vector will be
#'   returned with one value for each factor/loadings pair.
#' @export
ff.KL <- function(f, k, n) {
  if (inherits(f, "flash.fit")) {
    return(get.KL.k(f, k, n))
  } else {
    return(get.KL(f, n))
  }
}

#' @describeIn flash.fit The estimated prior \eqn{\hat{g}} for loadings
#'   \eqn{\ell_k} (when \code{n = 1}) or for factor \eqn{f_k} (when
#'   \code{n = 2}). During extrapolated backfits, a list of priors will be
#'   returned with one element for each factor/loadings pair.
#' @export
ff.g <- function(f, k, n) {
  if (inherits(f, "flash.fit")) {
    return(get.g.k(f, k, n))
  } else {
    return(get.g(f, n))
  }
}

#' @describeIn flash.fit The data matrix, with missing data (NAs) replaced by
#'   zeros.
#' @export
ff.data <- function(f) get.Y(f)

#' @describeIn flash.fit A matrix of 0s and 1s indicating whether data
#'   is missing (0) or present (1). If there is no missing data, then a scalar
#'   might be stored rather than a full matrix of 1s.
#' @export
ff.nonmissing <- function(f) get.nonmissing(f)

#' @describeIn flash.fit The sum of expected squared residuals (total, by row,
#'   or by column, depending on the variance type).
#' @export
ff.R2 <- function(f) get.R2(f)
