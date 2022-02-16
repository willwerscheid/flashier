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
