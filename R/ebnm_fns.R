#' EBNM functions
#'
#' \code{as.ebnm.fn} is a helper function that provides readable syntax for
#'   constructing \code{\link[ebnm]{ebnm}} functions that can serve as
#'   arguments to parameter \code{ebnm.fn} in functions \code{\link{flash}},
#'   \code{\link{flash.add.greedy}}, and \code{\link{flash.init.factors}} (see
#'   \strong{Examples} below). It is also possible to write a custom function
#'   from scratch: see \strong{Details} below for a simple example. A more
#'   involved example can be found in the "Advanced flashier" vignette.
#'
#' As input to parameter \code{ebnm.fn}, it should suffice for many purposes to
#'   provide functions from package \code{ebnm} as is (for example, one might
#'   set \code{ebnm.fn = ebnm::ebnm_point_laplace}). To use non-default
#'   arguments, function \code{as.ebnm.fn} may be used (see \strong{Examples}).
#'   Custom functions may also be written. In general, any function that is
#'   used as an argument to parameter \code{ebnm.fn} must accept parameters:
#'   \describe{
#'   \item{\code{x}}{A vector of observations.}
#'   \item{\code{s}}{A vector of standard errors, or a scalar if all standard
#'     errors are equal.}
#'   \item{\code{g_init}}{The prior \eqn{g}. Usually, this is left unspecified
#'     (\code{NULL}) and estimated from the data. If it is supplied and
#'     \code{fix_g = TRUE}, then the prior is fixed at \code{g_init}; if
#'     \code{fix_g = FALSE}, then \code{g_init} gives the
#'     initial value of \eqn{g} used during optimization. In \code{flashier},
#'     \eqn{g} is fixed during the wrap-up phase when estimating local false
#'     sign rates and constructing a sampler; \code{g_init} is used
#'     with \code{fix_g = FALSE} to "warmstart" backfits
#'     (see \code{\link{flash.backfit}}). If none of these features (local
#'     false sign rates, samplers, or warmstarts) are needed,
#'     then both \code{g_init} and \code{fix_g} can be ignored (the EBNM
#'     function must still accept them as parameters, but it need not do
#'     anything with their arguments).}
#'   \item{\code{fix_g}}{If \code{TRUE}, the prior is fixed at \code{g_init}
#'     instead of estimated. See the description of \code{g_init} above.}
#'   \item{\code{output}}{A character vector indicating which values are to be
#'     returned. Custom EBNM functions can safely ignore this parameter (again,
#'     it must accept it as a parameter, but it does not need to do anything
#'     with its argument).}
#'   }
#'   The return object must be a list that includes fields:
#'   \describe{
#'   \item{\code{posterior}}{A data frame that includes columns \code{mean}
#'     and \code{second_moment} (the first and second moments for each
#'     posterior distribution
#'     \eqn{p(\theta_i \mid s_i, \hat{g}), i = 1, ..., n}). Optionally,
#'     a column \code{lfsr} giving local false sign rates may also be
#'     included.}
#'   \item{\code{fitted_g}}{The estimated prior \eqn{\hat{g}}. Within
#'     \code{flashier}, \code{fitted_g} is only ever used as an argument to
#'     \code{g_init} in subsequent calls to the same EBNM function, so the
#'     manner in which it is represented is unimportant.}
#'   \item{\code{log_likelihood}}{The optimal log likelihood
#'     \eqn{L(\hat{g}) := \sum_i \log p(x_i \mid \hat{g}, s_i)}.}
#'   \item{\code{posterior_sampler}}{An optional field containing a function
#'     that samples from the posterior distributions of the "means"
#'     \eqn{\theta_i}. If included, the function should take a single parameter
#'     \code{nsamp} and return a matrix where rows correspond to samples and
#'     columns correspond to observations (that is, there should be
#'     \code{nsamp} rows and \eqn{n} columns).}
#'   }
#'
#' @param ... Parameters to be passed to function \code{\link[ebnm]{ebnm}}
#'   in package \code{ebnm}. An argument to \code{prior_family} should be
#'   provided unless the default family of point-normal priors is desired.
#'   Arguments to parameters \code{x}, \code{s}, or \code{output} must not be
#'   included. Finally, if \code{g_init} is included, then \code{fix_g = TRUE}
#'   must be as well. To fix a prior grid, use parameter \code{scale} rather
#'   than \code{g_init}.
#'
#' @seealso \code{\link[ebnm]{ebnm}}, \code{\link{flash}},
#'   \code{\link{flash.add.greedy}}, \code{\link{flash.init.factors}}.
#'
#' @examples
#' # A custom EBNM function might be written as follows:
#'
#' my.ebnm.fn <- function(x, s, g_init, fix_g, output) {
#'   ebnm.res <- ebnm::ebnm_point_laplace(
#'     x = x,
#'     s = s,
#'     g_init = g_init,
#'     fix_g = fix_g,
#'     output = output,
#'     control = list(iterlim = 10)
#'   )
#'   return(ebnm.res)
#' }
#'
#' # The following are equivalent:
#'
#' fl1 <- flash(
#'   gtex,
#'   ebnm.fn = my.ebnm.fn,
#'   greedy.Kmax = 2
#' )
#'
#' fl2 <- flash(
#'   gtex,
#'   ebnm.fn = as.ebnm.fn(
#'     prior_family = "point_laplace",
#'     control = list(iterlim = 10)
#'   ),
#'   greedy.Kmax = 2
#' )
#'
#' @importFrom ebnm ebnm
#'
#' @export
#'
as.ebnm.fn <- function(...) {
  if (any(c("x", "s", "output") %in% names(list(...)))) {
    stop("'x', 's', and 'output' must not be supplied as arguments to ",
         "'as.ebnm.fn'.")
  }

  if (any(c("g_init", "fix_g") %in% names(list(...)))) {
    args <- list(...)
    if (is.null(args$g_init)) {
      stop("If 'fix_g' is supplied as an argument to as.ebnm.fn, then ",
           "'g_init' must be as well.")
    }
    if (is.null(args$fix_g) || !args$fix_g) {
      stop("If 'g_init' is supplied as an argument to as.ebnm.fn, then it ",
           "must be fixed (i.e., 'fix_g = TRUE' should be included as an ",
           "argument). If you were trying to fix the grid of the prior, then ",
           "use parameter 'scale' rather than 'g_init'.")
    }

    ebnm.fn <- function(x, s, g_init, fix_g, output) {
      ebnm::ebnm(x, s, output = output, ...)
    }
  } else {
    ebnm.fn <- function(x, s, g_init, fix_g, output) {
      ebnm::ebnm(x, s, g_init = g_init, fix_g = fix_g, output = output, ...)
    }
  }

  return(ebnm.fn)
}

handle.ebnm.fn <- function(ebnm.fn, data.dim) {
  error.msg <- paste(
    "Argument to ebnm.fn is incorrectly specified. See ?as.ebnm.fn for details."
  )

  if (length(ebnm.fn) == 1)
    ebnm.fn <- rep(list(ebnm.fn), data.dim)
  if (length(ebnm.fn) != data.dim)
    stop(error.msg)

  for (i in 1:length(ebnm.fn))
    test.ebnm.fn(ebnm.fn[[i]], error.msg)

  return(ebnm.fn)
}

test.ebnm.fn <- function(ebnm.fn, error.msg) {
  tryCatch(
    tmp <- ebnm.fn(
      x = c(-1, 0, 1),
      s = rep(1, 3),
      g_init = NULL,
      fix_g = FALSE,
      output = default.ebnm.output
    ),
    error = function(e) stop(error.msg))
}

# Output required for usual flash updates.
#
default.ebnm.output <- c("posterior_mean",
                         "posterior_second_moment",
                         "fitted_g",
                         "log_likelihood")

# Since ebnm will typically be called many times, it suffices to do a small
#   number of optimization iterations each time.
#
mixsqp.defaults <- list(maxiter.sqp = 10, eps = 1e-15)
nlm.defaults <- list(iterlim = 10)

# Since maxiter.sqp is intentionally set to be small, ignore the ashr warning
#   about reaching the maximum number of iterations.
# Also ignore the ebnm warning about setting mode and scale when g is fixed.
#   This pops up when calculating LFSR and posterior samplers using
#   nonzero.mode prior families.
#
ignored.warnings <- c("Optimization failed to converge. Results",
                      "mode and scale parameters are ignored")

#' @importFrom ebnm ebnm
#'
ebnm.nowarn <- function(...) {
  withCallingHandlers(res <- ebnm(...),
                      warning = function(w) {
                        if (any(startsWith(conditionMessage(w), ignored.warnings)))
                          invokeRestart("muffleWarning")
                      })
  return(res)
}
