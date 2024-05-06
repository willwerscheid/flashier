#' Construct an EBNM function
#'
#' \code{flash_ebnm} is a helper function that provides readable syntax for
#'   constructing \code{\link[ebnm]{ebnm}} functions that can serve as
#'   arguments to parameter \code{ebnm_fn} in functions \code{\link{flash}},
#'   \code{\link{flash_greedy}}, and \code{\link{flash_factors_init}} (see
#'   \strong{Examples} below). It is also possible to write a custom function
#'   from scratch: see \strong{Details} below for a simple example. A more
#'   involved example can be found in the "Extending ebnm with custom ebnm-style
#'   functions" vignette in the \code{\link[ebnm]{ebnm}} package.
#'
#' As input to parameter \code{ebnm_fn} in functions \code{\link{flash}},
#'   \code{\link{flash_greedy}}, and \code{\link{flash_factors_init}},
#'   it should suffice for many purposes to
#'   provide functions from package \code{ebnm} as is (for example, one might
#'   set \code{ebnm_fn = ebnm_point_laplace}). To use non-default
#'   arguments, function \code{flash_ebnm} may be used (see \strong{Examples}).
#'   Custom functions may also be written; for details, see the "Extending
#'   ebnm with custom ebnm-style functions" vignette in the
#'   \code{\link[ebnm]{ebnm}} package.
#'
#' @param ... Parameters to be passed to function \code{\link[ebnm]{ebnm}}
#'   in package \code{ebnm}. An argument to \code{prior_family} should be
#'   provided unless the default family of point-normal priors is desired.
#'   Arguments to parameters \code{x}, \code{s}, or \code{output} must not be
#'   included. Finally, if \code{g_init} is included, then \code{fix_g = TRUE}
#'   must be as well. To fix a prior grid, use parameter \code{scale} rather
#'   than \code{g_init}.
#'
#' @return A function that can be passed as argument to parameter
#'   \code{ebnm_fn} in functions \code{\link{flash}},
#'   \code{\link{flash_greedy}}, and \code{\link{flash_factors_init}}.
#'
#' @seealso \code{\link[ebnm]{ebnm}}
#'
#' @examples
#' # A custom EBNM function might be written as follows:
#' my_ebnm_fn <- function(x, s, g_init, fix_g, output) {
#'   ebnm_res <- ebnm_point_laplace(
#'     x = x,
#'     s = s,
#'     g_init = g_init,
#'     fix_g = fix_g,
#'     output = output,
#'     control = list(iterlim = 10)
#'   )
#'   return(ebnm_res)
#' }
#'
#' # The following are equivalent:
#' fl1 <- flash(
#'   gtex,
#'   ebnm_fn = my_ebnm_fn,
#'   greedy_Kmax = 2
#' )
#' fl2 <- flash(
#'   gtex,
#'   ebnm_fn = flash_ebnm(
#'     prior_family = "point_laplace",
#'     control = list(iterlim = 10)
#'   ),
#'   greedy_Kmax = 2
#' )
#'
#' @importFrom ebnm ebnm ebnm_group
#'
#' @export
#'
flash_ebnm <- function(...) {
  args <- list(...)
  if (any(c("x", "s", "output") %in% names(args))) {
    stop("'x', 's', and 'output' must not be supplied as arguments to ",
         "'flash_ebnm'.")
  }

  if (any(c("g_init", "fix_g") %in% names(args))) {
    if (is.null(args$g_init)) {
      stop("If 'fix_g' is supplied as an argument to flash_ebnm, then ",
           "'g_init' must be as well.")
    }
    if (is.null(args$fix_g) || !args$fix_g) {
      stop("If 'g_init' is supplied as an argument to flash_ebnm, then it ",
           "must be fixed (i.e., 'fix_g = TRUE' should be included as an ",
           "argument). If you were trying to fix the grid of the prior, then ",
           "use parameter 'scale' rather than 'g_init'.")
    }

    if (!is.null(args$group)) {
      ebnm.fn <- function(x, s, g_init, fix_g, output) {
        # Workaround to pass test.ebnm.fn:
        if (identical(x, c(-10, 0, 10))) {
          x <- rep(x, length.out = length(args$group))
          s <- rep(s, length.out = length(args$group))
        }
        ebnm_group(x, s, output = output, ...)
      }
    } else {
      ebnm.fn <- function(x, s, g_init, fix_g, output) {
        ebnm(x, s, output = output, ...)
      }
    }
  } else {
    if (!is.null(args$group)) {
      ebnm.fn <- function(x, s, g_init, fix_g, output) {
        # Workaround to pass test.ebnm.fn:
        if (identical(x, c(-10, 0, 10))) {
          x <- rep(x, length.out = length(args$group))
          s <- rep(s, length.out = length(args$group))
        }
        ebnm_group(
          x, s, g_init = g_init, fix_g = fix_g, output = output, ...
        )
      }
    } else {
      ebnm.fn <- function(x, s, g_init, fix_g, output) {
        suppressWarnings(
          ebnm(x, s, g_init = g_init, fix_g = fix_g, output = output, ...)
        )
      }
    }
  }

  return(ebnm.fn)
}

handle.ebnm.fn <- function(ebnm.fn, data.dim) {
  error.msg <- paste(
    "Argument to ebnm.fn is incorrectly specified. See ?flash_ebnm for details."
  )

  if (length(ebnm.fn) == 1)
    ebnm.fn <- rep(list(ebnm.fn), data.dim)
  if (length(ebnm.fn) != data.dim)
    stop(error.msg)

  dim.signs <- sapply(ebnm.fn, test.ebnm.fn, error.msg)

  return(list(ebnm.fn = ebnm.fn, dim.signs = dim.signs))
}

#' @importFrom stats coef
#'
test.ebnm.fn <- function(ebnm.fn, error.msg) {
  tryCatch(
    tmp <- ebnm.fn(
      x = c(-10, 0, 10),
      s = rep(1, 3),
      g_init = NULL,
      fix_g = FALSE,
      output = default.ebnm.output
    ),
    error = function(e) stop(error.msg))

  # Check for nonnegativity/nonpositivity:
  if (all(coef(tmp) >= 0)) {
    return(1)
  } else if (all(coef(tmp) <= 0)) {
    return(-1)
  } else {
    return(0)
  }
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
