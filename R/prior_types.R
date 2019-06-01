#' Flashier prior types
#'
#' Classes of distributions from which priors on loadings are to be estimated.
#'
#' A \code{normal} prior is estimated from the class of normal distributions
#'   \eqn{N(0, 1/a)}. A \code{point.normal} prior is estimated from the
#'   two-parameter class \eqn{\pi \delta(0) + (1-\pi) N(0, 1/a)}. A
#'   \code{point.laplace} prior replaces the normal slab with a Laplace slab,
#'   and a \code{nonzero.mode} prior adds a third parameter \eqn{\mu},
#'   so that the prior is estimated from the class \eqn{\pi \delta(\mu) +
#'   (1-\pi) N(\mu, 1/a)}.
#'
#' A \code{normal.mix} prior is estimated from the class of mixtures of
#'   normals \eqn{\pi_0 \delta(0) + \pi_1 N(0, 1/a_1) + \ldots +
#'   \pi_m N(0, 1/a_m)}. The \code{uniform.mix} class of priors replaces
#'   \eqn{N(0, 1/a_j)} with \eqn{Unif(-a_j, a_j)}; the \code{nonnegative}
#'   class uses mixture components \eqn{Unif(0, a_j)}; and the
#'   \code{nonpositive} class uses components \eqn{Unif(-a_j, 0)}.
#'
#' Custom prior types can also be created. These should be lists containing
#'   exactly one element, which itself must be a list containing fields
#'   \code{ebnm.fn}, \code{ebnm.param}, and \code{sign}. That is, they should
#'   be of form \code{custom.prior = list(list(ebnm.fn = my.ebnm.fn,
#'   ebnm.param = my.ebnm.param, sign = my.sign))}. \code{ebnm.fn} gives
#'   the function used to solve the empirical Bayes normal means problem.
#'   Typically, this will be \code{ebnm.pn}, which is a wrapper to
#'   \code{ebnm::ebnm}, or \code{ebnm.ash}, a wrapper to \code{ashr::ash},
#'   but custom functions may also be used. For details, see
#'   \code{\link{ebnm.pn}}.
#'   \code{ebnm.param} lists any additional parameters to be passed to \code{ebnm.fn} (this
#'   field is useful when multiple prior types share the same \code{ebnm.fn}).
#'   \code{sign} should be set to +1 for classes of distributions with
#'   nonnegative support, -1 for classes with nonpositive support, and 0 for
#'   all other classes. It is only used when initializing new factors.
#'
#' @param ... Additional parameters to be passed to the function used to solve
#'   the empirical Bayes normal means problem.
#'
#' @rdname prior.type
#'
#' @export
#'
normal <- function(...) {
  param <- list(prior_type = "normal")
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.pn,
                   ebnm.param = modifyList(param, list(...)))))
}

#' @rdname prior.type
#'
#' @export
#'
point.normal <- function(...) {
  param <- list(prior_type = "point_normal")
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.pn,
                   ebnm.param = modifyList(param, list(...)))))
}

#' @rdname prior.type
#'
#' @export
#'
point.laplace <- function(...) {
  param <- list(prior_type = "point_laplace")
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.pn,
                   ebnm.param = modifyList(param, list(...)))))
}

#' @rdname prior.type
#'
#' @export
#'
nonzero.mode <- function(...) {
  param <- list(prior_type = "point_normal", fix_mu = FALSE)
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.pn,
                   ebnm.param = modifyList(param, list(...)))))
}

#' @rdname prior.type
#'
#' @export
#'
normal.mix <- function(...) {
  param <- list(mixcompdist = "normal")
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.ash,
                   ebnm.param = modifyList(param, list(...)))))
}

#' @rdname prior.type
#'
#' @export
#'
uniform.mix <- function(...) {
  param <- list(mixcompdist = "uniform")
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.ash,
                   ebnm.param = modifyList(param, list(...)))))
}

#' @rdname prior.type
#'
#' @export
#'
nonnegative <- function(...) {
  param <- list(mixcompdist = "+uniform")
  return(list(list(sign = 1,
                   ebnm.fn = ebnm.ash,
                   ebnm.param = modifyList(param, list(...)))))
}

#' @rdname prior.type
#'
#' @export
#'
nonpositive <- function(...) {
  param <- list(mixcompdist = "-uniform")
  return(list(list(sign = -1,
                   ebnm.fn = ebnm.ash,
                   ebnm.param = modifyList(param, list(...)))))
}
