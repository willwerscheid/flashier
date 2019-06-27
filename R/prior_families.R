#' Flashier prior families
#'
#' Families of distributions from which priors on loadings are to be estimated.
#'
#' \code{prior.normal} is estimated from the class of normal distributions
#'   \eqn{N(0, 1/a)}. \code{prior.point.normal} is estimated from the
#'   two-parameter class \eqn{\pi \delta(0) + (1-\pi) N(0, 1/a)}.
#'   \code{prior.point.laplace} replaces the normal slab with a Laplace slab,
#'   and \code{prior.nonzero.mode} adds a third parameter \eqn{\mu},
#'   so that the prior is estimated from the class \eqn{\pi \delta(\mu) +
#'   (1-\pi) N(\mu, 1/a)}.
#'
#' \code{prior.normal.mix} is estimated from the class of scale mixtures of
#'   normals \eqn{\pi_0 \delta(0) + \pi_1 N(0, 1/a_1) + \ldots +
#'   \pi_m N(0, 1/a_m)}. \code{prior.unimodal} replaces
#'   \eqn{N(0, 1/a_j)} with \eqn{Unif(-a_j, a_j)}; \code{prior.nonnegative}
#'   uses mixture components \eqn{Unif(0, a_j)}; and
#'   \code{prior.nonpositive} uses components \eqn{Unif(-a_j, 0)}.
#'
#' Custom prior classes can be created using the function \code{as.prior}.
#'
#' @param ebnm.fn The function used to solve the empirical Bayes normal means
#'   problem. Typically, this will be \code{ebnm::ebnm}, but custom functions
#'   may also be used as long as they have the same signature as
#'   \code{ebnm::ebnm}.
#'
#' @param sign Should be set to +1 for classes of distributions with
#'   nonnegative support and -1 for classes with nonpositive support. Only used
#'   when initializing new factors.
#'
#' @param ... Additional parameters to be passed to \code{ebnm::ebnm}.

#' @rdname prior.families
#' @importFrom ebnm ebnm
#' @export
as.prior <- function(ebnm.fn = ebnm, sign = 0, ...) {
  return(list(list(sign = sign,
                   ebnm.fn = function(x, s, g, fixg, output) {
                     ebnm.fn(x, s, g_init = g, fix_g = fixg, output = output, ...)
                   })))
}

#' @rdname prior.families
#' @export
prior.normal <- function(...) {
  return(as.prior(prior_family = "normal", ...))
}

#' @rdname prior.families
#' @export
prior.point.normal <- function(...) {
  return(as.prior(prior_family = "point_normal", ...))
}

#' @rdname prior.families
#' @export
prior.point.laplace <- function(...) {
  return(as.prior(prior_family = "point_laplace", ...))
}

#' @rdname prior.families
#' @export
prior.nonzero.mode <- function(...) {
  return(as.prior(prior_family = "point_normal", mode = "estimate", ...))
}

#' @rdname prior.families
#' @export
prior.normal.scale.mix <- function(...) {
  return(as.prior(prior_family = "normal_scale_mixture", ...))
}

#' @rdname prior.families
#' @export
prior.unimodal <- function(...) {
  return(as.prior(prior_family = "unimodal", ...))
}

#' @rdname prior.families
#' @export
prior.nonnegative <- function(...) {
  return(as.prior(prior_family = "unimodal_nonnegative", sign = 1, ...))
}

#' @rdname prior.families
#' @export
prior.nonpositive <- function(...) {
  return(as.prior(prior_family = "unimodal_nonpositive", sign = -1, ...))
}
