#' Flashier prior families
#'
#' Families of distributions from which priors on loadings are to be estimated.
#'
#' \code{prior.normal} is estimated from the family of normal distributions
#'   \eqn{N(0, 1/a)}. \code{prior.point.normal} is estimated from the
#'   two-parameter family \eqn{\pi \delta(0) + (1-\pi) N(0, 1/a)}.
#'   \code{prior.point.laplace} replaces the normal slab with a Laplace slab,
#'   and \code{prior.nonzero.mode} adds a third parameter \eqn{\mu},
#'   so that the prior is estimated from the family \eqn{\pi \delta(\mu) +
#'   (1-\pi) N(\mu, 1/a)}.
#'
#' \code{prior.scale.normal.mix} is estimated from the family of scale mixtures
#'   of normals \eqn{\pi_0 \delta(0) + \pi_1 N(0, 1/a_1) + \ldots +
#'   \pi_m N(0, 1/a_m)}. \code{prior.unimodal.symmetric} replaces
#'   \eqn{N(0, 1/a_j)} with \eqn{Unif(-a_j, a_j)}; \code{prior.nonnegative}
#'   uses mixture components \eqn{Unif(0, a_j)}; and
#'   \code{prior.nonpositive} uses components \eqn{Unif(-a_j, 0)}. Finally,
#'   \code{prior.unimodal} is (approximately) estimated from the general family
#'   of unimodal distributions with mode zero by using both components
#'   \eqn{Unif(-a_j, 0)} and components \eqn{Unif(0, a_j)}.
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
#' @export
as.prior <- function(ebnm.fn = ebnm.nowarn, sign = 0, ...) {
  return(list(list(sign = sign,
                   ebnm.fn = function(x, s, g, fixg, output) {
                     ebnm.fn(x, s, g_init = g, fix_g = fixg, output = output, ...)
                   })))
}

#' @importFrom ebnm ebnm
#' @export
ebnm.nowarn <- function(...) {
  withCallingHandlers(res <- ebnm(...),
                      warning = function(w) {
                        if (startsWith(conditionMessage(w),
                                       "Optimization failed to converge. Results"))
                          invokeRestart("muffleWarning")
                      })
  return(res)
}

#' @rdname prior.families
#' @export
prior.normal <- function(...) {
  args <- as.prior.args(prior.family = "normal",
                        optmethod = "optimize", ...)
  return(do.call(as.prior, args))
}

#' @rdname prior.families
#' @export
prior.point.normal <- function(...) {
  args <- as.prior.args(prior.family = "point_normal",
                        optmethod = "nlm", ...)
  return(do.call(as.prior, args))
}

#' @rdname prior.families
#' @export
prior.point.laplace <- function(...) {
  args <- as.prior.args(prior.family = "point_laplace",
                        optmethod = "nlm", ...)
  return(do.call(as.prior, args))
}

#' @rdname prior.families
#' @export
prior.nonzero.mode <- function(...) {
  args <- as.prior.args(prior.family = "point_normal",
                        optmethod = "nlm",
                        mode = "estimate", ...)
  return(do.call(as.prior, args))
}

#' @rdname prior.families
#' @export
prior.normal.scale.mix <- function(...) {
  args <- as.prior.args(prior.family = "normal_scale_mixture",
                        optmethod = "mixsqp", ...)
  return(do.call(as.prior, args))
}

#' @rdname prior.families
#' @export
prior.unimodal <- function(...) {
  args <- as.prior.args(prior.family = "unimodal",
                        optmethod = "mixsqp", ...)
  return(do.call(as.prior, args))
}

#' @rdname prior.families
#' @export
prior.unimodal.symmetric <- function(...) {
  args <- as.prior.args(prior.family = "unimodal_symmetric",
                        optmethod = "mixsqp", ...)
  return(do.call(as.prior, args))
}

#' @rdname prior.families
#' @export
prior.nonnegative <- function(...) {
  args <- as.prior.args(prior.family = "unimodal_nonnegative",
                        optmethod = "mixsqp",
                        sign = 1, ...)
  return(do.call(as.prior, args))
}

#' @rdname prior.families
#' @export
prior.nonpositive <- function(...) {
  args <- as.prior.args(prior.family = "unimodal_nonpositive",
                        optmethod = "mixsqp",
                        sign = -1, ...)
  return(do.call(as.prior, args))
}

# Add default control parameters to priors.
as.prior.args <- function(prior.family, optmethod, ...) {
  args <- list(...)

  args$prior_family <- prior.family

  if (is.null(args$control)) {
    args$control <- list()
  }

  if (identical(optmethod, "nlm")) {
    args$control <- modifyList(nlm.defaults(), args$control)
  } else if (identical(optmethod, "mixsqp")) {
    args$control <- modifyList(mixsqp.defaults(), args$control)
  }

  return(args)
}

mixsqp.defaults <- function() {
  return(list(maxiter.sqp = 10))
}

nlm.defaults <- function() {
  return(list(iterlim = 10))
}
