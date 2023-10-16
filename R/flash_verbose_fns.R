#' Display the current ELBO
#'
#' Displays the value of the variational lower bound (ELBO) at the current
#'   iteration.
#'
#' @details This function is an example of a function that may be passed to
#'   parameter \code{fns} in function \code{\link{flash_set_verbose}} to
#'   customize the output that is printed after each greedy or backfitting
#'   iteration. See \code{\link{flash_set_verbose}} for details and examples.
#'
#' @param curr The \code{\link{flash_fit}} object from the current iteration.
#'
#' @param prev The \code{\link{flash_fit}} object from the previous iteration.
#'
#' @param k Only used during sequential backfits (that is, calls to
#'   \code{\link{flash_backfit}} where \code{extrapolate = FALSE}). It then
#'   takes the index of the factor/loadings pair currently being optimized.
#'
#' @return A character string, suitable for printing progress updates.
#'
#' @seealso \code{\link{flash_verbose_elbo_diff}},
#'   \code{\link{flash_verbose_max_chg}},
#'   \code{\link{flash_verbose_max_chg_L}},
#'   \code{\link{flash_verbose_max_chg_F}}
#'
#' @export
#'
flash_verbose_elbo <- function(curr, prev, k) {
  obj <- get.curr.obj(curr, prev, k)
  if (is.na(obj))
    return("NA")
  return(formatC(obj, format = "f", digits = 2))
}

#' Display the difference in ELBO
#'
#' Displays the difference in the variational lower bound (ELBO) from one
#'   iteration to the next.
#'
#' @inherit flash_verbose_elbo details
#'
#' @inheritParams flash_verbose_elbo
#'
#' @inherit flash_verbose_elbo return
#'
#' @seealso \code{\link{flash_verbose_elbo}}, \code{\link{flash_verbose_max_chg}},
#'   \code{\link{flash_verbose_max_chg_L}}, \code{\link{flash_verbose_max_chg_F}}
#'
#' @export
#'
flash_verbose_elbo_diff <- function(curr, prev, k) {
  obj.diff <- calc.obj.diff(curr, prev, k)
  if (is.infinite(obj.diff))
    obj.diff <- NA
  return(obj.diff)
}

#' Display the maximum difference in scaled loadings and factors
#'
#' Displays the maximum (absolute) change over all (posterior expected values for)
#'   loadings \eqn{\ell_{ik}} and factors \eqn{f_{jk}}. At each iteration, the
#'   loadings vectors \eqn{\ell_{\cdot 1}, \ldots, \ell_{\cdot K}} and factors
#'   \eqn{f_{\cdot 1}, \ldots, f_{\cdot K}} are \eqn{L^2}-normalized.
#'
#' @inherit flash_verbose_elbo details
#'
#' @inheritParams flash_verbose_elbo
#'
#' @inherit flash_verbose_elbo return
#'
#' @seealso \code{\link{flash_verbose_elbo}}, \code{\link{flash_verbose_elbo_diff}},
#'   \code{\link{flash_verbose_max_chg_L}}, \code{\link{flash_verbose_max_chg_F}}
#'
#' @export
#'
flash_verbose_max_chg <- function(curr, prev, k) {
  return(calc.max.abs.chg.EF(curr, prev, k, n = NULL))
}

#' Display the maximum difference in scaled loadings
#'
#' Displays the maximum (absolute) change over all (posterior expected values for)
#'   loadings \eqn{\ell_{ik}}. At each iteration, the loadings vectors
#'   \eqn{\ell_{\cdot 1}, \ldots, \ell_{\cdot K}} are \eqn{L^2}-normalized.
#'
#' @inherit flash_verbose_elbo details
#'
#' @inheritParams flash_verbose_elbo
#'
#' @inherit flash_verbose_elbo return
#'
#' @seealso \code{\link{flash_verbose_elbo}}, \code{\link{flash_verbose_elbo_diff}},
#'   \code{\link{flash_verbose_max_chg}}, \code{\link{flash_verbose_max_chg_F}}
#'
#' @export
#'
flash_verbose_max_chg_L <- function(curr, prev, k) {
  return(calc.max.abs.chg.EF(curr, prev, k, n = 1))
}

#' Display the maximum difference in scaled factors
#'
#' Displays the maximum (absolute) change over all (posterior expected values for)
#'   factors \eqn{f_{jk}}. At each iteration, the factors
#'   \eqn{f_{\cdot 1}, \ldots, f_{\cdot K}} are \eqn{L^2}-normalized.
#'
#' @inherit flash_verbose_elbo details
#'
#' @inheritParams flash_verbose_elbo
#'
#' @inherit flash_verbose_elbo return
#'
#' @seealso \code{\link{flash_verbose_elbo}}, \code{\link{flash_verbose_elbo_diff}},
#'   \code{\link{flash_verbose_max_chg}}, \code{\link{flash_verbose_max_chg_L}}
#'
#' @export
#'
flash_verbose_max_chg_F <- function(curr, prev, k) {
  return(calc.max.abs.chg.EF(curr, prev, k, n = 2))
}

calc.update.info <- function(curr, prev, conv.crit.fn, verbose.fns, k = NULL) {
  if (length(verbose.fns) == 0) {
    all.fns <- list(conv.crit.fn)
  } else {
    all.fns <- verbose.fns
    conv.crit.idx <- which(sapply(verbose.fns, identical, conv.crit.fn))
    if (length(conv.crit.idx) == 0)
      all.fns <- c(all.fns, conv.crit.fn)
  }

  update.info <- lapply(all.fns, do.call, list(curr, prev, k))

  # Put a copy of the convergence criterion at the end of the list so that
  #   it's easy to find.
  if ((length(verbose.fns) > 0) && (length(conv.crit.idx) > 0))
    update.info <- c(update.info, update.info[conv.crit.idx])

  return(update.info)
}

get.verbose.info <- function(update.info) {
  return(update.info[-length(update.info)])
}

get.curr.obj <- function(curr, prev, k) {
  if (!is.obj.valid(curr))
    return(NA)
  return(get.obj(curr))
}

calc.obj.diff <- function(curr, prev, k) {
  if (!is.obj.valid(prev) || !is.obj.valid(curr))
    return(Inf)
  return(get.obj(curr) - get.obj(prev))
}

calc.max.abs.chg.EF <- function(curr, prev, k, n = NULL) {
  if (!is.null(k))
    return(calc.max.abs.chg(get.EF.k(curr, k, n), get.EF.k(prev, k, n)))
  return(calc.max.abs.chg(get.EF(curr, n), get.EF(prev, n)))
}

calc.max.chg.EF <- function(curr, prev, k, n = NULL) {
  if (!is.null(k))
    return(calc.max.chg(get.EF.k(curr, k, n), get.EF.k(prev, k, n)))
  return(calc.max.chg(get.EF(curr, n), get.EF(prev, n)))
}

which.max.chg.EF <- function(curr, prev, k, n = NULL) {
  if (!is.null(k))
    return(which.max.chg(get.EF.k(curr, k, n), get.EF.k(prev, k, n)))
  return(which.max.chg(get.EF(curr, n), get.EF(prev, n)))
}

get.sparsity <- function(curr, prev, k, n) {
  if (!is.null(k)) {
    g <- get.g.k(curr, k, n)
  } else {
    g <- get.g(curr, n)
  }
  return(g$pi[1])
}

get.exclusion.count <- function(curr, prev, k, n) {
  if (!is.null(k)) {
    return(length(get.exclusions(curr)[[k]][[n]]))
  } else {
    return(length(get.exclusions(curr)[[n]]))
  }
}

calc.max.abs.chg <- function(curr, prev) {
  curr <- l2.normalize.and.stack(curr)
  prev <- l2.normalize.and.stack(prev)
  return(max(abs(curr - prev)))
}

calc.max.chg <- function(curr, prev) {
  curr <- l2.normalize.and.stack(curr)
  prev <- l2.normalize.and.stack(prev)
  max.increase <- max(curr - prev)
  max.decrease <- max(prev - curr)
  if (max.increase > max.decrease)
    return(max.increase)
  return(-max.decrease)
}

which.max.chg <- function(curr, prev) {
  curr <- l2.normalize.and.stack(curr)
  prev <- l2.normalize.and.stack(prev)
  return(which.max(apply(abs(curr - prev), 1, max)))
}

l2.normalize.and.stack <- function(x) {
  if (is.list(x)) {
    norm.x <- lapply(x, l2.normalize)
    return(do.call(rbind, norm.x))
  }
  return(l2.normalize(x))
}

l2.normalize <- function(x) {
  if (is.matrix(x)) {
    norm <- sqrt(colSums(x^2))
  } else {
    norm <- sqrt(sum(x^2))
  }
  norm[norm == 0] <- 1
  if (is.matrix(x))
    return(x / matrix(norm, nrow = nrow(x), ncol = ncol(x), byrow = TRUE))
  return(matrix(x / norm, ncol = 1))
}
