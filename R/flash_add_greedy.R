#' Greedily add factors to a flash object
#'
#' Adds factors to a flash object in a "greedy" manner. Up to \code{Kmax}
#'   factors are added one at a time. At each step, \code{flash.add.greedy}
#'   attempts to find an optimal additional (rank-one) factor given all
#'   previously added factors. The additional factor is retained if it
#'   increases the ELBO; otherwise, fitting terminates.
#'
#' @inheritParams flash
#'
#' @param flash A \code{flash} or \code{flash.fit} object to which factors are
#'   to be added.
#'
#' @param Kmax The maximum number of factors to be added. This will not
#'   necessarily be the total number of factors added by
#'   \code{flash.add.greedy}, since factors are only added as long as they
#'   increase the variational lower bound on the log likelihood for the model.
#'
#' @param init.fn The function used to initialize factors. Functions
#'   \code{\link{init.fn.default}}, \code{\link{init.fn.softImpute}}, and
#'   \code{\link{init.fn.irlba}} have been supplied, but custom initialization
#'   functions may also be used. In particular, it is important to use an
#'   appropriate initialization function when loadings must be constrained in
#'   some fashion (otherwise, the greedy algorithm can stop adding factors
#'   prematurely). Custom initialization functions should accept a single parameter
#'   \code{flash} (the name of the parameter is unimportant) and should output
#'   a list consisting of two vectors, which will be used as initial values for
#'   the new loadings \eqn{\ell_k} and the new factor \eqn{f_k}. Typically,
#'   a custom initialization function will extract the matrix of residuals from
#'   \code{flash} using the method \code{resid(flash)} and then return a
#'   (possibly constrained) rank-one approximation to the matrix of residuals.
#'   See \strong{Examples} below.
#'
#' @param extrapolate Whether to use an extrapolation technique
#'   inspired by Ang and Gillis (2019) to accelerate the fitting process.
#'   Control parameters are handled via global options and can be set by
#'   calling \code{options("extrapolate.control") <- control.param}.
#'
#' @param conv.crit.fn The function used to determine whether convergence has
#'   occurred. TODO: details.
#'
#' @param tol The convergence tolerance.
#'
#' @param maxiter The maximum number of iterations when optimizing a greedily
#'   added factor.
#'
#' @examples
#' # Increase the maximum number of iterations in the default initialization
#' #   method.
#' fl <- flash.init(gtex) %>%
#'   flash.add.greedy(init.fn = function(f) init.fn.default(f, maxiter = 500))
#'
#' # Fit a semi-nonnegative matrix factorization.
#' snmf.fl <- flash.init(gtex) %>%
#'   flash.add.greedy(
#'     ebnm.fn = c(ebnm::ebnm_unimodal_nonnegative, ebnm::ebnm_point_normal),
#'     init.fn = function(f) init.fn.default(f, dim.signs = c(1, 0))
#'   )
#'
#' # Use a custom initialization function that wraps function nnmf from
#' #   package NNLM.
#' nnmf.init.fn <- function(flash) {
#'   nnmf.res <- NNLM::nnmf(resid(flash), verbose = FALSE)
#'   return(list(as.vector(nnmf.res$W), as.vector(nnmf.res$H)))
#' }
#' fl.nnmf <- flash.init(gtex) %>%
#'   flash.add.greedy(ebnm.fn = ebnm::ebnm_unimodal_nonnegative,
#'                    init.fn = nnmf.init.fn)
#'
#' @seealso \code{\link{init.fn.default}}, \code{\link{init.fn.softImpute}},
#'   \code{\link{init.fn.irlba}}
#'
#' @importFrom ebnm ebnm_point_normal
#'
#' @export
#'
flash.add.greedy <- function(flash,
                             Kmax = 1,
                             ebnm.fn = ebnm::ebnm_point_normal,
                             init.fn = init.fn.default,
                             extrapolate = FALSE,
                             conv.crit.fn = calc.obj.diff,
                             tol = set.default.tol(flash),
                             maxiter = 500,
                             verbose.lvl = get.verbose.lvl(flash)) {
  flash <- get.fit(flash)

  must.be.integer(Kmax, lower = 1, allow.null = FALSE)
  must.be.integer(maxiter, lower = 1, allow.null = FALSE)
  must.be.integer(verbose.lvl, lower = -1, upper = 3, allow.null = FALSE)

  ebnm.fn <- handle.ebnm.fn(ebnm.fn, get.dim(flash))

  if (missing(tol)) {
    report.tol.setting(verbose.lvl, tol)
  } else {
    must.be.numeric(tol, allow.infinite = FALSE, allow.null = FALSE)
  }

  if (uses.R(flash)) {
    if (!missing(extrapolate) && extrapolate) {
      stop("Greedy extrapolation has not yet been implemented for the 'zero' ",
           "variance structure.")
    } else {
      extrapolate <- FALSE
    }
  }

  extrapolate.control <- getOption("extrapolate.control", list())
  extrapolate.param <- set.extrapolate.param(extrapolate.control)

  verbose.fns <- get.verbose.fns(flash)
  verbose.colnames <- get.verbose.colnames(flash)
  verbose.colwidths <- get.verbose.colwidths(flash)

  if (verbose.lvl == -1) {
    print.tab.delim.table.header(verbose.colnames)
  }

  factors.added <- 0
  greedy.failed <- FALSE

  while (factors.added < Kmax && !greedy.failed) {
    announce.add.factor(verbose.lvl, k = get.next.k(flash))

    factor <- init.factor(flash, init.fn)
    factor <- set.ebnm.fn(factor, ebnm.fn)

    announce.factor.opt(verbose.lvl)
    print.table.header(verbose.lvl,
                       verbose.colnames,
                       verbose.colwidths,
                       backfit = FALSE)

    iter <- 0
    conv.crit <- Inf
    if (extrapolate) {
      old.f <- factor
      extrapolate.param <- init.beta(extrapolate.param)
    }
    while (conv.crit > tol && iter < maxiter) {
      iter <- iter + 1

      if (extrapolate) {
        proposed.factor <- extrapolate.f(factor, old.f, extrapolate.param)
        proposed.factor <- update.factor(proposed.factor, flash)
      }

      old.f <- factor
      if (!extrapolate) {
        factor <- update.factor(factor, flash)
      } else if (get.obj(proposed.factor) - get.obj(factor) < tol) {
        factor <- update.factor(factor, flash)
        extrapolate.param <- decelerate(extrapolate.param)
      } else {
        factor <- proposed.factor
        extrapolate.param <- accelerate(extrapolate.param)
      }

      obj.diff <- get.obj(factor) - get.obj(old.f)
      if (is.obj.valid(old.f) && obj.diff < 0) {
        report.greedy.obj.decrease(verbose.lvl, obj.diff)
        factor <- old.f
        break
      }

      info <- calc.update.info(factor, old.f, conv.crit.fn, verbose.fns)
      conv.crit <- get.conv.crit(info)
      print.table.entry(verbose.lvl,
                        verbose.colwidths,
                        iter,
                        info,
                        get.next.k(flash),
                        backfit = FALSE)
    }

    if (iter == maxiter) {
      report.maxiter.reached(verbose.lvl)
    }

    if (!is.zero(factor) &&
        (get.obj(factor) > get.obj(flash) + tol
         || !is.obj.valid(flash, factor))) {
      flash <- add.new.factor.to.flash(factor, flash)
      factors.added <- factors.added + 1
    } else {
      greedy.failed <- TRUE
    }

    report.add.factor.result(verbose.lvl, greedy.failed, get.obj(flash))
  }

  announce.wrapup(verbose.lvl)
  flash <- wrapup.flash(flash, output.lvl = 3L)

  report.completion(verbose.lvl)

  return(flash)
}
