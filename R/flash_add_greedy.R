#' Greedily add factors to a flash object
#'
#' Adds factor/loadings pairs to a flash object in a "greedy" manner. Up to
#'   \code{Kmax} pairs are added one at a time. At each step, \code{flash.add.greedy}
#'   attempts to find an optimal additional (rank-one) factor given all
#'   previously added factors. The additional factor is retained if it
#'   increases the variational lower bound (ELBO); otherwise, fitting terminates.
#'
#' @inheritParams flash
#'
#' @param flash A \code{flash} or \code{flash.fit} object to which factors are
#'   to be added.
#'
#' @param Kmax The maximum number of factors to be added. This will not
#'   necessarily be the total number of factors added by
#'   \code{flash.add.greedy}, since factors are only added as long as they
#'   increase the ELBO.
#'
#' @param init.fn The function used to initialize factor/loadings pairs. Functions
#'   \code{\link{flash_init_greedy_default}}, \code{\link{flash_init_greedy_softImpute}}, and
#'   \code{\link{flash_init_greedy_irlba}} have been supplied
#'   (\code{\link{flash_init_greedy_softImpute}} can yield better results than the
#'   default initialization function when there is missing data). Custom
#'   initialization functions may also be used. If \code{init.fn = NULL} then
#'   \code{flash_init_greedy_default} will be used; \code{flash.add.greedy} will
#'   attempt to set argument \code{dim.signs} appropriately via test calls to
#'   the EBNM function(s) specified by \code{ebnm.fn}. If factors or loadings
#'   are constrained in some other fashion (e.g., bounded support), then the
#'   initialization function should be modified to account for the constraints
#'   --- otherwise, the greedy algorithm can stop adding factor/loadings pairs
#'   too early. Custom initialization functions should accept a single
#'   parameter referring to a \code{\link{flash.fit}} object and should output
#'   a list consisting of two vectors, which will be used as initial values for
#'   the new loadings \eqn{\ell_k} and the new factor \eqn{f_k}. Typically,
#'   a custom initialization function will extract the matrix of residuals from
#'   the \code{flash.fit} object using the method \code{residuals.flash.fit} and
#'   then return a (possibly constrained) rank-one approximation to the matrix
#'   of residuals. See \strong{Examples} below.
#'
#' @param warmstart Whether to use "warmstarts" when solving the EBNM
#'   subproblems by initializing solutions at the previous value of the fitted
#'   prior \eqn{\hat{g}}. An important side effect of warmstarts for
#'   \code{ashr}-like prior families is to fix the grid at its initial setting.
#'   Fixing the grid can lead to poor fits if there
#'   are large changes in the scale of the estimated prior over the
#'   course of the fitting process. However, allowing the grid to
#'   vary can occasionally result in decreases in ELBO.
#'
#' @param extrapolate Whether to use an extrapolation technique
#'   inspired by Ang and Gillis (2019) to accelerate the fitting process.
#'   Control parameters are handled via global options and can be set by
#'   calling \code{options("extrapolate.control") <- control.param}.
#'
#' @param maxiter The maximum number of iterations when optimizing a greedily
#'   added factor.
#'
#' @examples
#' # Increase the maximum number of iterations in the default initialization
#' #   method.
#' fl <- flash.init(gtex) %>%
#'   flash.add.greedy(init.fn = function(f) flash_init_greedy_default(f, maxiter = 500))
#'
#' # Use a custom initialization function that wraps function nmf from
#' #   package RcppML.
#' nmf.init.fn <- function(f) {
#'   nmf.res <- RcppML::nmf(resid(f), k = 1, verbose = FALSE)
#'   return(list(as.vector(nmf.res$w), as.vector(nmf.res$h)))
#' }
#' fl.nmf <- flash.init(gtex) %>%
#'   flash.add.greedy(ebnm.fn = ebnm_unimodal_nonnegative,
#'                    init.fn = nmf.init.fn)
#'
#' @seealso \code{\link{flash_init_greedy_default}}, \code{\link{flash_init_greedy_softImpute}},
#'   \code{\link{flash_init_greedy_irlba}}
#'
#' @return A \code{\link{flash}} object.
#'
#' @importFrom ebnm ebnm_point_normal
#'
#' @export
#'
flash.add.greedy <- function(flash,
                             Kmax = 1,
                             ebnm.fn = ebnm_point_normal,
                             init.fn = NULL,
                             extrapolate = FALSE,
                             warmstart = FALSE,
                             maxiter = 500,
                             tol = NULL,
                             verbose = NULL) {
  flash <- get.fit(flash)

  tol <- handle.tol.param(tol, flash)
  verbose.lvl <- handle.verbose.param(verbose, flash)

  must.be.integer(Kmax, lower = 1, allow.null = FALSE)
  must.be.integer(maxiter, lower = 1, allow.null = FALSE)
  must.be.integer(verbose.lvl, lower = -1, upper = 3, allow.null = FALSE)

  ebnm.chk <- handle.ebnm.fn(ebnm.fn, get.dim(flash))
  ebnm.fn <- ebnm.chk$ebnm.fn
  if (is.null(init.fn)) {
    init.fn <- function(f) flash_init_greedy_default(f, dim.signs = ebnm.chk$dim.signs)
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

  flash <- set.gwarmstart(flash, warmstart)

  verbose.fns <- get.verbose.fns(flash)
  verbose.colnames <- get.verbose.colnames(flash)
  verbose.colwidths <- get.verbose.colwidths(flash)

  conv.crit.fn <- get.conv.crit.fn(flash)

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
