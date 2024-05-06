#' Greedily add factors to a flash object
#'
#' Adds factor/loadings pairs to a flash object in a "greedy" manner. Up to
#'   \code{Kmax} pairs are added one at a time. At each step, \code{flash_greedy}
#'   attempts to find an optimal additional (rank-one) factor given all
#'   previously added factors. The additional factor is retained if it
#'   increases the variational lower bound (ELBO); otherwise, fitting terminates.
#'   See \code{\link{flash}} for examples of usage.
#'
#' @inheritParams flash
#'
#' @param flash A \code{flash} or \code{flash_fit} object to which factors are
#'   to be added.
#'
#' @param Kmax The maximum number of factors to be added. This will not
#'   necessarily be the total number of factors added by
#'   \code{flash_greedy}, since factors are only added as long as they
#'   increase the ELBO.
#'
#' @param init_fn The function used to initialize factor/loadings pairs. Functions
#'   \code{\link{flash_greedy_init_default}}, \code{\link{flash_greedy_init_softImpute}}, and
#'   \code{\link{flash_greedy_init_irlba}} have been supplied; note, in particular, that
#'   \code{\link{flash_greedy_init_softImpute}} can yield better results than the
#'   default initialization function when there is missing data. Custom
#'   initialization functions may also be used. If \code{init_fn = NULL} then
#'   \code{\link{flash_greedy_init_default}} will be used, with an attempt made to set
#'   argument \code{sign_constraints} appropriately via test calls to
#'   the EBNM function(s) specified by parameter \code{ebnm_fn}. If factors or loadings
#'   are constrained in some other fashion (e.g., bounded support), then the
#'   initialization function should be modified to account for the constraints
#'   --- otherwise, the greedy algorithm can stop adding factor/loadings pairs
#'   too early. Custom initialization functions should accept a single
#'   parameter referring to a \code{\link{flash_fit}} object and should output
#'   a list consisting of two vectors, which will be used as initial values for
#'   the new loadings \eqn{\ell_{\cdot k}} and the new factor \eqn{f_{\cdot k}}. Typically,
#'   a custom initialization function will extract the matrix of residuals from
#'   the \code{flash_fit} object using method \code{residuals.flash_fit} and
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
#'   added factor/loadings pair.
#'
#' @param tol The convergence tolerance parameter. At each iteration, the fit
#'   is compared to the fit from the previous iteration using a convergence
#'   criterion function (by default, the difference in ELBO, but the criterion
#'   can be changed via \code{\link{flash_set_conv_crit}}). When
#'   the value returned by this function is less than or equal to \code{tol},
#'   the newly added factor/loadings pair is considered to have "converged,"
#'   so that \code{flash_greedy} moves on and attempts to add another new
#'   factor (or, if the maximum number of factors \code{Kmax} has been reached,
#'   the process terminates). Note that
#'   specifying \code{tol} here will override any value set by
#'   \code{flash_set_conv_crit}; to use the "global" tolerance parameter,
#'   \code{tol} must be left unspecified (\code{NULL}).
#'   If \code{tol = NULL} and a global tolerance
#'   parameter has not been set, then the default
#'   tolerance used is \eqn{np\sqrt{\epsilon}}, where \eqn{n} is the
#'   number of rows in the dataset, \eqn{p} is the number of columns, and
#'   \eqn{\epsilon} is equal to \code{\link{.Machine}$double.eps}.
#'
#' @examples
#' # The following are examples of advanced usage. See ?flash for basic usage.
#'
#' # Increase the maximum number of iterations in the default initialization
#' #   method.
#' my_init_fn <- function(f) flash_greedy_init_default(f, maxiter = 500)
#' fl <- flash_init(gtex) |>
#'   flash_greedy(init_fn = my_init_fn)
#'
#' # Use a custom initialization function that wraps function nmf from
#' #   package RcppML.
#' nmf_init_fn <- function(f) {
#'   nmf_res <- RcppML::nmf(resid(f), k = 1, verbose = FALSE)
#'   return(list(as.vector(nmf_res$w), as.vector(nmf_res$h)))
#' }
#' fl.nmf <- flash_init(gtex) |>
#'   flash_greedy(ebnm_fn = ebnm_unimodal_nonnegative,
#'                init_fn = nmf_init_fn)
#'
#' @seealso \code{\link{flash_greedy_init_default}},
#'   \code{\link{flash_greedy_init_softImpute}},
#'   \code{\link{flash_greedy_init_irlba}}
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with up
#'   to \code{Kmax} new factor/loadings pairs "greedily" added.
#'
#' @importFrom ebnm ebnm_point_normal
#'
#' @export
#'
flash_greedy <- function(flash,
                         Kmax = 1,
                         ebnm_fn = ebnm_point_normal,
                         init_fn = NULL,
                         extrapolate = FALSE,
                         warmstart = FALSE,
                         maxiter = 500,
                         tol = NULL,
                         verbose = NULL) {
  flash <- get.fit(flash)

  tol <- handle.tol.param(tol, flash)
  verbose.lvl <- handle.verbose.param(verbose, flash)

  if (is.timed.out(flash)) {
    report.timeout.no.greedy(verbose.lvl)
    verbose.lvl <- 0
  }

  must.be.integer(Kmax, lower = 1, allow.null = FALSE)
  must.be.integer(maxiter, lower = 1, allow.null = FALSE)
  must.be.integer(verbose.lvl, lower = -1, upper = 3, allow.null = FALSE)

  ebnm.chk <- handle.ebnm.fn(ebnm_fn, get.dim(flash))
  ebnm.fn <- ebnm.chk$ebnm.fn
  if (is.null(init_fn)) {
    init_fn <- function(f) flash_greedy_init_default(
      f, sign_constraints = ebnm.chk$dim.signs
    )
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

  while (factors.added < Kmax && !greedy.failed && !is.timed.out(flash)) {
    announce.add.factor(verbose.lvl, k = get.next.k(flash))

    factor <- init.factor(flash, init_fn)
    factor <- set.ebnm.fn(factor, ebnm.fn)

    announce.factor.opt(verbose.lvl)
    print_table.header(verbose.lvl,
                       verbose.colnames,
                       verbose.colwidths,
                       backfit = FALSE)

    iter <- 0
    conv.crit <- Inf
    if (extrapolate) {
      old.f <- factor
      extrapolate.param <- init.beta(extrapolate.param)
    }
    while (conv.crit > tol && iter < maxiter && !is.timed.out(flash)) {
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
      print_table.entry(verbose.lvl,
                        verbose.colwidths,
                        iter,
                        info,
                        get.next.k(flash),
                        backfit = FALSE)
    }

    if (is.timed.out(flash)) {
      t.diff <- Sys.time() - get.timeout.set.time(flash)
      report.timeout.reached(verbose.lvl, t.diff)
      flash <- set.timeout.reached.flag(flash)
    }

    if (iter == maxiter) {
      report.maxiter.reached(verbose.lvl)
    }

    if (!is.zero(factor) &&
        (get.obj(factor) > get.obj(flash) + tol # TODO: remove + tol here
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
