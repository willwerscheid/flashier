#' Nullcheck flash factors
#'
#' Sets factor/loadings pairs to zero if doing so improves the variational
#'   lower bound (ELBO). See \code{\link{flash}} for examples of usage.
#'
#' @param flash A \code{flash} or \code{flash_fit} object.
#'
#' @param kset A vector of integers specifying which factors to nullcheck.
#'   If \code{kset = NULL}, then all existing factors will be checked.
#'
#' @param remove Whether to remove factors that have been set to zero from the
#'   \code{flash} object. Note that this might change the indices of existing
#'   factors.
#'
#' @param tol The "tolerance" parameter: if a factor does not improve the ELBO
#'   by at least \code{tol}, then it will be set to zero. Note that
#'   \code{flash_nullcheck} does not respect "global" tolerance parameters set
#'   by \code{\link{flash_set_conv_crit}} (which only affects the convergence
#'   tolerance for greedy fits and backfits). The default tolerance is
#'   \eqn{np\sqrt{\epsilon}}, where \eqn{n} is the
#'   number of rows in the dataset, \eqn{p} is the number of columns, and
#'   \eqn{\epsilon} is equal to \code{\link{.Machine}$double.eps}.
#'
#' @param verbose When and how to display progress updates. For nullchecks,
#'   updates are only displayed when \code{verbose} > 0.
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with
#'   factors that do not improve the ELBO by at least \code{tol} either set
#'   to zero or removed (depending on the argument to parameter \code{remove}).
#'
#' @seealso \code{\link{flash_factors_remove}},
#'   \code{\link{flash_factors_set_to_zero}}
#'
#' @export
#'
flash_nullcheck <- function(flash,
                            kset = NULL,
                            remove = TRUE,
                            tol = NULL,
                            verbose = NULL) {
  fit <- get.fit(flash)

  must.be.numeric(tol, allow.infinite = TRUE, allow.null = TRUE)
  if (is.null(tol)) {
    tol <- sqrt(.Machine$double.eps) * prod(get.dims(fit))
  }
  verbose.lvl <- handle.verbose.param(verbose, fit)

  if (is.null(kset)) {
    if (get.n.factors(fit) > 0) {
      kset <- 1:get.n.factors(fit)
    } else {
      announce.no.nullchk(verbose.lvl)
      verbose.lvl <- 0
    }
  }
  must.be.valid.kset(fit, kset)

  must.be.integer(verbose.lvl, lower = -1, upper = 3)

  announce.nullchk(verbose.lvl, length(kset), sum(is.zero(fit)))

  for (k in kset) {
    fit <- nullcheck.factor(fit, k, verbose.lvl, tol)
  }

  if (length(kset) > 0 && !nullchk.failed(fit)) {
    report.nullchk.success(verbose.lvl)
  }

  if (sum(is.zero(fit)) > 0) {
    announce.wrapup(verbose.lvl)
    if (remove) {
      flash <- flash_factors_remove(fit, which(is.zero(fit)))
      report.factor.removal(verbose.lvl, sum(is.zero(fit)))
    } else {
      flash <- wrapup.flash(fit, output.lvl = 3L)
    }
  }

  report.completion(verbose.lvl)
  return(flash)
}

nullcheck.factor <- function(flash, k, verbose.lvl, tol) {
  if (!is.valid(flash, k) || is.zero(flash, k))
    return(flash)

  factor <- zero.factor(flash, k)
  factor <- update.R2.tau.and.obj(factor, flash)
  obj.diff <- get.obj(factor) - get.obj(flash)

  if (obj.diff >= -tol) {
    flash <- insert.factor(flash, factor)
    flash <- set.nullchk.fail.flag(flash)
    report.nullchk.failure(verbose.lvl, obj.diff, k)
  }

  return(flash)
}

zero.factor <- function(flash, k) {
  factor            <- list()
  factor$k          <- k
  factor$EF         <- r1.zeros(flash)
  factor$EF2        <- factor$EF
  factor$g          <- list()
  factor$KL         <- rep(0, get.dim(flash))
  factor$exclusions <- list()
  factor$is.zero    <- TRUE

  return(factor)
}
