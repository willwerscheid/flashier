#' Nullcheck flash factors
#'
#' Sets factors to zero if doing so improves the overall fit.
#'
#' @param flash A \code{flash} or \code{flash.fit} object.
#'
#' @param kset A vector of integers specifying which factors to nullcheck.
#'   If \code{kset = NULL}, then all existing factors will be checked.
#'
#' @param remove Whether to remove factors that have been set to zero from the
#'   flash object. Note that this might change the indices of existing factors.
#'
#' @param tol The tolerance parameter: if a factor does not improve the ELBO
#'   by at least \code{tol}, then it will be set to zero.
#'
#' @param verbose.lvl When and how to display progress updates. For nullchecks,
#'   updates are only displayed when \code{verbose.lvl} > 0.
#'
#' @export
#'
flash.nullcheck <- function(flash,
                            kset = NULL,
                            remove = TRUE,
                            tol = set.default.tol(flash),
                            verbose.lvl = get.verbose.lvl(flash)) {
  fit <- get.fit(flash)

  if (is.null(kset)) {
    if (get.n.factors(fit) > 0) {
      kset <- 1:get.n.factors(fit)
    } else {
      announce.no.nullchk(verbose.lvl)
      verbose.lvl <- 0
    }
  }
  must.be.valid.kset(fit, kset)

  must.be.numeric(tol, allow.infinite = TRUE, allow.null = FALSE)
  must.be.integer(verbose.lvl, lower = -1, upper = 3)

  announce.nullchk(verbose.lvl, n.factors = length(kset))

  for (k in kset) {
    fit <- nullcheck.factor(fit, k, verbose.lvl, tol)
  }

  if (length(kset) > 0 && !nullchk.failed(fit)) {
    report.nullchk.success(verbose.lvl)
  }

  if (sum(is.zero(fit)) > 0) {
    announce.wrapup(verbose.lvl)
    if (remove) {
      flash <- flash.remove.factors(fit, which(is.zero(fit)))
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
