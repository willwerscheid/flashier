#' Nullcheck flash factors
#'
#' Checks to see whether setting factors to zero improves the overall fit. If
#'   so, the factors are deleted.
#'
#' @inheritParams flash
#'
#' @param flash A \code{flash} or \code{flash.fit} object.
#'
#' @param kset A vector of integers specifying which factors to nullcheck.
#'   If \code{kset = NULL}, then all existing factors will be checked.
#'
#' @param tol The tolerance parameter: if a factor does not improve the ELBO
#'   by at least \code{tol}, then it will be set to zero.
#'
#' @export
#'
flash.nullcheck <- function(flash,
                            kset = NULL,
                            tol = set.default.tol(flash),
                            verbose.lvl = get.verbose.lvl(flash)) {
  flash <- get.fit(flash)

  if (is.null(kset)) {
    kset <- 1:get.n.factors(flash)
  }
  must.be.valid.kset(flash, kset)

  must.be.numeric(tol, allow.infinite = TRUE, allow.null = FALSE)
  must.be.integer(verbose.lvl, lower = -1, upper = 3)
  must.be.integer(output.lvl, lower = 0, upper = 3)

  announce.nullchk(verbose.lvl, n.factors = length(kset))

  for (k in kset) {
    flash <- nullcheck.factor(flash, k, verbose.lvl, tol)
  }

  if (length(kset) > 0 && !nullchk.failed(flash)) {
    report.nullchk.success(verbose.lvl)
  }

  announce.wrapup(verbose.lvl)
  flash <- wrapup.flash(flash, output.lvl = 3L)

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
